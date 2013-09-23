classdef TrialData
% TrialData represents a collection of trials, whose data is accessed via
% a TrialDataInterface. 
%
% TrialData is not a handle class, meaning that methods which modify this instance
% will return the new TrialData. Any changes made to the underlying TrialDataStore
% will be done via copy-on-write, so that changes will not propagate to another
% TrialData instance.
    
    % Obtained from TrialDataInterface
    properties% (SetAccess=protected)
        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = '';

        data = struct();  % standardized format nTrials x 1 struct with all trial data  
        
        valid

        datasetName = ''; % string describing entire collection of trials dataset

        datasetMeta
        
        timeUnitName
        
        timeUnitsPerSecond
        
        channelDescriptorsByName % struct with ChannelDescriptor for each channel, by name
    end

    properties(Dependent) 
        nTrials
        nTrialsValid
        nChannels
        channelNames % cell array of channel names
    end
    
    % Initializing and building
    methods
        function td = TrialData(varargin)
            if ~isempty(varargin)
                td = td.initialize(varargin{:});
            end
        end

        % copy everything over from the TrialDataInterface
        function td = initialize(td, varargin)
            p = inputParser();
            p.addRequired('trialDataInterface', @(tdi) isa(tdi, 'TrialDataInterface'));
            p.parse(varargin{:});
            tdi = p.Results.trialDataInterface;
            
            % copy over basic details from the TrialDataInterface
            td.trialDataInterfaceClass = class(tdi);
            td.datasetName = tdi.getDatasetName();
            td.datasetMeta = tdi.getDatasetMeta();
            td.timeUnitName = tdi.getTimeUnitName();
            td.timeUnitsPerSecond = tdi.getTimeUnitsPerSecond();
            
            % request channel descriptors for both special params and regular channels
            specialParams = makecol(tdi.getSpecialParamChannelDescriptors());
            specialNames = {specialParams.name};
            regularChannels = makecol(tdi.getChannelDescriptors());
            regularNames = {regularChannels.name};

            % check for reserved channel names
            overlap = intersect(specialNames, regularNames);
            if ~isempty(overlap)
                error('getChannelDescriptors() returned the following reserved channel names: %s', strjoin(overlap, ', '));
            end

            % combine all channelDescriptors together
            channelDescriptors = [specialParams; regularChannels];

            nTrials = tdi.getTrialCount();
            nChannels = numel(channelDescriptors);

            % request all channel data at once
            data = tdi.getChannelData(channelDescriptors); 

            % loop over channels and verify
            fprintf('Validating channel data...\n');
            ok = falsevec(nChannels);
            msg = cellvec(nChannels);
            for iChannel = 1:nChannels
                chd = channelDescriptors(iChannel); 

                % check to make sure all fields were provided as expected 
                [ok(iChannel), msg{iChannel}] = chd.checkData(data);
            end

            if any(~ok)
                msg = msg(~ok);
                fprintf('Error: some data was missing from getChannelData():\n');
                for i = 1:length(msg)
                    fprintf('\t%s\n', msg{i});
                end
                error('Missing data provided by getChannelData');
            end

            fprintf('Repairing and converting channel data...\n');
            for iChannel = 1:nChannels
                chd = channelDescriptors(iChannel); 
                data = chd.repairData(data);
                %data = chd.convertData(data);
            end

            td.data = data;
        
            % build channelDescriptorByName
            td.channelDescriptorsByName = struct();
            for i = 1:numel(channelDescriptors)
                td.channelDescriptorsByName.(channelDescriptors(i).name) = channelDescriptors(i);
            end
            
            td.valid = truevec(td.nTrials);
            td.initialized = true;
        end

        function printDescriptionShort(td)
            tcprintf('inline', '{yellow}%s: {bright blue}%d trials (%d valid) with %d channels\n', ...
                class(td), td.nTrials, td.nTrialsValid, td.nChannels);
        end
        
        function printChannelInfo(td)
            tcprintf('inline', '{bright cyan}Analog: {none}%s\n', strjoin(td.listAnalogChannels(), ', '));
            tcprintf('inline', '{bright cyan}Event : {none}%s\n', strjoin(td.listEventChannels(), ', '));
            tcprintf('inline', '{bright cyan}Param : {none}%s\n', strjoin(td.listParamChannels(), ', '));
        end

        function disp(td)
            td.printDescriptionShort();
            td.printChannelInfo();
            fprintf('\n');
        end
    end
    
    % Dependent property implementations
    methods % get. accessors for above properties which simply refer to tdi.?
        function nTrials = get.nTrials(td)
            nTrials = numel(td.data);
        end

        function nTrials = get.nTrialsValid(td)
            nTrials = nnz(td.valid);
        end

        function nChannels = get.nChannels(td)
            nChannels = numel(td.channelNames);
        end
        
        function names = get.channelNames(td)
            names = fieldnames(td.channelDescriptorsByName);
        end
    end

    methods % Trial selection
        % All masking or selecting of the .data array must go through this method
        % not much to do here, but be sure to override this in subclasses
        % so that they can update derivative computations
        function td = selectTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            td.data = td.data(mask);
            td.valid = td.valid(mask);
        end
        
        function td = selectValidTrials(td)
             td.warnIfNoArgOut(nargout);
             td = td.selectTrials(td.valid);
        end

        function td = markInvalid(td, mask)
            td.warnIfNoArgOut(nargout);
            td.valid(mask) = false;
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                message = sprintf('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
                expr = sprintf('debug(''%s'')', message);
                evalin('caller', expr); 
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy();
            end
        end
    end

    methods % Data access methods
        function tf = hasChannel(td, name)
            tf = ismember(name, td.channelNames);
        end
        
        function cds = getChannelDescriptorArray(td)
            fields = fieldnames(td.channelDescriptorsByName);
            for iF = 1:length(fields)
                cds(iF) = td.channelDescriptorsByName.(fields{iF});
            end
        end
    end
    
    methods % Analog channel methods
        function names = listAnalogChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'AnalogChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name};
        end
        
        function [data, time] = getAnalog(td, name)
            cd = td.channelDescriptorsByName.(name);
            data = {td.data.(name)}';
            time = {td.data.(cd.timeField)}';
        end
    end
    
    methods % Event channel methods
        function names = listEventChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'EventChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name};
        end

        function timesCell = getEvent(td, name)
            timesCell = {td.data.(name)}';
        end
        
        % used mainly by AlignInfo to make sure it can access unaligned
        % event info
        function eventStruct = getRawEventStruct(td)
            eventStruct = copyStructField(td.data, [], td.listEventChannels());
        end
        
        function timesCell = getEvents(td, nameCell)
            nEvents = numel(nameCell);
            timesCell = cell(td.nTrials, nEvents); 
            for iEv = 1:nEvents
                timesCell(:, iEv) = td.getEvent(nameCell{iEv});
            end
        end

        function counts = getEventCount(td, name)
            counts = cellfun(@numel, td.getEvent(name));
        end

        function tfList = getEventOccurred(td, name)
            tfList = ~cellfun(@isempty, td.getEvent(name));
        end

        function times = getEventNth(td, name, n)
            if strcmp(n, 'end')
                times = td.getEventLast(name);
                return;
            end
            
            times = cellfun(@getNth, td.getEvent(name));
            
            function t = getNth(times)
                if numel(times) >= n
                    t = times(n);
                else
                    t = NaN;
                end
            end
        end

        function times = getEventFirst(td, name)  
            times = getEventNthOccurrence(td, name, 1);
        end

        function times = getEventLast(td, name)  
            times = cellfun(@getLast, td.getEvent(name));

            function t = getLast(times)
                if ~isempty(times)
                    t = times(end);
                else
                    t = NaN;
                end
            end
        end
    end
    
    methods % Param channel methods
        function names = listParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name};
        end
        
        % Basic access methods, very fast
        function values = getParam(td, name)
            values = {td.data.(name)}';
            % if this channel is marked as a scalar, convert to a numeric array
            if ~td.channelDescriptorsByName.(name).collectAsCell
                values = cell2mat(values);
            end
        end
    end

    methods % Spike channel methods
        function names = listSpikeChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'SpikeChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name};
        end

        function timesCell = getSpikeTimesForUnit(td, unitName, varargin) 
            timesCell = td.getRawSpikeTimesForUnit(td, unitRaw);
        end
        
        function timesCell = getRawSpikeTimesForUnit(td, unitName)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitName);
            timesCell = {td.data.(name)}';
        end
    end

    methods % Add data methods
        function td = updatePostDataChange(td)
            td.warnIfNoArgOut(nargout);
        end

        function td = addParam(td, name, values, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('name', @ischar);
            p.addRequired('values', @isvector);
            p.addParamValue('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.parse(name, values, varargin{:});
            name = p.Results.name;
            values = p.Results.values;
            cd = p.Results.channelDescriptor;

            assert(~tf.hasChannel(name), 'TrialData already has channel with name %s', name);
            assert(numel(values) == td.nTrials, 'Values must be vector with length %d', td.nTrials);

            if isempty(cd)
                cd = ParamChannelDescriptor.inferFromValues(values);
            end
            cd.name = name;

            td.channelDescriptors(end+1) = cd;
            td.channelDescriptorsByName.(name) = cd;

            td.data = assignIntoStructArray(td.data, name, values);

            td = td.updatePostDataChange();
        end
    end

    methods % Plotting functions
        function str = getAxisLabelForChannel(td, name)
            str = td.channelDescriptorsByName.(name).getAxisLabel();
        end

        function str = getTimeAxisLabel(td)
           str = sprintf('Time (%s)',  td.timeUnitName);
        end

        % general utility to send plots to the correct axis
        function axh = getRequestedPlotAxis(td, varargin)
            p = inputParser();
            p.addParamValue('axh', [], @(x) isempty(x) || isscalar(x));
            p.addParamValue('cla', false, @islogical); 
            p.KeepUnmatched;
            p.parse(varargin{:});

            % default to gca
            if isempty(p.Results.axh)
                axh = gca;
            else
                axh = p.Results.axh;
            end

            % optionally clear axis
            if p.Results.cla
                cla(axh);
            end
        end

        function plotAnalogEachTrial(td, name, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            [dataCell timeCell] = td.getAnalog(name);     

            dataCell = dataCell(td.valid);
            timeCell = timeCell(td.valid);

            for i = 1:td.nTrialsValid
                plot(axh, double(timeCell{i}), dataCell{i}, '-', 'Color', 0.5*ones(3,1), ...
                    p.Results.plotOptions{:});
                if i == 1, hold(axh, 'on'); end
            end
            box(axh, 'off');
            
            xlabel(td.getTimeAxisLabel());
            ylabel(td.getAxisLabelForChannel(name));
        end
    end

end

