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

        datasetName = ''; % string describing entire collection of trials dataset

        datasetMeta
        
        timeUnitName
        
        timeUnitsPerSecond
        
        channelDescriptors % struct with ChannelDescriptor for each channel, by name
        
        analogChannelMask
        eventChannelMask
        paramChannelMask
        
        analogChannelNames
        eventChannelNames
        paramChannelNames
        
        channelDescriptorsByName
        
        channelNames % cell array of channel names
    end

    properties(Dependent) 
        nTrials
        nChannels
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
            p.addOptional('trialDataInterface', [], @(tdi) isa(tdi, 'TrialDataInterface'));
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
            td.channelDescriptors = [specialParams; regularChannels];

            nTrials = tdi.getTrialCount();
            nChannels = numel(td.channelDescriptors);

            % request all channel data at once
            td.channelNames = {td.channelDescriptors.name};
            data = tdi.getChannelData(td.channelDescriptors); 

            % loop over channels and verify
            for iChannel = 1:nChannels
                chd = td.channelDescriptors(iChannel); 
                name = chd.name;

                % check for main field
                assert(isfield(data, name), 'getChannelData missing field %s', name);

                % replace empty values with default
                if ~isempty(chd.defaultValue)
                    emptyMask = cellfun(@isempty, {data.(name)});
                    [data(emptyMask).(name)] = deal(chd.defaultValue);
                 end

                % convert to appropriate type
                data = structConvertFieldValues(data, chd.dataClass, chd.name);
                
                % these are the data subfields for this channel
                chExtraFields = chd.getExtraDataFields();
                for iEx = 1:numel(chExtraFields)
                    subfield = sprintf('%s_%s', name, chExtraFields{iEx});
                    assert(isfield(data, subfield), 'getChannelData missing field %s', subfield);
                end
            end

            td.data = data;
            
            % build masks for each channel type
            td.analogChannelMask = arrayfun(@(cd) isa(cd, 'AnalogChannelDescriptor'), td.channelDescriptors);
            td.eventChannelMask = arrayfun(@(cd) isa(cd, 'EventChannelDescriptor'), td.channelDescriptors);
            td.paramChannelMask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor'), td.channelDescriptors);
            
            % collect names of channels by type
            td.analogChannelNames = {td.channelDescriptors(td.analogChannelMask).name};
            td.eventChannelNames = {td.channelDescriptors(td.eventChannelMask).name};
            td.paramChannelNames = {td.channelDescriptors(td.paramChannelMask).name};
        
            % build channelDescriptorByName
            td.channelDescriptorsByName = struct();
            for i = 1:numel(td.channelDescriptors)
                td.channelDescriptorsByName.(td.channelDescriptors(i).name) = td.channelDescriptors(i);
            end
            
            td.initialized = true;
        end
    end
    
    % Dependent property implementations
    methods % get. accessors for above properties which simply refer to tdi.?
        function nTrials = get.nTrials(td)
            nTrials = numel(td.data);
        end

        function nChannels = get.nChannels(td)
            nChannels = numel(td.channelDescriptors);
        end
    end

    methods % Trial selection
        % All masking or selecting of the .data array must go through this method
        % not much to do here, but be sure to override this in subclasses
        % so that they can update derivative computations
        function td = selectTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            td.data = td.data(mask);
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
            
        % LIST CHANNEL NAMES BY TYPE
        function names = listAnalog(td)
            names = td.analogChannelNames;
        end

        function names = listEvents(td)
            names = td.eventChannelNames;
        end

        function names = listParams(td)
            names = td.paramChannelNames;
        end

        % Basic access methods, very fast
        function values = getParam(td, name)
            values = {td.data.(name)}';
            % if this channel is marked as a scalar, convert to a numeric array
            if td.channelDescriptorsByName.(name).scalar
                values = cell2mat(values);
            end
        end
        
        function [data time] = getAnalog(td, name)
            data = {td.data.(name)}';
            time = {td.data.([name '_time'])}';
        end
        
        function [timesCell tags] = getEvent(td, name)
            timesCell = {td.data.(name)}';
            tags = {td.data.([name '_tags'])}';
        end
        
        function timesCell = getEvents(td, nameCell)
            nEvents = numel(nameCell);
            timesCell = cell(td.nTrials, nEvents); 
            for iEv = 1:nEvents
                timesCell(:, iEv) = td.getEvent(nameCell{iEv});
            end
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

    methods % Special event access
        function counts = getEventCount(td, name)
            counts = cellfun(@numel, td.getEvent(name));
        end

        function tfList = getEventOccurred(td, name)
            tfList = ~cellfun(@isempty, td.getEvent(name));
        end

        function [times] = getEventNth(td, name, n)
            timesCell = td.getEvent(name);
            
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
            timesCell = td.getEvent(name);
            
            function t = getLast(times)
                if ~isempty(times)
                    t = times(end);
                else
                    t = NaN;
                end
            end
        end
    end
    

end

