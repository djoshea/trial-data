classdef TrialData
% TrialData represents a collection of trials, whose data is accessed via
% a TrialDataInterface. 
%
% TrialData is not a handle class, meaning that methods which modify this instance
% will return the new TrialData. Any changes made to the underlying TrialDataStore
% will be done via copy-on-write, so that changes will not propagate to another
% TrialData instance.
    
    % Metadata
    properties
        datasetName = ''; % string describing entire collection of trials dataset

        datasetMeta % arbitrary, user determined 
    end
    
    % Internal data storage
    properties(SetAccess=protected)
        data = struct();  % standardized format nTrials x 1 struct with all trial data  
        
        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = ''; % how did we access the original data?
        
        manualValid

        timeUnitName
        
        timeUnitsPerSecond
        
        channelDescriptorsByName % struct with ChannelDescriptor for each channel, by name
    end
    
    properties(Access=protected, Hidden)
        odc % TrialDataOnDemandCache instance
    end

    % Convenience dependent properties
    properties(Dependent) 
        valid
        invalidCause % cellstr of explanations
        nTrials
        nTrialsValid
        nChannels
        channelNames % cell array of channel names
    end
    
    % Initializing and building
    methods
        function td = TrialData(varargin)
            td.odc = TrialDataOnDemandCache();
            
            if ~isempty(varargin)
                if isa(varargin{1}, 'TrialData')
                    td = td.initializeFromTrialData(varargin{1});
                elseif isa(varargin{1}, 'TrialDataInterface')
                    td = td.initializeFromTrialDataInterface(varargin{:});
                else
                    error('Unknown initializer');
                end
            end
        end

        function td = initializeFromTrialData(td, tdOther)
            % this is used by subclasses, so can't copy to output
            meta = ?TrialData;
            props = meta.PropertyList;

            for iProp = 1:length(props)
                prop = props(iProp);
                if prop.Dependent || prop.Constant || prop.Transient
                    continue;
                else
                    name = prop.Name;
                    td.(name) = tdOther.(name);
                end
            end
        end
        
        % copy everything over from the TrialDataInterface
        function td = initializeFromTrialDataInterface(td, varargin)
            td.warnIfNoArgOut(nargout);
            
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

            %nTrials = tdi.getTrialCount();
            nChannels = numel(channelDescriptors);

            % request all channel data at once
            data = tdi.getChannelData(channelDescriptors);  %#ok<PROP>

            % loop over channels and verify
            %fprintf('Validating channel data...\n');
            ok = falsevec(nChannels);
            msg = cellvec(nChannels);
            for iChannel = 1:nChannels
                chd = channelDescriptors(iChannel); 
                
                % check to make sure all fields were provided as expected 
                [ok(iChannel), msg{iChannel}] = chd.checkData(data); %#ok<PROP>
                
                if ~chd.required
                    % fill in missing values for optional channels but
                    % issue a warning
                    chd.addMissingFields(data); %#ok<PROP>
                    fprintf('Warning: No data provided by interface for channel %s\n', chd.describe());
                    ok(iChannel) = false;
                end
            end

            if any(~ok)
                msg = msg(~ok);
                fprintf('Error: data not provided by interface for channels:\n');
                for i = 1:length(msg)
                    fprintf('\t%s\n', msg{i});
                end
                error('Required channel data not provided by getChannelData');
            end

            %fprintf('Repairing and converting channel data...\n');
            for iChannel = 1:nChannels
                chd = channelDescriptors(iChannel); 
                data = chd.repairData(data); %#ok<PROP>
                data = chd.convertDataToMemoryClass(data); %#ok<PROP>
            end

            td.data = data; %#ok<PROP>
        
            % build channelDescriptorByName
            td.channelDescriptorsByName = struct();
            for i = 1:numel(channelDescriptors)
                td.channelDescriptorsByName.(channelDescriptors(i).name) = channelDescriptors(i);
            end
            
            td.manualValid = truevec(td.nTrials);
            td.initialized = true;
        end
    end

    % General utilities
    methods
        function printDescriptionShort(td)
            tcprintf('inline', '{yellow}%s: {none}%d trials (%d valid) with %d channels\n', ...
                class(td), td.nTrials, td.nTrialsValid, td.nChannels);
            if ~isempty(td.datasetName)
                tcprintf('inline', '{yellow}Dataset: {none}%s\n', td.datasetName);
            end
        end
        
        function printChannelInfo(td)
            tcprintf('inline', '{yellow}Analog: {none}%s\n', strjoin(td.listAnalogChannels(), ', '));
            tcprintf('inline', '{yellow}Event: {none}%s\n', strjoin(td.listEventChannels(), ', '));
            tcprintf('inline', '{yellow}Param: {none}%s\n', strjoin(td.listParamChannels(), ', '));
            tcprintf('inline', '{yellow}Spike: {none}%s\n', strjoin(td.listSpikeChannels(), ', '));
        end

        function disp(td)
            td.printDescriptionShort();
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end
    end
    
    % Dependent property implementations
    methods % get. accessors for above properties which simply refer to tdi.?
        function v = get.valid(td)
            if isempty(td.odc), v = td.buildValid(); return; end
            v = td.odc.valid;            
            if isempty(v)
                td.odc.valid = td.buildValid();
                v = td.odc.valid;
            end
        end
        
        function td = set.valid(td, v)
            if isempty(td.odc), td.odc = TrialDataOnDemandCache(); end
            td.odc = td.odc.copy();
            td.odc.valid = v;
        end
        
        function valid = buildValid(td)
            % compute the valid flag considering only trials marked as
            % manually invalid to be invalid. This will be overriden in
            % TDCA to consider the condition and align invalid as well
            valid = td.getManualValid();
        end
        
        function valid = getManualValid(td)
            if isempty(td.manualValid)
                valid = truevec(td.nTrials);
            else
                valid = makecol(td.manualValid);
            end
        end
        
        function cause = get.invalidCause(td)
            cause = td.buildInvalidCause();
        end
        
        function cause = buildInvalidCause(td)
            cause = cell(td.nTrials, 1);
            cause(~td.valid) = {'marked invalid manually'};
            cause(td.valid) = {''};
        end
        
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
            if isempty(td.manualValid)
                td.manualValid = truevec(numel(td.data));
            else
                td.manualValid = td.manualValid(mask);
            end
        end
        
        function td = selectValidTrials(td)
             td.warnIfNoArgOut(nargout);
             td = td.selectTrials(td.valid);
        end

        function td = markInvalid(td, mask)
            td.warnIfNoArgOut(nargout);
            td.manualValid(mask) = false;
            td = td.updateValid();
        end
        
        function td = setInvalid(td, mask)
            td.manualValid = true(td.nTrials, 1);
            td.manualValid(mask) = false;
            td = td.updateValid();
        end
        
        function td = setAllValid(td, mask)
            td.manualValid = true(td.nTrials, 1);
            td = td.updateValid();
        end
        
        function td = updateValid(td)
            td.warnIfNoArgOut(nargout);
            td.valid = [];
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect', ...
                    class(obj));
            end
        end
        
        function vals = replaceInvalidMaskWithValue(td, vals, value)
            if iscell(vals)
                [vals{~td.valid}] = deal(value);
            else
                if isempty(value)
                    value = NaN;
                end
                vals(~td.valid) = value;
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end

    methods % Channel metadata access and manipulation
        function tf = hasChannel(td, name)
            tf = ismember(name, td.channelNames);
        end
        
        function assertHasChannel(td, name)
            assert(td.hasChannel(name), 'TrialData does not have channel %s', name);
        end
        
        function cd = getChannelDescriptor(td, name)
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
        end
        
%         function cd = setChannelDescriptor(td, name, cd)
%             td.assertHasChannel(name);
%             assert(isa(cd, 'ChannelDescriptor'));
%             td.channelDescriptorsByName.(name) = cd;
%         end
        
        function type = getChannelType(td, name)
            type = td.getChannelDescriptor(name).getType();
        end
        
        function units = getChannelUnitsPrimary(td, name)
            % return a string describing the units of a given channel
            if td.hasSpikeChannel(name)
                units = 'spikes / sec';
            else
                units = td.getChannelDescriptor(name).unitsPrimary;
            end
        end
        
        function td = setChannelUnitsPrimary(td, name, units)
            td.warnIfNoArgOut(nargout);
            td.channelDescriptorsByName.(name) = td.getChannelDescriptor(name).setPrimaryUnits(units);
        end
        
        function names = listChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            names = {channelDescriptors.name}';
        end
        
        function names = listSpecialChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) cd.special, channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end

        function names = listNonSpecialChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) ~cd.special, channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function durations = getValidDurations(td)
            % get the time window for each trial
            starts = td.getEventFirst('TrialStart');
            ends = td.getEventLast('TrialEnd');
            durations = ends - starts;
            durations = td.replaceInvalidMaskWithValue(durations, NaN);
        end
        
        function cds = getChannelDescriptorArray(td)
            fields = fieldnames(td.channelDescriptorsByName);
            for iF = 1:length(fields)
                cds(iF) = td.channelDescriptorsByName.(fields{iF}); %#ok<AGROW>
            end
        end
        
        % channels may reference common data fields in .data to prevent
        % data duplication. This function returns the list of channels
        % which reference a particular data field. fields is a cellstr
        % and nameCell is a cell of cellstr
        function namesByField = getChannelsReferencingFields(td, fields)
            if ischar(fields)
                wasCell = false;
                fields = {fields};
            else
                wasCell = true;
            end
            channels = fieldnames(td.channelDescriptorsByName);
            fieldsByChannel = structfun(@(cd) cd.dataFields, ...
                td.channelDescriptorsByName, 'UniformOutput', false);
            
            namesByField = cellvec(numel(fields));
            for iF = 1:numel(fields)
                assert(isfield(td.data, fields{iF}), 'TrialData does not have data field %s', fields{iF});
                mask = structfun(@(chFields) ismember(fields{iF}, chFields), fieldsByChannel);
                namesByField{iF} = channels(mask);
            end
            if ~wasCell
                namesByField = namesByField{1};
            end
        end
        
        % drops all channel data except specified and special
        % parameter channels
        function td = dropNonSpecialChannelsExcept(td, names)
            td.warnIfNoArgOut(nargout);
            removeNames = setdiff(td.listNonSpecialChannels(), names);
            td = td.dropChannels(removeNames);
        end
        
        % drop all channels except specified and 
        function td = dropChannelsExcept(td, names)
            td.warnIfNoArgOut(nargout);
            removeNames = setdiff(td.listChannels(), names);
            td = td.dropChannels(removeNames);
        end
        
        function td = dropChannels(td, names)
            td.warnIfNoArgOut(nargout);
            
            if isempty(names)
                return;
            end
            
            if ischar(names)
                names = {names};
            end
            
            % first hold onto the to-be-removed channel descriptors
            cds = cellvec(numel(names));
            for i = 1:numel(names)
                cds{i} = td.channelDescriptorsByName.(names{i});
            end
            
            % remove the channel descriptors
            td.channelDescriptorsByName = rmfield(td.channelDescriptorsByName, names);
            
            % for the removed data channels' fields, figure out which ones 
            % aren't referenced by any other channels
            fieldsRemoveByChannel = cellfun(@(cd) makecol(cd.dataFields), ...
                cds, 'UniformOutput', false);
            fieldsRemove = unique(cat(1, fieldsRemoveByChannel{:}));
            
            otherChannelsReferencingFields = td.getChannelsReferencingFields(fieldsRemove);
            maskRemove = cellfun(@isempty, otherChannelsReferencingFields);
            fieldsRemove = fieldsRemove(maskRemove);
            
            td.data = rmfield(td.data, fieldsRemove);
            td = td.updatePostDataChange();
        end
    end
    
    methods % Analog channel methods
        function td = addAnalog(td, name, varargin)
            p = inputParser();
            p.addOptional('values', {}, @(x) iscell(x) || ismatrix(x));
            p.addOptional('times', {}, @(x) ischar(x) || iscell(x) || isvector(x));
            p.addParamValue('units', '', @ischar);
            p.parse(varargin{:});
            times = p.Results.times;
            values = p.Results.values;
            units = p.Results.units;

            td.warnIfNoArgOut(nargout);
            
            % times can either be a field/channel name, or it can be raw
            % time values
            if ischar(times)
                if td.hasChannel(times)
                    % treat times as analog channel name
                    % share that existing channel's time field 
                    cd = td.channelDescriptorsByName.(times);
                    assert(isa(cd, 'AnalogChannelDescriptor'), ...
                        'Channel %s is not an analog channel', times);
                    timeField = cd.timeField;
                    times = {td.data.(timeField)};
                    
                elseif isfield(td.data, times)
                    % use directly specified time field in .data
                    timeField = times;
                    times = {td.data.(timeField)};
                    
                else
                    error('%s is not a channel or data field name', times);
                end
                times = [];

            elseif iscell(times) || isnumeric(times)
                % pass specified times along directly as cell
                % and generate unique time field name
                timeField = genvarname(sprintf('%s_time', name), fieldnames(td.data));
            else
                error('Times must be channel/field name or cell array');
            end
                
            % build a channel descriptor for the data
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units);
            
            % AnalogChannelDescriptor declares data fields as data, times
            td = td.addChannel(cd);
            
            if ~isempty(values)
                td = td.setAnalog(name, values, times, ...
                    'clearForInvalid', false);
            end
        end
        
        function td = setAnalog(td, name, values, varargin)
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);
            
            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) ||  isnumeric(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            times = p.Results.times;
            
            % check the values and convert to nTrials cellvec
            if ismatrix(values) && isnumeric(values)
                % values must be nTrials x nTimes
                assert(size(values, 1) == td.nTrials, 'Values as matrix must be nTrials along dimension 1');
                values = mat2cell(values', size(values, 2), onesvec(td.nTrials))';
            elseif iscell(values)
                assert(numel(values) == td.nTrials, 'Values as cell must have numel == nTrials');
            else
                error('Values must be numeric matrix or cell array');
            end
                
            if ~isempty(times)
                if iscell(times)
                    assert(numel(times) == td.nTrials, 'numel(times) must match nTrials'); 
                elseif ismatrix(times)
                    if isvector(times)
                        times = repmat(makerow(times), td.nTrials, 1);
                    end
                    assert(size(times,1) == td.nTrials, 'size(times, 1) must match nTrials');  
                    times = mat2cell(times', size(times, 2), onesvec(td.nTrials))';

                else
                    error('Times must be numeric matrix or cell array');
                end
            else
                % pass along the current times since the data is coming in with the
                % existing alignment
                [~, times] = td.getAnalog(name);
            end
            
            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);

            td = td.setChannelData(name, {values, times}, p.Unmatched);
        end

        function tf = hasAnalogChannel(td, name) 
            if td.hasChannel(name)
                tf = isa(td.getChannelDescriptor(name), 'AnalogChannelDescriptor');
            else
                tf = false;
            end
        end
        
        function names = listAnalogChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'AnalogChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function td = selectAnalogChannels(td, names)
            td.warnIfNoArgOut(nargout);
            full = td.listAnalogChannels();
            assert(all(ismember(names, full)), 'Missing analog channels %s', ...
                strjoin(setdiff(names, full), ', ')); 
            td = td.dropChannels(setdiff(full, names));
        end

        function delta = getAnalogTimeDelta(td, name)
            % compute the median delta betwen successive samples of an
            % analog channel(s), returns the minimum timeDelta across all channels
            
            if ischar(name)
                name = {name};
            end

            delta = nanvec(numel(name));
            for i = 1:numel(name)
                [~, time] = td.getAnalog(name{i});
                % median of medians is faster and close enough
                delta(i) = nanmedian(cellfun(@(x) nanmedian(diff(x)), time));
            end

            delta = nanmin(delta);
        end
        
        function [data, time] = getAnalogRaw(td, name)
            cd = td.channelDescriptorsByName.(name);
            
            data = {td.data.(cd.dataFields{1})}';
            time = {td.data.(cd.dataFields{2})}';
            for i = 1:numel(data)
                if numel(data{i}) == numel(time{i}) - 1
                    time{i} = makecol(time{i}(1:end-1));
                else
                    time{i} = makecol(time{i});
                end
                data{i} = makecol(data{i});
            end
        end
        
        % same as raw, except empty out invalid trials
        function [data, time] = getAnalog(td, name)
            [data, time] = td.getAnalogRaw(name);
            data = td.replaceInvalidMaskWithValue(data, []);
            time = td.replaceInvalidMaskWithValue(time, []);
        end
        
        function [dataVec, timeVec] = getAnalogSample(td, name, varargin)
            [dataCell, timeCell] = td.getAnalog(name);
            dataVec = cellfun(@(v) v(1), dataCell, 'ErrorHandler', @(varargin) NaN);
            if nargout > 1
                timeVec = cellfun(@(v) v(1), timeCell, 'ErrorHandler', @(varargin) NaN);
            end
        end 
    end
    
    methods % Event channel methods
        function td = addEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addRequired('name', @ischar);
            p.addRequired('times', @isvector);
            %p.addParamValue('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.parse(name, times, varargin{:});
            %cd = p.Results.channelDescriptor;
            
            assert(~td.hasChannel(name), 'TrialData already has channel with name %s', name);
            assert(numel(times) == td.nTrials, 'Times must be vector with length %d', td.nTrials);
            times = makecol(times);
                
            if iscell(times)
                % multiple occurrence event
                cd = EventChannelDescriptor.buildMultipleEvent(name, td.timeUnitName);
                
                % check contents
                assert(all(cellfun(@(x) isvector(x) && isnumeric(x), times)), ...
                    'Cell elements must be numeric vector of event occurrence times');
            elseif isnumeric(times)
                % single occurrence event
                cd = EventChannelDescriptor.buildSingleEvent(name, td.timeUnitName);
            else
                error('Times must be numeric vector or cell vector of numeric vectors');
            end
         
            % for TDCA, assume events come in aligned to the current 'zero' time
            % we add the zero offset to the times so that they are stored
            % as absolute time points
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            
            if iscell(times)
                times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);
            else
                times = times + offsets;
            end
            
            td = td.addChannel(cd, {times});
        end
        
        function tf = hasEventChannel(td, name) 
            if td.hasChannel(name)
                tf = isa(td.getChannelDescriptor(name), 'EventChannelDescriptor');
            else
                tf = false;
            end
        end
        
        function names = listEventChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'EventChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        % used mainly by AlignInfo to make sure it can access unaligned
        % event info
        function eventStruct = getRawEventFlatStruct(td, chList)
            if nargin < 2
                chList = td.listEventChannels();
            end
            
            for iCh = 1:numel(chList)
                ch = chList{iCh};
                if td.channelDescriptorsByName.(ch).isScalarByField(1)
                    eventStruct.(ch) = makecol([td.data.(ch)]);
                else
                    eventStruct.(ch) = makecol({td.data.(ch)});
                end
            end
        end 
            
        function eventStruct = getRawEventStructArray(td)
            eventStruct = copyStructField(td.data, [], td.listEventChannels());
        end
                
        function timesCell = getEventRaw(td, name)
            timesCell = {td.data.(name)}';
        end
        
        % replace invalid trials with []
        function timesCell = getEvent(td, name)
            timesCell = td.getEventRaw(name);
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
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
            times = getEventNth(td, name, 1);
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
        
        function td = setEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);
            
            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            if iscell(times)
                times = cellfun(@plus, makecol(times), num2cell(offsets), 'UniformOutput', false);
            else
                times = makecol(times) + offsets;
            end
            td = td.setChannelData(name, {times}, varargin{:});
            
        end
    end
    
    methods % Param channel methods
        function td = addParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('name', @ischar);
            p.addOptional('values', @isvector);
            p.addParamValue('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.parse(name, varargin{:});
            
            name = p.Results.name;
            values = p.Results.values;
            cd = p.Results.channelDescriptor;

%             if td.hasChannel(name)
%                 warning('Overwriting existing param channel with name %s', name);
%             end
            if ~isempty(values)
                assert(numel(values) == td.nTrials, 'Values must be vector with length %d', td.nTrials);
            end
            
            if isempty(cd)
                cd = ParamChannelDescriptor.buildFromValues(name, values);
            end

            td = td.addChannel(cd, {values});
        end
        
        function td = addScalarParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.addParamValue('units', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cd = ParamChannelDescriptor.buildScalarParam(name, p.Results.units);
            values = p.Results.values;
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end
        
        function td = addStringParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cd = ParamChannelDescriptor.buildStringParam(name);
            values = p.Results.values;
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end
        
        function td = addBooleanParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cd = ParamChannelDescriptor.buildBooleanParam(name);
            values = p.Results.values;
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end
        
        function tf = hasParamChannel(td, name) 
            if td.hasChannel(name)
                tf = isa(td.getChannelDescriptor(name), 'ParamChannelDescriptor');
            else
                tf = false;
            end
        end
        
        function names = listParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function names = listScalarParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray(); 
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor') && cd.isScalar, ...
                channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end

        function names = getStringParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray(); 
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor') && cd.isString, ...
                channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function values = getParamRaw(td, name)
            % grab the raw value for a parameter without considering
            % validity. if name is an analog channel, grab the first sample
            % of that channel
            cd = td.channelDescriptorsByName.(name);
            
            if isa(cd, 'ParamChannelDescriptor')
                values = {td.data.(cd.dataFields{1})}';
                % if this channel is marked as a scalar, convert to a numeric array
                if ~cd.collectAsCellByField(1)
                    values = cellfun(@double, values);
                end
            elseif isa(cd, 'AnalogChannelDescriptor')
                values = td.getAnalogSample(name);
            else
                error('Only valid for parameter or analog channels');
            end
        end
        % Basic access methods, very fast
        function values = getParam(td, name)
            cd = td.channelDescriptorsByName.(name);
            values = td.getParamRaw(name);
            values = td.replaceInvalidMaskWithValue(values, cd.missingValueByField{1});
        end
        
        function values = getParamUnique(td, name)
            vals = td.getParam(name);
            if ~iscell(vals)
                vals =removenan(vals);
            end
            values = unique(vals);
        end

        function paramStruct = getRawChannelDataAsStruct(td, names)
            paramStruct = copyStructField(td.data, [], names);
        end
        
        function paramStruct = getParamStruct(td)
            paramStruct = copyStructField(td.data, [], td.listParamChannels());
        end
        
        function td = setParam(td, name, vals, varargin)
            td.warnIfNoArgOut(nargout);
            cd = td.channelDescriptorsByName.(name);
            % convert to strcell
            if cd.isStringByField(1) && ischar(vals)
                vals = {vals};
            end
            % clone scalar to size
            if isscalar(vals)
                vals = repmat(vals, td.nTrials, 1);
            end  
            td = td.setChannelData(name, {vals}, varargin{:});
        end
    end

    methods % Spike channel methods
        function tf = hasSpikeChannel(td, name)
            if td.hasChannel(name)
                tf = isa(td.getChannelDescriptor(name), 'SpikeChannelDescriptor');
            else
                tf = false;
            end
        end
        
        function names = listSpikeChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'SpikeChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function td = selectSpikeChannels(td, names)
            td.warnIfNoArgOut(nargout);
            full = td.listSpikeChannels();
            assert(all(ismember(names, full)), 'Missing spike channels %s', ...
                strjoin(setdiff(names, full), ', ')); 
            td = td.dropChannels(setdiff(full, names));
        end
        
        function timesCell = getRawSpikeTimes(td, unitName)
%             if td.hasSpikeChannel(unitName)
                name = unitName;
%             else
%                 name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitName);
%             end
            timesCell = {td.data.(name)}';
        end

        function timesCell = getSpikeTimes(td, unitName, varargin) 
            timesCell = td.getRawSpikeTimes(unitName);
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
        end
            
        function counts = getSpikeCounts(td, unitName)
            counts = cellfun(@numel, td.getSpikeTimes(unitName));
        end
        
        function rates = getSpikeMeanRate(td, unitName)
            counts = td.getSpikeCounts(unitName);
            durations = td.getValidDurations();
            rates = counts ./ durations * td.timeUnitsPerSecond;
        end
    end

    methods % Generic add data methods
        function td = updatePostDataChange(td)
            td.warnIfNoArgOut(nargout);
        end
        
        function offsets = getTimeOffsetsFromZeroEachTrial(td)
            % when adding new data to the trial, all times are stored relative
            % to the current time zero. This will be overridden in 
            % TrialDataConditionAlign. This will be added automatically
            % to all new channel time data to match the offsets produced when
            % getting data
            offsets = zerosvec(td.nTrials);
        end
        
        function td = addChannel(td, cd, varargin)
            % adds a new channel described by ChannelDescriptor cd.
            % valueCell must be nDataFields x 1 cell each with nTrials x 1
            % cell within. Note that if any trials are marked invalid, this
            % data will be cleared when adding the channel.
            p = inputParser();
            p.addOptional('valueCell', {}, @(x) true);
            p.addParamValue('ignoreOverwriteChannel', false, @islogical);
            p.parse(varargin{:});
            valueCell = p.Results.valueCell;
            
            td.warnIfNoArgOut(nargout);
            
            % check for overwrite if requested
            assert(isa(cd, 'ChannelDescriptor'), 'Argument cd must be ChannelDescriptor');
%             if ~p.Results.ignoreOverwriteChannel && td.hasChannel(cd.name)
%                 warning('Overwriting existing channel %s', cd.name);
%             end

            td.channelDescriptorsByName.(cd.name) = cd;
            
            % touch each of the value fields to make sure they exist
            for iF = 1:cd.nFields
                if ~isfield(td.data, cd.dataFields{iF});
                    td.data(end).(cd.dataFields{iF}) = [];
                end
            end     
            
            if isempty(valueCell)
                td = td.clearChannelData(cd.name, 'fieldMask', ~cd.isShareableByField);
            else
                % clear on fields where no values provided and it's not shared, 
                % set on fields where values are provided
                nonEmptyMask = ~cellfun(@isempty, valueCell);
                td = td.clearChannelData(cd.name, 'fieldMask', ~nonEmptyMask && ~cd.isShareableByField);
                td = td.setChannelData(cd.name, valueCell, 'fieldMask', nonEmptyMask);
            end
        end
        
        function td = copyRenameSharedChannelFields(td, name, fieldMask)
            % copy all fields which belong to this channel which are shared
            % with another channel, rename those fields, and mark the
            % updated fields inside that fields channelDescriptor
           td.warnIfNoArgOut(nargout); 
           
           cd = td.channelDescriptorsByName.(name);
           if nargin < 3 || isempty(fieldMask)
               fieldMask = true(cd.nFields, 1);
           end
           isShareable = cd.isShareableByField; 
           
            dataFields = cd.dataFields;
            for iF = 1:cd.nFields
                if ~isShareable(iF) || ~fieldMask(iF), continue; end
                
                otherChannels = setdiff(td.getChannelsReferencingFields(dataFields{iF}), cd.name);
                if isempty(otherChannels), continue; end
                
                % being used by other channels, rename and copy
                newName = genvarname([dataFields{iF} '_' name 'Copy'], fieldnames(td.data));
                td.data = copyStructField(td.data, td.data, dataFields{iF}, newName);
                
                cd = cd.renameDataField(iF, newName);
            end
            
            td.channelDescriptorsByName.(name) = cd;
        end
        
        function td = clearChannelData(td, name, varargin)
            p = inputParser();
            p.addParamValue('fieldMask', [], @islogical);
            p.parse(varargin{:});
            
            cd = td.channelDescriptorsByName.(name);
            fieldMask = p.Results.fieldMask;
            if isempty(fieldMask)
                fieldMask = true(cd.nFields, 1);
            end
            
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);
            
            % avoid clearing data that other channels reference 
            td = td.copyRenameSharedChannelFields(name, p.Results.fieldMask);

            % get the channel descriptor again after performing any renames above
            cd = td.channelDescriptorsByName.(name);
            
            % update the values of a given channel's data fields on 
            dataFields = cd.dataFields;
            nFields = numel(dataFields);
            for iF = 1:nFields
                if ~fieldMask(iF), continue; end
                val = cd.missingValueByField{iF};
                for iT = 1:td.nTrials
                    td.data(iT).(dataFields{iF}) = val;
                end
            end
        end
        
        function td = setChannelData(td, name, valueCell, varargin)
            % note that by default, updateValidOnly is true, meaning that
            % the values on invalid trials will not be updated
            p = inputParser();
            p.addParamValue('fieldMask', [], @islogical);
            p.addParamValue('clearForInvalid', false, @islogical);
            p.addParamValue('updateValidOnly', true, @islogical);
            p.addParamValue('updateMask', [], @isvector);
            p.parse(varargin{:});
            
            updateMaskManual = p.Results.updateMask;
            if isempty(updateMaskManual)
                updateMaskManual = true(td.nTrials, 1);
            end
            assert(numel(updateMaskManual) == td.nTrials, ...
                'Size of updateMask must match nTrials');
            
            % update the values of a given channel's data fields on 
            td.warnIfNoArgOut(nargout);
            
            cd = td.channelDescriptorsByName.(name);
            fieldMask = p.Results.fieldMask;
            if isempty(fieldMask)
                fieldMask = true(cd.nFields, 1);
            else
                assert(numel(fieldMask) == cd.nFields);
            end
            
            % avoid overwriting data shared by other channels
            td = td.copyRenameSharedChannelFields(cd.name, fieldMask);
            cd = td.channelDescriptorsByName.(name);
            
            dataFields = cd.dataFields;
            nFields = numel(dataFields);
            
            % check that one value list was provided for each data field
            % referenced by the ChannelDescriptor
            assert(numel(valueCell) == nFields, ...
                'Channel Descriptor references %d fields but only %d field value lists provided', ...
                numel(nFields), numel(valueCell));
            
            for iF = 1:nFields
                % only touch specified fields
                if ~fieldMask(iF), continue; end
                
                % IMPORTANT: DO NOT ACCEPT DATA ON INVALID TRIALS
                updateMask = updateMaskManual;
                if p.Results.updateValidOnly
                    updateMask = updateMask & td.valid;
                elseif p.Results.clearForInvalid
                    valueCell{iF} = td.replaceInvalidMaskWithValue(valueCell{iF}, cd.missingValueByField{iF});
                end
               
                td.data = assignIntoStructArray(td.data, dataFields{iF}, valueCell{iF}, updateMask);
            end 

            td.channelDescriptorsByName.(cd.name) = cd;
            td = td.updatePostDataChange();
        end
    end 

    methods % Plotting functions
        function str = getAxisLabelForChannel(td, name)
            str = td.channelDescriptorsByName.(name).getAxisLabelPrimary();
        end

        function str = getTimeAxisLabel(td)
           str = sprintf('Time (%s)',  td.timeUnitName);
        end

        % general utility to send plots to the correct axis
        function [axh, unmatched] = getRequestedPlotAxis(td, varargin) %#ok<INUSL>
            p = inputParser();
            p.addParamValue('axh', [], @(x) isempty(x) || isscalar(x));
            p.addParamValue('cla', false, @islogical); 
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % default to gca
            if isempty(p.Results.axh)
                axh = newplot();
            else
                axh = p.Results.axh;
                if ~ishold(axh)
                    cla(axh);
                end
            end

            % optionally clear axis
            if p.Results.cla
                cla(axh);
            end
            
            unmatched = p.Unmatched;
        end

        function [pan, unmatched] = getRequestedPlotPanel(td, varargin) %#ok<INUSL>
            p = inputParser();
            p.addParamValue('figh', [], @(x) isempty(x) || ishandle(x));
            p.addParamValue('panel', [], @(x) isempty(x) || ishandle(x));
            p.addParamValue('clf', true, @islogical); 
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % default to gcf
            if isempty(p.Results.figh)
                figh = gcf;
            else
                figh = p.Results.figh;
            end
           
            % optionally clear figure
            if p.Results.clf
                clf(figh);
            end
            
            pan = OuterPanel(figh);
            
            unmatched = p.Unmatched;
        end
        
        function plotAnalogEachTrial(td, name, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            [dataCell, timeCell] = td.getAnalog(name);     

            dataCell = dataCell(td.valid);
            timeCell = timeCell(td.valid);

            cla(axh);
            hold(axh, 'on');
            
            for i = 1:td.nTrialsValid
                if isempty(dataCell{i}), continue, end;
                
                % allow for derivative channels to have one fewer time
                % point
                if numel(dataCell{i}) == numel(timeCell{i}) - 1
                    timeCell{i} = timeCell{i}(1:end-1);
                end
                plot(axh, double(timeCell{i}), dataCell{i}, '-', 'Color', 0.5*ones(3,1), ...
                    p.Results.plotOptions{:});
                
            end
            box(axh, 'off');
            
            xlabel(td.getTimeAxisLabel());
            ylabel(td.getAxisLabelForChannel(name));
        end
    end
end

