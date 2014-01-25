classdef TrialData
% TrialData represents a collection of trials, whose data is accessed via
% a TrialDataInterface. 
%
% TrialData is not a handle class, meaning that methods which modify this instance
% will return the new TrialData. Any changes made to the underlying TrialDataStore
% will be done via copy-on-write, so that changes will not propagate to another
% TrialData instance.
    
    % Properties stored on disk
    properties
        datasetName = ''; % string describing entire collection of trials dataset

        datasetMeta % arbitrary
    end
    
    properties(SetAccess=protected)
        data = struct();  % standardized format nTrials x 1 struct with all trial data  
        
        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = ''; % how did we access the original data?
        
        manualValid

        timeUnitName
        
        timeUnitsPerSecond
        
        channelDescriptorsByName % struct with ChannelDescriptor for each channel, by name
    end

    properties(Dependent) 
        valid
        nTrials
        nTrialsValid
        nChannels
        channelNames % cell array of channel names
    end
    
    % Initializing and building
    methods
        function td = TrialData(varargin)
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
            end

            if any(~ok)
                msg = msg(~ok);
                fprintf('Error: some data was missing from getChannelData():\n');
                for i = 1:length(msg)
                    fprintf('\t%s\n', msg{i});
                end
                error('Missing data provided by getChannelData');
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

        function printDescriptionShort(td)
            tcprintf('inline', '{yellow}%s: {none}%d trials (%d valid) with %d channels\n', ...
                class(td), td.nTrials, td.nTrialsValid, td.nChannels);
            if ~isempty(td.datasetName)
                tcprintf('inline', '{yellow}Dataset: {none}%s\n', td.datasetName);
            end
        end
        
        function printChannelInfo(td)
            tcprintf('inline', '{bright blue}Analog: {none}%s\n', strjoin(td.listAnalogChannels(), ', '));
            tcprintf('inline', '{bright blue}Event: {none}%s\n', strjoin(td.listEventChannels(), ', '));
            tcprintf('inline', '{bright blue}Param: {none}%s\n', strjoin(td.listParamChannels(), ', '));
            tcprintf('inline', '{bright blue}Spike: {none}%s\n', strjoin(td.listSpikeUnits(), ', '));
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
        function valid = get.valid(td)
            valid = td.buildValid();
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
        
        function vals = replaceInvalidMaskWithValue(td, vals, value)
            if iscell(vals)
                [vals{~td.valid}] = deal(value);
            else
                vals(~td.valid) = value;
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end

    methods % Data access methods
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
        
        function type = getChannelType(td, name)
            type = td.getChannelDescriptor(name).getType();
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
        
        % get the time window for each trial
        function durations = getValidDurations(td)
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
        
        function td = setChannelData(td, name, varargin)
            td.warnIfNoArgOut(nargout);
            fields = td.channelDescriptorsByName.(name).dataFields;
            % for now just set first field
            fld = fields{1};
            vals = varargin{1};
            if ~iscell(vals)
                vals = num2cell(vals);
            end
            
            [td.data.(fld)] = deal(vals{:});
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
            % analog channel
            [~, time] = td.getAnalog(name);
            
            % median of medians is faster and close enough
            delta = nanmedian(cellfun(@(x) nanmedian(diff(x)), time));
        end
        
        function [data, time] = getAnalogRaw(td, name)
            cd = td.channelDescriptorsByName.(name);
            
            data = {td.data.(cd.dataFields{1})}';
            time = {td.data.(cd.dataFields{2})}';
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
        
        function td = setAnalog(td, name, vals)
            td.warnIfNoArgOut(nargout);
            td = td.setChannelData(name, vals);
        end
    end
    
    methods % Event channel methods
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
        function eventStruct = getRawEventStruct(td)
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
        
        function td = setEvent(td, name, vals)
            td.warnIfNoArgOut(nargout);
            td = td.setChannelData(name, vals);
        end
    end
    
    methods % Param channel methods
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
            cd = td.channelDescriptorsByName.(name);
            values = {td.data.(cd.dataFields{1})}';
            % if this channel is marked as a scalar, convert to a numeric array
            if ~cd.collectAsCellByField(1)
                values = cell2mat(values);
            end
        end
        % Basic access methods, very fast
        function values = getParam(td, name)
            cd = td.channelDescriptorsByName.(name);
            values = td.getParamRaw(name);
            values = td.replaceInvalidMaskWithValue(values, cd.missingValueByField{1});
        end

        function paramStruct = getRawChannelDataAsStruct(td, names)
            paramStruct = copyStructField(td.data, [], names);
        end
        
        function paramStruct = getParamStruct(td)
            paramStruct = copyStructField(td.data, [], td.listParamChannels());
        end
        
        function td = setParam(td, name, vals)
            td.warnIfNoArgOut(nargout);
            td = td.setChannelData(name, vals);
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
        
        function tf = hasSpikeUnit(td, unitName)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitName);
            tf = td.hasSpikeChannel(name);
        end
        
        function tf = hasSpikeChannelOrUnit(td, name)
            tf = td.hasSpikeChannel(name);
            if ~tf
                tf = td.hasSpikeUnit(name);
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
        
        function td = selectSpikeUnits(td, names)
            td.warnIfNoArgOut(nargout);
            chNames = cellfun(@SpikeChannelDescriptor.convertUnitNameToChannelName, names, ...
                'UniformOutput', false);  
            td = td.selectSpikeChannels(chNames);
        end
        
        function names = listSpikeUnits(td)
            chNames = td.listSpikeChannels();
            names = cellfun(@SpikeChannelDescriptor.convertChannelNameToUnitName, chNames, ...
                'UniformOutput', false);
        end
        
        function timesCell = getRawSpikeTimes(td, unitName)
            if td.hasSpikeChannel(unitName)
                name = unitName;
            else
                name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitName);
            end
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

    methods % Add data methods
        function td = updatePostDataChange(td)
            td.warnIfNoArgOut(nargout);
        end
        
        % when adding new data to the trial, all times are stored relative
        % to the current time zero. This will be overridden in 
        % TrialDataConditionAlign. This will be added automatically
        % to all new channel time data to match the offsets produced when
        % getting data
        function offsets = getTimeOffsetsFromZeroEachTrial(td)
            offsets = zerosvec(td.nTrials);
        end
        
        function td = addChannel(td, cd, valueCell, varargin)
            p = inputParser();
            p.addParamValue('ignoreOverwriteChannel', false, @islogical);
            p.parse(varargin{:});
            
            td.warnIfNoArgOut(nargout);
            
            % check for overwrite if requested
            assert(isa(cd, 'ChannelDescriptor'), 'Argument cd must be ChannelDescriptor');
            if ~p.Results.ignoreOverwriteChannel && td.hasChannel(cd.name)
                warning('Overwriting existing channel %s', cd.name);
            end
            
            dataFields = cd.dataFields;
            nFields = numel(dataFields);
            
            % check that one value list was provided for each data field
            % referenced by the ChannelDescriptor
            assert(numel(valueCell) == nFields, ...
                'Channel Descriptor references %d fields but only %d field value lists provided', ...
                numel(nFields), numel(valueCell));
            
            for iF = 1:nFields
                % IMPORTANT: DO NOT ACCEPT DATA ON INVALID TRIALS
                if any(~td.valid)
                    warning('New data will be cleared on invalid trials. Recommended to call filterValidTrials() before manipulating data.');
                    valueCell{iF} = td.replaceInvalidMaskWithValue(valueCell{iF}, cd.missingValueByField{iF});
                end
                
                % if valueCell{iF} is a string, copy that field's values
                % (if necessary)
                if ischar(valueCell{iF})
                    if ~strcmp(valueCell{iF}, dataFields{iF})
                        % copy field values since the field names differ
                        td.data = copyStructField(td.data, td.data, valueCell{iF}, dataFields{iF});
                    end
                else
                    % check whether we're overwriting data used by other
                    % channels (besides this one) and throw an error if we are
                    if isfield(td.data, dataFields{iF})
                        otherChannels = setdiff(td.getChannelsReferencingFields(dataFields{iF}), cd.name);
                        if any(otherChannels)
                            error('Refusing to overwrite data field %s referenced by other channels %s', ...
                                strjoin(otherChannels, ', '));
                        end
                    end
                end
                td.data = assignIntoStructArray(td.data, dataFields{iF}, valueCell{iF});
            end 

            td.channelDescriptorsByName.(cd.name) = cd;
            td = td.updatePostDataChange();
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

            assert(~td.hasChannel(name), 'TrialData already has channel with name %s', name);
            assert(numel(values) == td.nTrials, 'Values must be vector with length %d', td.nTrials);

            if isempty(cd)
                cd = ParamChannelDescriptor.buildFromValues(name, values);
            end
            cd.name = name;

            td = td.addChannel(cd, {values});
        end
        
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
        
        function td = addAnalog(td, name, values, times, units)
            td.warnIfNoArgOut(nargout);
            
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
            elseif iscell(times) || isnumeric(times)
                % pass specified times along directly as cell
                % and generate unique time field name
                timeField = genvarname(sprintf('%s_time', name), fieldnames(td.data));
                
            else
                error('Times must be channel/field name or cell array');
            end
                
            % build a channel descriptor for the data
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units);
            
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
            
            % check each trial's time vector against its value vector
            sizeMatch = cellfun(@(time, vals) numel(time) == numel(vals), times, values);
            assert(all(sizeMatch), 'Sizes of times vs. values vectors do not match on %d trials', nnz(sizeMatch));
            
            % add the zero offset to the time vector for each trial
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);
            
            % AnalogChannelDescriptor declares data fields as data, times
            td = td.addChannel(cd, {values, times});
        end
    end 

    methods % Plotting functions
        function units = getChannelUnitsPrimary(td, name)
            if td.hasSpikeChannelOrUnit(name)
                units = 'spikes / sec';
            else
                units = td.getChannelDescriptor(name).unitsPrimary;
            end
        end
        
        function str = getAxisLabelForChannel(td, name)
            str = td.channelDescriptorsByName.(name).getAxisLabelPrimary();
        end

        function str = getTimeAxisLabel(td)
           str = sprintf('Time (%s)',  td.timeUnitName);
        end

        % general utility to send plots to the correct axis
        function axh = getRequestedPlotAxis(td, varargin) %#ok<INUSL>
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

            [dataCell, timeCell] = td.getAnalog(name);     

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

