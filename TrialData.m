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
    end

    % Internal data storage
    properties(SetAccess=protected)
        trialDataVersion = 20160313; % used for backwards compatibility
        
        datasetMeta = struct(); % arbitrary, user determined data, use setMetaByKey and getMetaByKey to access

        data = struct([]);  % standardized format nTrials x 1 struct with all trial data  

        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = ''; % how did we access the original data?

        manualValid = true(0, 1);
        manualInvalidCause = {};
        
        temporaryValid = [];
        temporaryInvalidCause = {};

        timeUnitName = 'ms';

        timeUnitsPerSecond = 1000;

        channelDescriptorsByName  = struct(); % struct with ChannelDescriptor for each channel, by name
        
        trialDescriptionExtraParams = {} % used in generating descriptions of each trial
    end

    properties(Access=protected, Hidden)
        odc % TrialDataOnDemandCache instance
    end
    
    % properties which are stored in odc, see get/set below
    properties(Dependent)
        valid
        invalidCause % cellstr of explanations
    end

    % Convenience dependent properties
    properties(Dependent)
        ch % a struct which contains all channel names as fields, whose values are the channel names. This is a hack to facilitate tab completion instead of strings.
        nTrials
        nTrialsValid
        
        nTrialsPermanentlyInvalid
        
        nChannels
        channelNames % cell array of channel names
    end

    % Initializing and building
    methods
        function td = TrialData(varargin)
            p = inputParser();
            p.addOptional('buildFrom', [], @(x) isa(x, 'TrialData') || isa(x, 'TrialDataInterface'));
%             p.addParamValue('quiet', false, @islogical); % completely silent output
            p.addParamValue('suppressWarnings', false, @islogical); % don't warn about any minor issues
            p.parse(varargin{:});
%             quiet = p.Results.quiet;
            suppressWarnings = p.Results.suppressWarnings;
            
            td = td.rebuildOnDemandCache();

            if ~isempty(varargin)
                if isa(p.Results.buildFrom, 'TrialData')
                    td = td.initializeFromTrialData(p.Results.buildFrom, ...
                        'suppressWarnings', suppressWarnings);
                elseif isa(p.Results.buildFrom, 'TrialDataInterface')
                    td = td.initializeFromTrialDataInterface(p.Results.buildFrom, ...
                        'suppressWarnings', suppressWarnings);
                else
                    error('Unknown initializer');
                end
            end
        end
        
        function td = rebuildOnDemandCache(td)
            td.warnIfNoArgOut(nargout);
            td.odc = TrialDataOnDemandCache(); 
        end

        function td = initializeFromTrialData(td, tdOther, varargin)
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
            p.addParamValue('suppressWarnings', false, @islogical);
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
            regularChannels = makecol(tdi.getChannelDescriptors('suppressWarnings', p.Results.suppressWarnings));
            regularNames = {regularChannels.name};

            % check for reserved channel names
            overlap = intersect(specialNames, regularNames);
            if ~isempty(overlap)
                error('getChannelDescriptors() returned the following reserved channel names: %s', strjoin(overlap, ', '));
            end

            % combine all channelDescriptors together
            channelDescriptors = [specialParams; regularChannels];

            %nTrials = tdi.getTrialCount();
            %nChannels = numel(channelDescriptors);

            % request all channel data at once
            data = tdi.getChannelData(channelDescriptors);

            % build channelDescriptorByName
            td.channelDescriptorsByName = struct();
            for i = 1:numel(channelDescriptors)
                td.channelDescriptorsByName.(channelDescriptors(i).name) = channelDescriptors(i);
            end
            
            % validate and replace missing values, update data classes
            % inside channelDescriptors
            [data, td.channelDescriptorsByName] = td.validateData(data, td.channelDescriptorsByName, 'suppressWarnings', p.Results.suppressWarnings);
            td.data = makecol(data);

            td.manualValid = truevec(td.nTrials);
            td.initialized = true;
        end

        function td = addNewTrialsFromTrialDataInterface(td, varargin)
            % this is a reduced form of initializeFromTDI that requests only new channel descriptors and new data
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addRequired('trialDataInterface', @(tdi) isa(tdi, 'TrialDataInterface'));
            p.addParamValue('suppressWarnings', false, @islogical);
            p.parse(varargin{:});
            tdi = p.Results.trialDataInterface;

            % check equivalence of basic details from the TrialDataInterface
            assert(strcmp(td.trialDataInterfaceClass, class(tdi)));
            assert(strcmp(td.timeUnitName, tdi.getTimeUnitName()));
            assert(td.timeUnitsPerSecond == tdi.getTimeUnitsPerSecond());

            % request channel descriptors for both special params and regular channels
            newChannels = makecol(tdi.getNewChannelDescriptors());
            maskRemove = ismember({newChannels.name}, td.listChannels());
            newChannels = newChannels(~maskRemove);
            td = td.addChannels(newChannels);

            % request all channel data at once
            newData = tdi.getNewChannelData(td.getChannelDescriptorArray());

            % defer to addNewTrials
            td = td.addNewTrialsRaw(newData, 'suppressWarnings', p.Results.suppressWarnings);
        end

        function td = addTrialsFromTrialData(td, trialData, varargin)
            td.warnIfNoArgOut(nargout);
            
            if ~iscell(trialData)
                trialData = {trialData};
            end
           
            if ~isempty(varargin) && isa(varargin{1}, 'TrialData')
                error('Pass all TrialData as a cell array to .addTrialsFromTrialData');
            end
            
            % add the channels
            for i = 1:numel(trialData)
                td = td.addChannelsMissing(trialData{i}.channelDescriptorsByName);
            end
            
            % add the data
            for i = 1:numel(trialData)
                td = td.addNewTrialsRaw(trialData{i}, varargin{:});
            end
        end

        function td = addNewTrialsRaw(td, dataOrTrialData, varargin)
            td.warnIfNoArgOut(nargout);
            p = inputParser();
            p.addParamValue('suppressWarnings', false, @islogical);
            p.parse(varargin{:});
            
            concatData = td.data;
            concatValid = td.manualValid;

            if isa(dataOrTrialData, 'TrialData')
                tdNew = dataOrTrialData;
                newData = tdNew.data;
            else
                tdNew = [];
                newData = dataOrTrialData;
            end

            % validate new data against all channel descriptors (old + new)
            %debug('Validating new channel data...\n');
            newData = td.validateData(newData, td.channelDescriptorsByName, 'addMissingFields', true, ...
                'suppressWarnings', p.Results.suppressWarnings);

            % concatenate onto the old data
            concatData = TrialDataUtilities.Data.structcat(concatData, newData); 

            % concatenate the valid array
            if isempty(tdNew)
                newValid = truevec(numel(newData));
            else
                newValid = tdNew.manualValid;
                assert(numel(newValid) == numel(newData), 'manualValid must be same length as data');
            end
            concatValid = cat(1, concatValid, newValid);
            
            td.data = concatData;
            td.manualValid = concatValid;
            
            td = td.postAddNewTrials(); % update valid here, also to allow TDCA to update itself before update valid since TDCA's override executes first
        end
        
        function td = postAddNewTrials(td)
            td.warnIfNoArgOut(nargout);
            td = td.invalidateValid();
        end

        % performs a validation of all channel data against the specified channelDescriptors,
        % also fixing empty values appropriately. We will also adjust
        % channelDescriptors whose data classes do not match the data in
        % data (having convert data to match that class)
        function [data, channelDescriptorsByName] = validateData(td, data, channelDescriptorsByName, varargin) %#ok<INUSL>
            p = inputParser();
            p.addParameter('addMissingFields', false, @islogical); % if true, don't complain about missing channels, just add the missing fields
            p.addParamValue('suppressWarnings', false, @islogical); % don't warn about any minor issues
            p.parse(varargin{:});
            suppressWarnings = p.Results.suppressWarnings;

            names = fieldnames(channelDescriptorsByName);
            nChannels = numel(names); %#ok<*PROPLC>

            % loop over channels and verify
            %fprintf('Validating channel data...\n');
            ok = falsevec(nChannels);
            required = falsevec(nChannels);
            missing = cellvec(nChannels);
            chDescs = cellvec(nChannels);
            for iChannel = 1:nChannels
                name = names{iChannel};
                chd = channelDescriptorsByName.(name);
                
                % check to make sure all fields were provided as expected 
                [ok(iChannel), missing{iChannel}] = chd.checkData(data); 
                chDescs{iChannel} = chd.describe();
                required(iChannel) = chd.required;
                
                if ~ok(iChannel)
                    data = chd.addMissingFields(data);
                    if ~suppressWarnings
                        if ~chd.required
                            % fill in missing values for optional channels but
                            % issue a warning
                            tcprintf('inline', '{yellow}Warning: {none}Missing optional channel {light blue}%s {none}fields {purple}%s\n', ...
                                chd.describe(), strjoin(missing{iChannel}, ', '));
                        else
                            tcprintf('inline', '{yellow}Warning: {none}Missing required channel {light blue}%s {none}fields {purple}%s\n', ...
                                chd.describe(), strjoin(missing{iChannel}, ', '));
                        end
                    end
                    if p.Results.addMissingFields || ~chd.required
                        ok(iChannel) = true; % mark as okay since not required
                    end 
                end
            end
            
            if any(~ok)
                missing = missing(~ok);
                chDescs = chDescs(~ok);
                for i = 1:length(missing)
                    tcprintf('inline', '{red}Error:   {none}Missing required channel {light blue}%s {none}fields {purple}%s\n', ...
                        chDescs{i}, strjoin(missing{i}, ', '));
                end
                error('Required channel data fields not provided by getChannelData');
            end

            % here we both change classes of data in memory and update
            % channelDescriptors depending on what we find
            prog = ProgressBar(nChannels, 'Repairing and converting channel data');
            for iChannel = 1:nChannels
                prog.update(iChannel);
                chd = channelDescriptorsByName.(names{iChannel}); 
                [data, chd] = chd.repairData(data); 
                data = chd.convertDataToMemoryClass(data);
                channelDescriptorsByName.(names{iChannel}) = chd;
            end
            prog.finish();
        end
        
        function v = getMetaByKey(td, key)
            if isstruct(td.datasetMeta) && isfield(td.datasetMeta, key)
                v = td.datasetMeta.(key);
            else
                v = [];
            end
        end
        
        function td = setMetaByKey(td, key, value)
            td.warnIfNoArgOut(nargout);
            assert(ischar(key), 'Key must be string');
            
            if ~isstruct(td.datasetMeta)
                td.datasetMeta = struct(key, value);
            else
                td.datasetMeta.(key) = value;
            end
        end
    end
    
    % Faster saving
    methods
        function saveFast(td, location)
            % saves in a custom hdf5 format that allows partial loading
            data = td.data;
            td.data = 'saved separately';
            td.odc = []; %#ok<MCHV2>
            
            mkdirRecursive(location);
            savefast(fullfile(location, 'td.mat'), 'td');
            
            % save elements of data
            msg = sprintf('Saving TrialData to %s', location);
            TrialDataUtilities.Data.SaveArrayIndividualized.saveArray(location, data, 'message', msg);
        end
    end
    
    methods(Static)
        function td = loadobj(s)
            if ~isa(s, 'TrialData')
                td = builtin('loadobj', s);
            else
                td = s;
                if isempty(td.odc)
                    td = td.rebuildOnDemandCache();
                end
            end
        end
           
        function td = loadFast(location)
            % strip extension
%             [path, name, ext] = fileparts(location);
%             location = fullfile(path, name);
            
            if ~exist(location, 'dir')
                error('Directory %s not found. Did you save with saveFast?', location);
            end
            loaded = load(fullfile(location, 'td.mat'));
            td = loaded.td;
            
            % load elements of data
            msg = sprintf('Loading TrialData from %s', location);
            td.data = TrialDataUtilities.Data.SaveArrayIndividualized.loadArray(location, 'message', msg);
            
            td = td.rebuildOnDemandCache();
        end

        % general utility to send plots to the correct axis
        function [axh, unmatched] = getRequestedPlotAxis(varargin)
            if isa(varargin{1}, 'TrialData') % used to be non-static method
                varargin = varargin(2:end);
            end
            
            p = inputParser();
            p.addParameter('figh', [], @(x) isempty(x) || ishandle(x));
            p.addParameter('axh', [], @(x) isempty(x) || ishandle(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if ~isempty(p.Results.axh)
                % if specified, use that axis
                axh = p.Results.axh;
                if ~ishold(axh)
                    cla(axh);
                end
            elseif ~isempty(p.Results.figh)
                % if figure specified, use that figure's current axis
                axh = get(p.Results.figh, 'CurrentAxes');
                if isempty(axh)
                    axh = axes('Parent', figh);
                end
            else
                % if neither, use newplot
                axh = newplot();
            end
            
            unmatched = p.Unmatched;
        end
        
        function tdCombined = mergeTrialsFromMultipleTrialData(varargin)
            if isempty(varargin)
                error('Please provide at least 1 argument');
            end
            tdCombined = varargin{1};
            if numel(varargin) > 1
                % pass as a cell array of trial data
                tdCombined = tdCombined.addTrialsFromTrialData(varargin(2:end));
            end
        end
    end

    % General utilities
    methods
        function printDescriptionShort(td)
            if td.nTrialsValid < td.nTrials
                tcprintf('inline', '{yellow}%s: {none}%d trials {bright red}(%d valid){none}, %d permanently invalid, with %d channels\n', ...
                class(td), td.nTrials, td.nTrialsValid, td.nTrialsPermanentlyInvalid, td.nChannels);
            else
                tcprintf('inline', '{yellow}%s: {none}%d trials (%d valid) with %d channels\n', ...
                    class(td), td.nTrials, td.nTrialsValid, td.nChannels);
            end
            if ~isempty(td.datasetName)
                tcprintf('inline', '{yellow}Dataset: {none}%s\n', td.datasetName);
            end
        end
        
        function printChannelInfo(td)
            % parse analog channels into grouped and non grouped
            analogCh = td.listNonContinuousNeuralAnalogChannels();
            continuousCh = td.listContinuousNeuralChannels();
            groups = td.listAnalogChannelGroups();
            groupChannels = cellvec(numel(groups));
            for iG = 1:numel(groups)
                groupChannels{iG} = td.listAnalogChannelsInGroup(groups{iG});
            end
            allGroupChannels = cat(1, groupChannels{:});

            analogChNonGroup = setdiff(analogCh, allGroupChannels);
            continuousChNonGroup = setdiff(continuousCh, allGroupChannels);
            
            tcprintf('inline', '{yellow}Param: {none}%s\n', strjoin(td.listParamChannels(), ', '));
            tcprintf('inline', '{yellow}Event: {none}%s\n', strjoin(td.listEventChannels(), ', '));
            tcprintf('inline', '{yellow}Analog: {none}%s\n', strjoin(analogChNonGroup, ', '));
            for iG = 1:numel(groups)
                isect = groupChannels{iG}(ismember(groupChannels{iG}, analogCh)); % want them in the group order
                if ~isempty(isect)
                    tcprintf('inline', '  {blue}%s: \{{none}%s\}\n', groups{iG}, strjoin(isect, ', '));
                end
            end
            tcprintf('inline', '{yellow}Spike: {none}%s\n', strjoin(td.listSpikeChannels(), ', '));
            tcprintf('inline', '{yellow}Continuous Neural: {none}%s\n', strjoin(continuousChNonGroup, ', '));
            for iG = 1:numel(groups)
                isect = groupChannels{iG}(ismember(groupChannels{iG}, continuousCh)); % want them in the group order
                if ~isempty(isect)
                    tcprintf('inline', '  {bright blue}%s: {none}\{%s\}\n', groups{iG}, strjoin(isect, ', '));
                end
            end
        end

        function disp(td)
            td.printDescriptionShort();
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end
        
        function td = setTrialDescriptionExtraParams(td, params)
            td.warnIfNoArgOut(nargout);
            if isempty(params)
                td.trialDescriptionExtraParams = {};
            else
                if ischar(params)
                    params = {params};
                end
                assert(iscellstr(params));
                params = makecol(params);
                td.trialDescriptionExtraParams = params;
            end
        end
    end
    
    % get / set accessors that read / copy-then-write through to ODC
    methods
        function v = get.valid(td)
            v = td.odc.valid;            
            if isempty(v)
                td.buildValid();
                v = td.odc.valid;
            end
        end
        
        function td = set.valid(td, v)
            td.odc = td.odc.copy();
            td.odc.valid = v;
        end
        
        function v = get.invalidCause(td)
            v = td.odc.invalidCause;            
            if isempty(v)
                td.buildValid();
                v = td.odc.invalidCause;
            end
        end
        
        function td = set.invalidCause(td, v)
            td.odc = td.odc.copy();
            td.odc.invalidCause = v;
        end
    end
    
    % builder methods for ODC: do not copy
    methods(Access=protected)
        function buildValid(td)
            % builds .valid and .invalidCause
            %
            % compute the valid flag considering only trials marked as
            % manually invalid to be invalid. This will be overriden in
            % TDCA to consider the condition and align invalid as well
            [manualValid, manualCause] = td.getManualValid();
            
            % then require the trials not be marked 'temporarily invalid'
            [tempValid, tempCause] = td.getTemporaryValid();
            
            valid = manualValid & tempValid;
            cause = manualCause;
            cause(manualValid & ~tempValid) = tempCause(manualValid & ~tempValid);
            
            % override, don't write
            c = td.odc;
            c.valid = valid;
            c.invalidCause = cause;
            td.odc = c;
        end
    end
    
    % Dependent property implementations
    methods % get. accessors for above properties which simply refer to tdi.?
        function ch = get.ch(td)
            names = fieldnames(td.channelDescriptorsByName);
            for iCh = 1:numel(names)
                ch.(names{iCh}) = names{iCh};
            end
        end
        
        function nTrials = get.nTrials(td)
            nTrials = numel(td.data);
        end

        function nTrials = get.nTrialsPermanentlyInvalid(td)
            nTrials = nnz(~td.getManualValid);
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
            if ~isempty(td.manualInvalidCause)
                td.manualInvalidCause = td.manualInvalidCause(mask);
            end
                
            if isempty(td.temporaryValid)
                td.temporaryValid = [];
            else
                td.temporaryValid = td.temporaryValid(mask);
            end
            if ~isempty(td.temporaryInvalidCause)
                td.temporaryInvalidCause = td.temporaryInvalidCause(mask);
            end
            td = td.invalidateValid();
        end
        
        function td = selectValidTrials(td, mask)
            % mask indexes into set of valid trials, omit to keep all valid
            % trials
             td.warnIfNoArgOut(nargout);
             if nargin < 2
                 maskFull = td.valid;
             else
                 maskFull = mi2ui(mask, td.valid);
             end
             
             td = td.selectTrials(maskFull);
        end

        function td = markTrialsInvalid(td, mask, reason)
            td.warnIfNoArgOut(nargout);
            td.manualValid(mask) = false;
            if nargin > 2
                assert(ischar(reason));
                if isempty(td.manualInvalidCause)
                    td.manualInvalidCause = cellvec(td.nTrials);
                end
                td.manualInvalidCause(mask) = {reason};
            end
            td = td.invalidateValid();
        end
        
        function td = setManualValidTo(td, mask)
            td.warnIfNoArgOut(nargout);
            assert(isvector(mask) && islogical(mask) && numel(mask) == td.nTrials);
            td.manualValid = makecol(mask);
            td = td.invalidateValid();
        end
        
        function td = setAllValid(td)
            td.manualValid = true(td.nTrials, 1);
            td = td.invalidateValid();
        end
        
        function td = markTrialsTemporarilyInvalid(td, mask, reason)
            td.warnIfNoArgOut(nargout);
            
            tempValid = td.getTemporaryValid();
            tempValid(mask) = false;
            if nargin > 2
                assert(ischar(reason));
                td.temporaryInvalidCause(mask) = {reason};
            end
            td.temporaryValid = makecol(tempValid);
            td = td.invalidateValid();
        end
        
        function td = withTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            mask = TensorUtils.vectorIndicesToMask(mask, td.nTrials);
            td = td.markTrialsTemporarilyInvalid(~mask);
        end
        
        function td = setTrialsTemporarilyInvalid(td, mask)
            td.warnIfNoArgOut(nargout);
            
            td.temporaryValid = TensorUtils.vectorIndicesToMask(mask, td.nTrials);
            td = td.invalideValid();
        end
        
        function td = clearTrialsTemporarilyInvalid(td)
            td.warnIfNoArgOut(nargout);
            td.temporaryValid = [];
            td = td.invalidateValid();
        end
        
        function td = withAllTrials(td)
            td.warnIfNoArgOut(nargout);
            td = td.clearTrialsTemporarilyInvalid();
        end
        
        function td = reset(td)
            td.warnIfNoArgOut(nargout);
            % don't touch .manualValid. this is not consistent with what reset means
            % for TrialDataConditionAlign
            
            td = td.clearTrialsTemporarilyInvalid();
        end
    end
    
    methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect', ...
                    class(obj));
            end
        end
    end

    methods(Access=protected) % Utility methods   
        function [valid, cause] = getManualValid(td)
            if isempty(td.manualValid)
                valid = truevec(td.nTrials);
                cause = cellvec(td.nTrials);
            else
                valid = makecol(td.manualValid);
                
                % generate the list of causes
                if isempty(td.manualInvalidCause)
                    cause = cellvec(td.nTrials);
                else
                    cause = td.manualInvalidCause;
                end
                
                emptyMask = cellfun(@isempty, cause);
                cause(~valid & emptyMask) = {'marked invalid manually'};
            end
        end
        
        function [valid, cause] = getTemporaryValid(td)
            if isempty(td.temporaryValid)
                valid = truevec(td.nTrials);
                cause = cellvec(td.nTrials);
            else
                valid = makecol(td.temporaryValid);
                
                % generate the list of causes, prefixed with (temporary)
                if isempty(td.temporaryInvalidCause)
                    cause = cellvec(td.nTrials);
                else
                    cause = td.temporaryInvalidCause;
                end
                
                for iT = 1:td.nTrials
                    if ~valid(iT)
                        if isempty(cause{iT})
                            cause{iT} = '(temporary) marked invalid temporarily';
                        else
                            cause{iT} = ['(temporary) ' cause{iT}];
                        end
                    end
                end
            end
        end
        
        function vals = replaceInvalidMaskWithValue(td, vals, value)
            % (valid, :) notation is to allow vals to be high dimensional
            if iscell(vals)
                [vals{~td.valid, :}] = deal(value);
            else
                if isempty(value)
                    value = NaN;
                end
                vals(~td.valid, :) = value;
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end

    methods % Channel metadata access and manipulation
                
        function td = invalidateValid(td)
            td.warnIfNoArgOut(nargout);
            
            % copy and flush valid
            td.odc = td.odc.copy();
            td.odc.flushValid();
            
            if isempty(td.manualValid)
                td.manualValid = truevec(td.nTrials);
            end
        end
        
        function tf = hasChannel(td, name)
            tf = ismember(name, td.channelNames);
        end
        
        function assertHasChannel(td, name)
            if ischar(name)
                assert(td.hasChannel(name), 'TrialData does not have channel %s', name);
            elseif iscellstr(name)
                tf = td.hasChannel(name);
                missing = name(~tf);
                assert(all(tf), 'Trial data does not have channels %s', strjoin(missing, ','));
            else
                error('Name must be string or cellstr');
            end
        end
        
        function cd = getChannelDescriptor(td, name)
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
        end
        
        function cdCell = getChannelDescriptorMulti(td, names)
            td.assertHasChannel(names);
            cdCell = cellfun(@(name) td.channelDescriptorsByName.(name), names, 'UniformOutput', false);
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
        
        function td = setChannelMeta(td, name, meta)
            td.warnIfNoArgOut(nargout);
            td.channelDescriptorsByName.(name).meta = meta;
        end
        
        function meta = getChannelMeta(td, name)
            meta = td.channelDescriptorsByName.(name).meta;
        end
        
        function names = listChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            names = {channelDescriptors.name}';
        end
        
        function names = listChannelsMatchingWildcard(td, varargin)
            % search for string matching 
            for i = 1:numel(varargin)
                varargin{i} = regexptranslate('wildcard', varargin{i});
            end
            names = td.listChannelsMatchingRegexp(varargin{:});
        end
        
        function names = listChannelsMatchingRegexp(td, varargin)
            chList = td.listChannels();
            mask = falsevec(numel(chList));
            for i = 1:numel(varargin)
                starts = regexp(chList, varargin{i});
                mask = mask | ~cellfun(@isempty, starts);
            end
            names = chList(mask);
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
        
        function durations = getDurationsRaw(td)
            starts = td.getEventRawFirst('TrialStart');
            ends = td.getEventRawLast('TrialEnd');
            durations = ends - starts;
        end
        
        function durations = getValidDurations(td)
            % get the time window for each trial
            starts = td.getEventFirst('TrialStart');
            ends = td.getEventLast('TrialEnd');
            durations = ends-starts;
            durations = td.replaceInvalidMaskWithValue(durations, NaN);
        end
             
        function durations = getValidDurationsForSpikeChannel(td, unitName)
            % similar to getValidDurations, except factors in the blanking
            % region
            durations = td.getValidDurations();
            
            blankIntervals = td.getSpikeBlankingRegions(unitName);
            
            % blankIntervals are guaranteed to be non-overlapping and lie
            % within the alignment window, so we can just add up the
            % individual intervals
            totalFn = @(mat) sum(mat(:, 2) - mat(:, 1), 1);
            blankDurations = zerosvec(td.nTrials);
            for iT = 1:td.nTrials
                if ~isempty(blankIntervals{iT})
                    blankDurations(iT) = totalFn(blankIntervals{iT});
                end
            end
            
            durations = durations - blankDurations;
            if any(durations < 0)
                error('Internal issue with blanking region durations');
            end
        end 
        
        function cds = getChannelDescriptorArray(td)
            % this sort operation is what puts channel names in order
            % during display, don't remove
            fields = sort(fieldnames(td.channelDescriptorsByName));
            if isempty(fields)
                cds = [];
            else
                for iF = 1:length(fields)
                    cds(iF) = td.channelDescriptorsByName.(fields{iF}); %#ok<AGROW>
                end
            end
            cds = makecol(cds);
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
                % don't need the assert since this might be called before
                % the field is written
                if ~isfield(td.data, fields{iF}), continue; end
%                 assert(isfield(td.data, fields{iF}), 'TrialData does not have data field %s', fields{iF});
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
        
        function td = dropChannel(td, name)
            td.warnIfNoArgOut(nargout);
            td = td.dropChannels(name);
        end
        
        function td = dropAnalogChannelGroup(td, groupName)
            td.warnIfNoArgOut(nargout);
            if ~td.hasAnalogChannelGroup(groupName)
                return;
            end
            
            chList = td.listAnalogChannelsInGroup(groupName);
            td = td.dropChannels(chList);
        end
        
        function td = dropChannels(td, names)
            td.warnIfNoArgOut(nargout);
            
            if isempty(names)
                return;
            end
            
            if ischar(names)
                names = {names};
            end
            
            % don't remove special channels
            names = setdiff(names, td.listSpecialChannels());
            
            % don't remove channels that don't exist
            names = intersect(names, td.listChannels());
            if isempty(names)
                return;
            end
            
            % first hold onto the to-be-removed channel descriptors
            cds = cellvec(numel(names));
            groupNames = cellvec(numel(names));
            inGroup = falsevec(numel(names));
            for i = 1:numel(names)
                cds{i} = td.channelDescriptorsByName.(names{i});
                if isa(cds{i}, 'AnalogChannelDescriptor') && cds{i}.isColumnOfSharedMatrix
                    groupNames{i} = cds{i}.primaryDataField;
                    inGroup(i) = true;
                end
            end
            
            groupsAffected = unique(groupNames(inGroup));
            
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
                     
            % clean the analog channel groups
            for iG = 1:numel(groupsAffected)
                if isfield(td.data, groupsAffected{iG}) % could be that the whole group was dropped
                    td = td.removeUnusedColumnsAnalogChannelGroup(groupsAffected{iG});
                end
            end
            
            td = td.updatePostDataChange(fieldsRemove);
        end
    end
    
    methods
        function saveTags = listSaveTags(td)
            saveTags = td.getParamUnique('saveTag');
        end
        
        function td = selectTrialsFromSaveTag(td, saveTags)
            td.warnIfNoArgOut(nargout);
            
            mask = ismember(td.getParam('saveTag'), saveTags);
            td = td.selectTrials(mask);
        end
    end
    
    methods % Analog channel methods
        function td = addAnalog(td, name, varargin)
            % td = td.addAnalog(td, name, values, times)
            % values may be cell of vectors or matrix with each row
            % corresponding to each trial.
            % times may be cell of vectors or single vector (for the matrix
            % values). Alternatively, times may be blank, and then the
            % parameter 'timeField' can specify which time field in which
            % to find the times for this channel
            p = inputParser();
            p.addOptional('values', {}, @(x) iscell(x) || ismatrix(x));
            p.addOptional('times', {}, @(x) ischar(x) || iscell(x) || isvector(x)); % char or time cell
            p.addParameter('timeField', '', @ischar);
            p.addParameter('units', '', @ischar);
            p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're identical
            p.addParameter('isAligned', true, @islogical);
            p.addParameter('clearForInvalid', false, @islogical);
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);
            p.parse(varargin{:});
            times = p.Results.times;
            values = p.Results.values;
            units = p.Results.units;

            td.warnIfNoArgOut(nargout);
            
            if ischar(times)
                % allow specifying char in times as well as timeField
                timeField = times;
                times = [];
            else
                timeField = p.Results.timeField;
            end
            
            % remove NaN values from analog signals?
%             if ~iscell(values) && ismatrix(values) && ~isempty(values)
%                 assert(size(values, 1) == td.nTrials, 'Data matrix must be size nTrials along first dimension');
%                 
%                 % strip nans from values and times together
%                 if ~isempty(times)
%                     if isvector(times)
%                         % same time vector for each row
%                         [values, times] = arrayfun(@(idx) removenanBoth(values(idx, :), times), makecol(1:td.nTrials), 'UniformOutput', false);
%                     else
%                         [values, times] = arrayfun(@(idx) removenanBoth(values(idx, :), times(idx, :)), makecol(1:td.nTrials), 'UniformOutput', false);
%                     end
%                 else
%                     values = arrayfun(@(idx) removenan(values(idx, :)), makecol(1:td.nTrials), 'UniformOutput', false);
%                 end
%             end
            
            % times can either be a field/channel name, or it can be raw
            % time values
            if isempty(times)
                % reference another time field or another analog channels
                % time field
                if ~isempty(timeField)
                    if td.hasChannel(timeField)
                        % treat timeField as analog channel name
                        % share that existing channel's time field 
                        cd = td.channelDescriptorsByName.(timeField);
                        assert(isa(cd, 'AnalogChannelDescriptor'), ...
                            'Channel %s is not an analog channel', timeField);
                        timeField = cd.timeField;
                        times = {td.data.(timeField)};

                    elseif isfield(td.data, timeField)
                        % use directly specified time field in .data
                        times = {td.data.(timeField)};
                    else
                        % timeField not found
                        error('%s is not a channel or data field name', timeField);
                    end
                    
                else
                    error('Time vector or cell of vectors must be passed in when not referencing an existing timeField');
                end
            else
                % times were provided
                % we'll set the field value for the new times field
                
                if isempty(timeField)
                    % generate unique time field name
                    timeField = matlab.lang.makeUniqueStrings(sprintf('%s_time', name), fieldnames(td.data));
                else
                    % check that we're not overwriting another channel's
                    % time field
                    if isfield(td.data, timeField)
                        otherChannels = setdiff(td.getChannelsReferencingFields(timeField), name);
                        if ~isempty(otherChannels)
                            error('Analog channel time field %s conflicts with existing channel(s) %s. If you meant to reference their timeField, specify ''timeField'' parameter and leave times argument empty', ...
                                name, strjoin(otherChannels, ','));
                        end
                    end
                end
            end
                
            % build a channel descriptor for the data
            if p.Results.isContinuousNeural
                cd = ContinuousNeuralChannelDescriptor.buildVectorAnalogFromValues(name, timeField, units, td.timeUnitName, values, times);
            else
                cd = AnalogChannelDescriptor.buildVectorAnalogFromValues(name, timeField, units, td.timeUnitName, values, times);
            end
            
            cd.scaleFromLims = p.Results.scaleFromLims;
            cd.scaleToLims = p.Results.scaleToLims;
            
            td = td.addChannel(cd);
            
            if ~isempty(values)
                td = td.setAnalog(name, values, times, ...
                    'clearForInvalid', p.Results.clearForInvalid, 'isAligned', p.Results.isAligned);
            end
            
%             function [t, v] = removenanBoth(t, v)
%                 m = ~isnan(t) & ~isnan(v);
%                 t = t(m);
%                 v = v(m);
%             end
        end
        
        function td = addContinuousNeural(td, name, varargin)
            % see addAnalog, same signature
            td.warnIfNoArgOut(nargout);
            td = td.addAnalog(name, varargin{:}, 'isContinuousNeural', true);
        end
        
        function td = setAnalog(td, name, values, varargin)
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);
            
            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) ||  isnumeric(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', false, @islogical);
            
            % for internal use mostly, used to indicate that the data provided is for the full trial, 
            % not just the aligned window. this is for when times are not
            % provided, then when this flag is true, we know that the times
            % are not being updated (since we're not losing times outside
            % the alignment window). Used currently by
            % setAnalogWithinAlignWindow
            p.addParameter('timesMatchFullTrialRaw', false, @islogical); 
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            times = makecol(p.Results.times);
          
            % check the values and convert to nTrials cellvec
            if ismatrix(values) && isnumeric(values)
                % values must be nTrials x nTimes
                assert(size(values, 1) == td.nTrials, 'Values as matrix must be nTrials along dimension 1');
                values = mat2cell(values', size(values, 2), onesvec(td.nTrials))';
                
            elseif iscell(values)
                assert(numel(values) == td.nTrials, 'Values as cell must have numel == nTrials');
                values = makecol(values);
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
                updateTimes = true;
            else
                % if no times are specified, then we're only updating the
                % timestamps if we're aligned to a subset of the trial
                % currently, in which case we're dropping the timepoints
                % outside the alignment window
                if p.Results.timesMatchFullTrialRaw
                    updateTimes = false;
                    times = td.getAnalogTimeRaw(name);
                else
                    if td.alignIncludesFullTrial()
                        % no need to update timestamps, but we do need to
                        % trim the timestamps to TrialStart:TrialEnd since
                        % that's what we expect to receive back 
                        updateTimes = false;
                        td = td.trimAnalogChannelToTrialStartEnd(name);
                    else
                        updateTimes = true;
                    end
                
                    % pass along the current times since the data is coming in with the
                    % existing alignment
                    times = td.getAnalogTime(name);
                end
            end
            
            % check that times have same length as data
            assert(numel(times) == numel(values), 'Times and values must have same number of trials');
            nTimes = cellfun(@numel, times);
            nValues = cellfun(@numel, values);
            assert(all(nTimes == nValues), 'Mismatch between number of times and values. If the number of times has changed be sure to specify times parameter');
            
            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            if p.Results.isAligned && ~p.Results.timesMatchFullTrialRaw
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);

            cd = td.channelDescriptorsByName.(name);
                        
            % for shared column channels where the scaling or time vectors
            % change, it needs to be separated from the shared column form
            if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix && ((cd.hasScaling && ~p.Results.keepScaling) || updateTimes)
                error('Setting analog channel while ''keepScaling'' is false or when specifying new sample times requires this channel to be separated from its analog channel group. Use separateAnalogChannelFromGroup if you want to do this. Or use setAnalogChannelGroup to set all channels in the group at once.');
%                 debug('Separating channel %s from column of shared matrix %s\n', name, cd.primaryDataField);
%                 td = td.separateAnalogChannelFromGroup(name, false, false); % no separate time field, no copy data
%                 cd = td.channelDescriptorsByName.(name);
            end
            
            if ~p.Results.keepScaling
                % data being passed in is now in original units
                % so change scaling factors
                cd = cd.withNoScaling();
            else
                % take new data back into scaled values
                values = cd.convertAccessDataCellToMemory(1, values);
            end
            
            % update the channel descriptor accordingly
            td.channelDescriptorsByName.(name) = cd;
            
            % setChannelData will call repairData which will update
            % memoryDataClassByField{1} to reflect the type of values
            if updateTimes
                td = td.setChannelData(name, {values, times}, p.Unmatched);
            else
                td = td.setChannelData(name, {values}, p.Unmatched);
            end
        end
        
        function td = convertAnalogChannelToNoScaling(td, name)
            % convert a single channel to using scaling (representing the
            % data in memory in a different scaling, and often a different
            % data class like uint16) to one not using scaling
            % (representing the data in memory in raw units)
            
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);

            cd = td.channelDescriptorsByName.(name);
            
            if ~cd.hasScaling
                % nothing to do
                return;
            end

            % if its part of a shared column group of analog channels,
            % separate it first
            if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix
                td = td.separateAnalogChannelFromGroup(name, false);
                cd = td.channelDescriptorsByName.(name);
            end
            
            [data, time] = td.getAnalogRaw(name);
            
            % data being passed in is now in original units
            % so change scaling factors
            cd = cd.withNoScaling();
            
            % update the channel descriptor accordingly
            td.channelDescriptorsByName.(name) = cd;
            
            % use setChannel data, which will call repairData to change
            % data classes of each field
            td = td.setChannelData(name, {data, time}, 'fieldMask', [true false]);
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
            if isempty(channelDescriptors)
                names = {};
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'AnalogChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function names = listNonContinuousNeuralAnalogChannels(td)
            names = setdiff(td.listAnalogChannels(), td.listContinuousNeuralChannels());
        end
        
        function names = listContinuousNeuralChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            if isempty(channelDescriptors)
                names = {};
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'ContinuousNeuralChannelDescriptor'), channelDescriptors);
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
        
        function fsHz = getAnalogSamplingRateHz(td, name)
            fsHz = td.timeUnitsPerSecond / td.getAnalogTimeDelta(name);
        end
        
        function time = getAnalogTimeRaw(td, name)
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            assert(isa(cd, 'AnalogChannelDescriptor'), 'Channel %s is not analog', name);
            time = {td.data.(cd.dataFields{2})}';
        end
        
        function [data, time] = getAnalogRaw(td, name, varargin)
            p = inputParser();
            p.addParameter('sort', false, @islogical);
            p.addParameter('applyScaling', true, @islogical);
            p.parse(varargin{:});
            
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            assert(isa(cd, 'AnalogChannelDescriptor'), 'Channel %s is not analog', name);
            
            if cd.isColumnOfSharedMatrix
                data = arrayfun(@(t) t.(cd.dataFields{1})(:, cd.primaryDataFieldColumnIndex), ...
                    td.data, 'UniformOutput', false, 'ErrorHandler', @(varargin) []);
            else
                data = {td.data.(cd.dataFields{1})}';
            end
            time = {td.data.(cd.dataFields{2})}';
            for i = 1:numel(data)
                if numel(data{i}) == numel(time{i}) - 1
                    time{i} = makecol(time{i}(1:end-1));
                else
                    time{i} = makecol(time{i});
                end
                data{i} = makecol(data{i});
            end
            
            if p.Results.applyScaling
                % do scaling and convert to double
                data = cd.convertDataCellOnAccess(1, data);
            end
            
            data = makecol(data);
            
            if p.Results.sort
                for iT = 1:td.nTrials
                    [time{iT}, idx] = sort(time{iT}, 'ascend');
                    data{iT} = data{iT}(idx);
                end
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
       
        function [dataUnif, timeUnif, delta] = getAnalogRawUniformlySampled(td, name, varargin)
            p = inputParser();
            p.addParameter('method', 'linear', @ischar);
            p.addParameter('delta', [], @(x) isscalar(x) || isempty(x));
            p.parse(varargin{:});
            
            delta = p.Results.delta;
            if isempty(delta)
                delta = td.getAnalogTimeDelta(name);
            end
            [data, time] = td.getAnalogRaw(name);  %#ok<*PROP>
            
            [dataUnif, timeUnif] = cellvec(td.nTrials);
            
            % get first and last valid sample for all trials
            % very important to use this in case data has NaN samples
            [tmins, tmaxs] = TrialDataUtilities.Data.getValidTimeExtents(time, data);
            for iT = 1:td.nTrials
                if isempty(time{iT}) || isempty(data{iT})
                    continue;
                end
                timeUnif{iT} = (tmins(iT):delta:tmaxs(iT))';
                if nnz(~isnan(data{iT})) > 5
                    dataUnif{iT} = interp1(time{iT}, data{iT}, timeUnif{iT}, p.Results.method);
                else
                    dataUnif{iT} = nan(size(timeUnif{iT}));
                end
            end
        end
        
        function [dataUnif, timeUnif, delta] = getAnalogUniformlySampled(td, name, varargin)
            [dataUnif, timeUnif, delta] = td.getAnalogRawUniformlySampled(name, varargin{:});
            dataUnif = td.replaceInvalidMaskWithValue(dataUnif, []);
            timeUnif = td.replaceInvalidMaskWithValue(timeUnif, []);
        end
        
        function [dataCell, timeCell] = getAnalogMulti(td, name, varargin)
            % [data, time] = getAnalogMulti(td, chNames, varargin)
            % data and time are cell(nTrials, nChannels)
            p = inputParser();
            p.addParameter('raw', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            if ischar(name)
                name = {name};
            end
            
            % build nTrials x nChannels cell of data/time vectors
            C = numel(name);
            [dataCell, timeCell] = deal(cell(td.nTrials, C));
            prog = ProgressBar(C, 'Fetching analog channels');
            for c = 1:C
                prog.update(c);
                if p.Results.raw
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalogRaw(name{c});
                else
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalog(name{c}, varargin{:});
                end
            end
            prog.finish();
        end
    end
    
    methods % Analog Group and multi-channel analog methods
        function groupNames = listAnalogChannelGroups(td)
            chList = td.listAnalogChannels();
            
            mask =falsevec(numel(chList));
            groupNames = cellvec(numel(chList));
            for iC = 1:numel(chList)
                [mask(iC), groupNames{iC}] = td.isAnalogChannelInGroup(chList{iC});
            end
            
            groupNames = unique(groupNames(mask));
        end
        
        function tf = checkAnalogChannelGroupHasUniformScaling(td, groupName)
            td.assertHasAnalogChannelGroup(groupName);
            
            chList = td.listAnalogChannelsInGroup(groupName);
            cdCell = td.getChannelDescriptorMulti(chList);
            
            tf = true;
            hasScaling = cdCell{1}.hasScaling;
            from = cdCell{1}.scaleFromLims;
            to = cdCell{1}.scaleToLims;
            
            for i = 2:numel(chList)
                cd = cdCell{i};
                if cd.hasScaling ~= hasScaling || ~isequal(cd.scaleFromLims, from) || ~isequal(cd.scaleToLims, to)
                    tf = false;
                    break;
                end
            end
        end
       
       function td = convertAnalogChannelGroupToNoScaling(td, groupName, varargin)
            % convert a set of channels sharing a single data field (columns in that matrix)
            % that use scaling (representing the data in memory in a different scaling)
            % and often a different data class like uint16) to one not using scaling
            % (representing the data in memory in raw units). Doing this
            % for all shared channels simultaneously allows us to preserve
            % the shared matrix structure
            
            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelGroup(groupName);

            chList = td.listAnalogChannelsInGroup(groupName);
            cdCell = td.getChannelDescriptorMulti(chList);
            
            hasScaling = cellfun(@(cd) cd.hasScaling, cdCell);
            if ~any(hasScaling)
                % nothing to do
                return;
            end
            
            % don't want to handle the more complicated edge case where the
            % channels don't have uniform scaling at this point
            uniform = td.checkAnalogChannelGroupHasUniformScaling(groupName);
            if ~uniform
                error('Not supported for analog channel groups whose channels are not uniformly scaled');
            end
            
            newCdCell = cellfun(@(cd) cd.withNoScaling(), cdCell, 'UniformOutput', false);
            
            % manually replace the .data field with the new values
            % we convert memory->access using the old cd, then
            % access->memory with the new cd.
            cdOld = cdCell{1};
            dataField = cdOld.primaryDataField;
            cdNew = newCdCell{1};
            prog = ProgressBar(td.nTrials, 'Converting %s data to unscaled form', groupName);
            for t = 1:td.nTrials
                prog.update(t);
                unscaled = cdOld.convertDataSingleOnAccess(1, td.data(t).(dataField));
                td.data(t).(dataField) = cdNew.convertAccessDataSingleToMemory(1, unscaled);
            end
            prog.finish();

            % and replace all the channel descriptors
            for i = 1:numel(newCdCell)
                td.channelDescriptorsByName.(newCdCell{i}.name) = newCdCell{i};
            end
        end
        
        function [tf, groupName] = checkAnalogChannelsInSameGroup(td, names)
            groupNames = cellfun(@(name) td.getAnalogChannelGroupName(name), names, 'UniformOutput', false);
            tf = numel(unique(groupNames)) == 1 && ~isempty(groupNames{1});
            if tf
                groupName = groupNames{1};
            else
                groupName = '';
            end
        end
        
        function cd = getAnalogChannelGroupSingleChannelDescriptor(td, groupName)
            chList = td.listAnalogChannels();
            for iC = 1:numel(chList)
                cd = td.channelDescriptorsByName.(chList{iC});  
                if cd.isColumnOfSharedMatrix && strcmp(cd.primaryDataField, groupName)
                    % in the group, return it
                    return;  
                end
            end
            error('Analog channel group not found');
        end
        
        function [names, colIndex] = listAnalogChannelsInGroup(td, groupName)
            td.assertHasAnalogChannelGroup(groupName);
            chList = td.listAnalogChannels();
            cdCell = cellvec(numel(chList));
            mask = falsevec(numel(chList));
            for iC = 1:numel(chList)
                cd = td.channelDescriptorsByName.(chList{iC});  
                if cd.isColumnOfSharedMatrix && strcmp(cd.primaryDataField, groupName)
                    mask(iC) = true;
                    cdCell{iC} = cd;     
                end
                
            end
            names = chList(mask);
            cdCell = cdCell(mask); 
            
            % order them by their column index
            colIndex = cellfun(@(cd) cd.primaryDataFieldColumnIndex, cdCell);
            [colIndex, idx] = sort(colIndex, 'ascend');
            names = names(idx);
        end
        
        function [chList, colIdx, groupName] = listAnalogChannelsInGroupWith(td, name)
            % chList = listAnalogChannelsInGroupWith(td, name)
            % chList includes name (or is empty if not in any group)
            td.assertHasChannel(name);
            groupName = td.getAnalogChannelGroupName(name);
            if isempty(groupName)
                chList = {};
                colIdx = [];
            else
                [chList, colIdx] = td.listAnalogChannelsInGroup(groupName);
            end
        end
        
        function [colIdx, groupName] = getAnalogChannelColumnIdxInGroup(td, names)
            [tf, groupName] = td.checkAnalogChannelsInSameGroup(names);
            assert(tf, 'Channels are not in the same analog channel group');
            
            cdCell = td.getChannelDescriptorMulti(names);
            colIdx = cellfun(@(cd) cd.primaryDataFieldColumnIndex, cdCell);
        end
        
        function [tf, groupName] = isAnalogChannelInGroup(td, name, groupName)
            % [tf, groupName] = isAnalogChannelInGroup(td, name) % in any group?
            % tf = isAnalogChannelInGroup(td, groupName) % in specific group?
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            if nargin > 2
                tf = cd.isColumnOfSharedMatrix & strcmp(cd.primaryDataField, groupName);
            else
                tf = cd.isColumnOfSharedMatrix;
            end
            if tf
                groupName = cd.primaryDataField;
            else
                groupName = '';
            end
        end
        
        function groupName = getAnalogChannelGroupName(td, name)
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            
            if ~cd.isColumnOfSharedMatrix
                groupName = '';
            else
                groupName = cd.primaryDataField;
            end
        end

        function td = separateAnalogChannelFromGroup(td, name, varargin)
            p = inputParser();
            p.addParameter('separateTimeField', false, @islogical);
            p.parse(varargin{:});
            
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);
            
            groupName = td.getAnalogChannelGroupName(name);
            
            if isempty(groupName)
                % not in group
                return;
            end
            
            cd = td.channelDescriptorsByName.(name);
            
            % fetch current data without scaling
            oldData = td.getAnalogRaw(name, 'applyScaling', false);

            assert(~isfield(td.data, name), 'Issue with creating field %s already found in td.data', name); 

            if p.Results.separateTimeField
                newTimeField = matlab.lang.makeUniqueStrings([name '_time'], fieldnames(td.data));
                td.data = copyStructField(td.data, td.data, cd.timeField, newTimeField);
                cdNew = cd.separateFromColumnOfSharedMatrix(newTimeField);
            else
                cdNew = cd.separateFromColumnOfSharedMatrix();
            end

            % we also need to copy the data field over, since we may
            % only be updating valid trials
            
            [td.data.(cdNew.primaryDataField)] = deal(oldData{:});
            
            td.channelDescriptorsByName.(name) = cdNew;
            
            if td.hasAnalogChannelGroup(groupName)
                % if this were the last channel in the group, the group
                % will no longer "exist"
                td = td.removeUnusedColumnsAnalogChannelGroup(groupName);
            end
        end
        
        function td = removeUnusedColumnsAnalogChannelGroup(td, groupName)
            td.warnIfNoArgOut(nargout);
            [chList, colIdx] = td.listAnalogChannelsInGroup(groupName);
            
            prog = ProgressBar(td.nTrials, 'Cleaning columns of analog channel group data');
            for t = 1:td.nTrials
                prog.update(t);
                
                if ~isempty(td.data(t).(groupName))
                    td.data(t).(groupName) = td.data(t).(groupName)(:, colIdx);
                end
            end
            prog.finish();
            
            % and update the channel descriptors
            for c = 1:numel(chList)
                td.channelDescriptorsByName.(chList{c}).primaryDataFieldColumnIndex = c;
            end
            
            td = td.updatePostDataChange(groupName);
        end
        
        function td = separateAnalogChannelsIntoSeparateGroup(td, names, newGroupName)
            td.warnIfNoArgOut(nargout);
            assert(td.checkAnalogChannelsInSameGroup(names), 'All channels must be part of a group already');
            
            [chAll, ~, groupName] = td.listAnalogChannelsInGroupWith(names{1});
            if ~strcmp(groupName, newGroupName)
                % name is changing
                conflicts = td.getChannelsReferencingFields(newGroupName);
                assert(isempty(conflicts), 'Channels %s already reference field %s', strjoin(conflicts, ', '), newGroupName);
            
            elseif ~isempty(setdiff(chAll, names))
                error('When group name is not changing, all channels must be included in the new ordering');
            end
            
            if isequal(names, chAll) && strcmp(groupName, newGroupName)
                % we only do nothing if the list matches in order, since
                % the user may expect the channels to be in this specific
                % order from now on
                return;
            end
            
            colIdx = td.getAnalogChannelColumnIdxInGroup(names);
            
            prog = ProgressBar(td.nTrials, 'Splitting analog channel group %s');
            for iT = 1:td.nTrials
                prog.update(iT);
                td.data(iT).(newGroupName) = td.data(iT).(groupName)(:, colIdx);
            end
            prog.finish();
            
            % update the channel descriptors to point to the new group
            for c = 1:numel(names)
                cd = td.channelDescriptorsByName.(names{c});
                cd.primaryDataField = newGroupName;
                cd.primaryDataFieldColumnIndex = c;
                td.channelDescriptorsByName.(names{c}) = cd;
            end
            
            % clear out the old channel
            if ~strcmp(groupName, newGroupName)
                if isempty(setdiff(chAll, names))
                    % no need for old field 
                    td.data = rmfield(td.data, groupName);
                else
                    % clean out the unused columns from the old field
                    td = td.removeUnusedColumnsAnalogChannelGroup(groupName);
                end
            end
            
            td = td.updatePostDataChange(newGroupName);
        end

        function td = renameAnalogChannelGroup(td, groupName, newGroupName)
            td.warnIfNoArgOut(nargout);
            chNames = td.listAnalogChannelsInGroup(groupName);
            td = td.separateAnalogChannelsIntoSeparateGroup(chNames, newGroupName);
        end
        
        function td = sortChannelsInAnalogChannelGroup(td, groupName, chNamesInOrder)
            td.warnIfNoArgOut(nargout);
            if nargin < 3
                % default to alpha sort
                chNames = td.listAnalogChannelsInGroup(groupName);
                chNamesInOrder = sort(chNames);
            end
            td = td.separateAnalogChannelsIntoSeparateGroup(chNamesInOrder, groupName);
        end
        
        function tf = hasAnalogChannelGroup(td, groupName)
            if ~isfield(td.data, groupName)
                tf = false;
            else
                % check if there are any channels in that group, do not
                % call listAnalogChannelsInGroup b/c of infinite recursion
                chList = td.listAnalogChannels();
                for iC = 1:numel(chList)
                    cd = td.channelDescriptorsByName.(chList{iC});  
                    if cd.isColumnOfSharedMatrix && strcmp(cd.primaryDataField, groupName)
                        tf = true;
                        return;
                    end
                end
                tf = false;
            end
        end 
        
        function assertHasAnalogChannelGroup(td, groupName)
            assert(td.hasAnalogChannelGroup(groupName), 'No analog channel group %s found', groupName);
        end
        
        function time = getAnalogChannelGroupTime(td, groupName)
            time = td.getAnalogChannelGroupTimeRaw(groupName);
            time = td.replaceInvalidMaskWithValue(time, []);
        end
        
        function time = getAnalogChannelGroupTimeRaw(td, groupName)
            td.assertHasAnalogChannelGroup(groupName);
            chList = td.listAnalogChannelsInGroup(groupName);
            cd = td.channelDescriptorsByName.(chList{1});
            timeField = cd.timeField;
            time = {td.data.(timeField)}';
        end
        
        function [data, time] = getAnalogChannelGroupRaw(td, groupName, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical);
            p.parse(varargin{:});
            
            td.assertHasAnalogChannelGroup(groupName);
            
            chList = td.listAnalogChannelsInGroup(groupName);
            cd = td.channelDescriptorsByName.(chList{1});
            timeField = cd.timeField;
            
            data = {td.data.(groupName)}';
            time = {td.data.(timeField)}';
            
            for i = 1:numel(data)
                time{i} = makecol(time{i});
                assert(size(data{i}, 1) == numel(time{i}), 'Number of timepoints in data on trial %d does not match time', i);
                assert(size(data{i}, 2) == numel(chList), 'Number of channels on trial %d does not match channel count', i);
            end
            
            % do scaling and convert to double
            if p.Results.applyScaling
                data = cd.convertDataCellOnAccess(1, data);
            end
        end
        
        function [dataMat, tvec] = getAnalogChannelGroupSingleTrial(td, groupName, trialInd, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical);
            p.parse(varargin{:});
            
            assert(isscalar(trialInd), 'TrialInd must be scalar single trial');
            
            td.assertHasAnalogChannelGroup(groupName);
            
            dataMat = td.data(trialInd).(groupName);
            if nargout > 1
                cd = td.getAnalogChannelGroupSingleChannelDescriptor(groupName);
                tvec = makecol(td.data(trialInd).(cd.timeField));
                assert(size(dataMat, 1) == numel(tvec), 'Number of timepoints in data on trial %d does not match time', trialInd);
            end
            
            % do scaling and convert to double
            if p.Results.applyScaling
                dataMat = cd.convertDataSingleOnAccess(1, dataMat);
            end
        end
        
        function [data, time] = getAnalogChannelGroup(td, groupName)
            [data, time] = td.getAnalogChannelGroupRaw(groupName);
            data = td.replaceInvalidMaskWithValue(data, []);
            time = td.replaceInvalidMaskWithValue(time, []);
        end
        
        function td = setAnalogChannelGroup(td, groupName, values, varargin)
            % values can be either
            %   - nTrials cellvec with nTime x nChannels matrices
            %     then 'times', should be nTrials cellvec of nTime vectors
            %   - nTrials x nTime x nChannels tensor
            %     then 'times' should be nTime vector
            %
            
            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelGroup(groupName);
            
            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) ||  ismatrix(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            times = makecol(p.Results.times);
            
            chList = td.listAnalogChannelsInGroup(groupName);
            nCh = numel(chList);
          
            % check the values and convert to nTrials cellvec
            if isnumeric(values) && isnumeric(values)
                % values must be nTrials x nTimes x nChannels
                assert(size(values, 1) == td.nTrials, 'Values as matrix must be nTrials along dimension 1');
                
                % convert to cellvec
                tensor = values;             
                values = cellvec(td.nTrials, 1);
                for r = 1:td.nTrials
                    values{r} = TensorUtils.squeezeDims(tensor(r, :, :), 1);
                end
                clear tensor;
                
            elseif iscell(values)
                assert(numel(values) == td.nTrials, 'Values as cell must have numel == nTrials');
                nCol = cellfun(@(x) size(x, 2), values);
                assert(all(nCol == nCh), 'All elements of values must have same number of columns as channels (%d)', nCh);
                values = makecol(values);
            else
                error('Values must be numeric tensor or cell array');
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
                updateTimes = true;
            else
                if td.alignIncludesFullTrial
                    % no need to update timestamps, but we do need to
                    % trim the timestamps to TrialStart:TrialEnd since
                    % that's what we expect to receive back 
                    updateTimes = false;
                    td = td.trimAnalogChannelGroupToTrialStartEnd(groupName);
                else
                    updateTimes = true;
                end
                
                % pass along the current times since the data is coming in with the
                % existing alignment
                times = td.getAnalogChannelGroupTime(groupName);
            end
            
            % check that times have same length as data
            assert(numel(times) == numel(values), 'Times and values must have same number of trials');
            nTimes = cellfun(@numel, times);
            nValues = cellfun(@(x) size(x, 1), values);
            assert(all(nTimes == nValues), 'Mismatch between number of times and values. If the number of times has changed be sure to specify times parameter');
            
            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);
          
            if ~p.Results.keepScaling
                % data being passed in is now in original units
                % so change scaling factors of the channel and the data
                td = td.convertAnalogChannelGroupToNoScaling(groupName);
                % and convert back to memory anyway in case the data class
                % has changed, although this shouldn't do any scaling
                values = td.channelDescriptorsByName.(chList{1}).convertAccessDataCellToMemory(1, values);
            else
                % take new data back into scaled values to match the
                % existing
                values = cd.convertAccessDataCellToMemory(1, values);
            end
            
            if updateTimes
                % avoid overwriting data shared by other channels,
                % including time field
                td = td.copyRenameSharedChannelFields(chList, 2);                
                timeField = td.channelDescriptorsByName.(chList{1}).timeField;
            end
            
            prog = ProgressBar(td.nTrials, 'Overwriting analog channel group data on valid trials');
            for t = 1:td.nTrials
                prog.update(t);
                if ~td.valid(t), continue; end
                
                td.data(t).(groupName) = values{t};
                if updateTimes
                    td.data(t).(timeField) = times{t};
                end
            end
            
            if updateTimes
                td = td.updatePostDataChange({groupName, timeField});
            else
                td = td.updatePostDataChange({groupName});
            end
            
            prog.finish();
        end
        
        function td = trimAnalogChannelGroupToTrialStartEnd(td, groupName)
            % Timepoints that lie outside of TrialStart and TrialStop will
            % never be accessible via getAnalog since they will be filtered
            % out by the AlignInfo. This can create problems when using
            % getAnalog to fetch data, and then setAnalog to change the
            % data in place without altering the timepoints, since the
            % timepoints outside of TrialStart:TrialStop will never be
            % considered. This method strips those "hidden" timepoints from
            % the analog channel group groupName
            
            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelGroup(groupName);
           
            chList = td.listAnalogChannelsInGroup(groupName);
            
            cd = td.channelDescriptorsByName.(chList{1});
            timeField = cd.timeField;
            
            td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeField);
        end
        
        function td = trimAnalogChannelToTrialStartEnd(td, chName)
            td.warnIfNoArgOut(nargout);
           
            td.assertHasChannel(chName);
            cd = td.channelDescriptorsByName.(chName);
            timeField = cd.timeField;
            
            td = td.trimAnalogChannelTimeFieldAndReferencingChannels(timeField);
        end
            
        function td = trimAnalogChannelTimeFieldAndReferencingChannelsRaw(td, timeField, varargin)
            % delete time points outside of a certain start stop interval
            % defaults to TrialStart:TrialStop
            p = inputParser();
            p.addOptional('startTimes', {}, @isvector);
            p.addOptional('stopTimes', {}, @isvector);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            times = {td.data.(timeField)}';
            
            % in TDCA this is the equivalent of getAlignedTimesCell, but we
            % don't have that infrastructure in TD, so we just do it here
            % directly
            if isempty(p.Results.startTimes)
                startTimes = cellfun(@(x) x(1), td.getEventRaw('TrialStart'));
            else
                startTimes = p.Results.startTimes;
            end
            if isempty(p.Results.stopTimes)
                stopTimes = cellfun(@(x) x(1), td.getEventRaw('TrialEnd'));
            else
                stopTimes = p.Results.stopTimes;
            end
            
            timesMask = cellvec(td.nTrials);
            for iT = 1:td.nTrials
                timesMask{iT} = times{iT} >= startTimes(iT) & times{iT} <= stopTimes(iT);
            end
            
            maskNeedsModification = ~cellfun(@all, timesMask);
            if ~any(maskNeedsModification)
                return;
            end
            
            % list channels that reference time field
            chList = td.getChannelsReferencingFields(timeField);
            cdCell = td.getChannelDescriptorMulti(chList);  
            
            % separate groups from non groups
            inGroup = cellfun(@(cd) cd.isColumnOfSharedMatrix, cdCell);
            groupList = unique(cellfun(@(cd) cd.primaryDataField, cdCell(inGroup), 'UniformOutput', false));
            chList = chList(~inGroup);
            
            prog = ProgressBar(td.nTrials, 'Stripping timepoints outside of trial for %s', timeField);
            for iT = 1:td.nTrials
                if ~maskNeedsModification(iT), continue; end
                prog.update(iT);
                td.data(iT).(timeField) = times{iT}(timesMask{iT});
                %
                % loop over channels
                for iC = 1:numel(chList)
                    td.data(iT).(chList{iC}) = td.data(iT).(chList{iC})(timesMask{iT});
                end
                % loop over channel groups
                for iG = 1:numel(groupList)
                    td.data(iT).(groupList{iG}) = td.data(iT).(groupList{iG})(timesMask{iT}, :);
                end
            end
            prog.finish();
        end
    end
    
    methods % Event channel methods
        function td = addEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addRequired('name', @ischar);
            p.addRequired('times', @(x) isempty(x) || isvector(x));
            p.addParameter('isAligned', true, @islogical);
            %p.addParamValue('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.parse(name, times, varargin{:});
            %cd = p.Results.channelDescriptor;
            
            if isempty(times)
                times = cellvec(td.nTrials);
            end
            if isscalar(times)
                times = repmat(times, td.nTrials, 1);
            end
            
            assert(numel(times) == td.nTrials, 'Times must be vector with length %d', td.nTrials);
            times = makecol(times);
            
            if iscell(times)
                % multiple occurrence event
                cd = EventChannelDescriptor.buildMultipleEvent(name, td.timeUnitName);
                
                % check contents
                assert(all(cellfun(@(x) isempty(x) || (isvector(x) && isnumeric(x)), times)), ...
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
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
                if iscell(times)
                    times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);
                else
                    times = times + offsets;
                end
            end
            
            if ~iscell(times)
                times = arrayfun(@(x) x(~isnan(x)), times, 'UniformOutput', false);
            end
            td = td.addChannel(cd, {times});
        end
        
        function td = addEventOccurrence(td, name, times, varargin)
            p = inputParser;
            p.addRequired('name', @ischar);
            p.addRequired('times', @(x) isempty(x) || isvector(x));
            p.addParameter('isAligned', true, @islogical);
            p.parse(name, times, varargin{:});
            
            assert(td.hasEventChannel(name), 'Event channel %s not found', name);
            
            % make times unaligned
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
                if iscell(times)
                    times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);
                else
                    times = times + offsets;
                end
            end
           
            % append the new time
            fullTimes = td.getEventRaw(name);
            for iT = 1:td.nTrials
                fullTimes{iT} = cat(1, fullTimes{iT}, times{iT});
            end
            
            td = td.setEvent(name, fullTimes, 'isAligned', false);
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
            if isempty(channelDescriptors)
                names = {};
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'EventChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        % used mainly by AlignInfo to make sure it can access unaligned
        % event info
        function eventStruct = getRawEventFlatStruct(td, chList)
            eventStruct = struct();
            
            if nargin < 2
                chList = td.listEventChannels();
            end
            
            for iCh = 1:numel(chList)
                ch = chList{iCh};
                if td.channelDescriptorsByName.(ch).isScalarByField(1)
                    eventStruct.(ch) = makecol([td.data.(ch)]);
                else
                    eventStruct.(ch) = cellRemoveNaN(makecol({td.data.(ch)}));
                end
            end
            
            function c = cellRemoveNaN(c)
                c = cellfun(@(x) x(~isnan(x)), c, 'UniformOutput', false);
            end
        end 
            
        function eventStruct = getRawEventStructArray(td)
            eventStruct = copyStructField(td.data, [], td.listEventChannels());
        end
                
        function timesCell = getEventRaw(td, name)
            timesCell = {td.data.(name)}';
            timesCell = cellfun(@(x) x(~isnan(x)), timesCell, 'UniformOutput', false);
        end
        
        function times = getEventRawFirst(td, name)
            timesCell = {td.data.(name)}';
            
            times = cellfun(@getFirst, timesCell);
            
            function t = getFirst(times)
                if ~isempty(times)
                    t = times(1);
                else
                    t = NaN;
                end
            end
        end
        
        function times = getEventRawLast(td, name)
            timesCell = {td.data.(name)}';
            
            times = cellfun(@getLast, timesCell);
            
            function t = getLast(times)
                if ~isempty(times)
                    t = times(end);
                else
                    t = NaN;
                end
            end
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
            p = inputParser();
            p.addOptional('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.parse(varargin{:});
            
            td.warnIfNoArgOut(nargout);
            
            if p.Results.isAligned
                % add the zero offset to the time vector for each trial
                % this is mostly for TDCA, so that alignments info is
                % preserved
                offsets = td.getTimeOffsetsFromZeroEachTrial();
                if iscell(times)
                    times = cellfun(@plus, makecol(times), num2cell(offsets), 'UniformOutput', false);
                else
                    times = makecol(times) + offsets;
                end
            end
            
            % sort and makecol
            if ~iscell(times)
                times = arrayfun(@(x) x(~isnan(x)), times, 'UniformOutput', false);
            end
            times = cellfun(@(x) makecol(sort(x)), times, 'UniformOutput', false);
            
            td = td.setChannelData(name, {times}, varargin{:});
        end
        
        function td = addOrSetEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);
            if td.hasEventChannel(name)
                td = td.setEvent(name, times, varargin{:});
            else
                td = td.addEvent(name, times, varargin{:});
            end
        end
            
    end
    
    methods % Param channel methods
        function td = addParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('name', @ischar);
            p.addOptional('values', '', @(x) true);
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.addParameter('like', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(name, varargin{:});
            
            name = p.Results.name;
            values = p.Results.values;

%             if td.hasChannel(name)
%                 warning('Overwriting existing param channel with name %s', name);
%             end

            if ~isempty(values)
                % expand scalar values to be nTrials x 1
                if ischar(values)
                    values = repmat({values}, td.nTrials, 1);
                elseif numel(values) == 1
                    values = repmat(values, td.nTrials, 1);
                end
                    
                assert(numel(values) == td.nTrials, 'Values must be vector with length %d', td.nTrials);
            end
            
            if ~isempty(p.Results.channelDescriptor)
                % use manually specified
                cd = p.Results.channelDescriptor;
            elseif ~isempty(p.Results.like)
                % copy from another channel
                cd = ParamChannelDescriptor.buildLike(td.channelDescriptorsByName.(p.Results.like), name);
            else
                % auto infer from data
                cd = ParamChannelDescriptor.buildFromValues(name, values);
            end
            cd = cd.rename(name);

            td = td.addChannel(cd, {values}, p.Unmatched);
        end
        
        function td = addScalarParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.addParameter('units', '', @ischar);
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
            if isempty(channelDescriptors)
                names = {};
                return;
            end
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
                % if this channel is marked as a scalar, convert to a double numeric array
                if ~cd.collectAsCellByField(1)
                    if cd.isBooleanByField(1)
                        values = cellfun(@logical, values);
                    else
                        values = cellfun(@double, values);
                    end
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
        
        % return nTrials x nParams cell of raw values
        function valueCell = getParamRawMultiAsCell(td, names)
            if ischar(names)
                names = {names};
            end
            
            nC = numel(names);
            valueCell = cell(td.nTrials, nC);
            for iC = 1:nC
                vals = td.getParamRaw(names{iC});
                if ~iscell(vals), vals = num2cell(vals); end
                valueCell(:, iC) = vals;
            end
        end
        
        % return nTrials x 1 cell with nParams fields
        function valueStruct = getParamRawMultiAsStruct(td, names)
            c = td.getParamRawMultiAsCell(names);
            valueStruct = cell2struct(c, names, 2);
        end
        
        % used primarily by condition info when fetching attribute data
        function paramStruct = getRawParamStruct(td)
            paramStruct = td.getParamRawMultiAsStruct(td.listParamChannels());
        end
        
        % return nTrials x nParams cell of valid values only
        function valueCell = getParamMultiAsCell(td, names)
            if ischar(names)
                names = {names};
            end
            
            nC = numel(names);
            valueCell = cell(td.nTrials, nC);
            for iC = 1:nC
                vals = td.getParam(names{iC});
                if ~iscell(vals), vals = num2cell(vals); end
                valueCell(:, iC) = vals;
            end
        end
        
        % return nTrials x 1 cell with nParams fields, valid values only
        function valueStruct = getParamMultiAsStruct(td, names)
            c = td.getParamMultiAsCell(names);
            valueStruct = cell2struct(c, names, 2);
        end
        
        function valueTable = getParamMultiAsTable(td, names, varargin)
            p = inputParser();
            p.addParameter('includeValidColumn', false, @islogical);
            p.parse(varargin{:});
            
            valueCell = td.getParamMultiAsCell(names);
            
            if p.Results.includeValidColumn
                valueCell = horzcat( num2cell(td.valid), valueCell );
                names = vertcat({'valid'}, makecol(names));
            end
            
            numWidth = ceil(log(td.nTrials)/log(10));
            trialNames = arrayfun(@(ind) sprintf('trial%0*d', numWidth, ind), (1:td.nTrials)', 'UniformOutput', false);
            
            valueTable = cell2table(valueCell, 'VariableNames', names, ...
                'RowNames', trialNames);
        end
        
        function valueStrings = getParamMultiAsStrings(td, names, varargin)
            p = inputParser();
            p.addParameter('separator', ' ', @ischar);
            p.addParameter('includeParamNames', true, @islogical);
            p.parse(varargin{:});
            
            valueStruct = td.getParamMultiAsStruct(names);
            
            valueStrings = TrialDataUtilities.Data.structArrayToStrings(valueStruct, ...
                p.Results.separator, 'includeFieldNames', p.Results.includeParamNames);
        end
        
        function values = getParamUnique(td, name)
            vals = td.getParam(name);
            vals = vals(td.valid);
%             if ~iscell(vals)
%                 vals = removenan(vals);
%             end
            values = unique(vals);
        end

        % this does no data transformations at all, just copies out fields
        % from td.data. this will not do any data class or scaling
        % conversions for you, so be careful. use getParamMultiAsStruct
        % instead
        function paramStruct = getRawChannelDataAsStruct(td, names)
            paramStruct = copyStructField(td.data, [], names);
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
        
        function values = getParamNBack(td, name, n)
            % values = getParamNBack(td, name, n)
            % retrieves the value of a parameter on the trial n trials
            % before each one. So values(j) is the value of the param on
            % trial j-n. values(1:n) will be the invalid value (NaN or
            % empty)
            valuesCurrent = td.getParam(name);
            
            cd = td.channelDescriptorsByName.(name);
            filler = cd.vectorWithMissingValue(n);
            values = cat(1, filler, valuesCurrent(1:(td.nTrials-n)));
            
            % flush values on currently invalid trials
            values = td.replaceInvalidMaskWithValue(values, cd.missingValueByField{1});
        end
        
        function td = addParamNBack(td, name, n, varargin)
            p = inputParser;
            p.addParamValue('as', '', @ischar); % new channel name
            p.parse(varargin{:});
            
            td.warnIfNoArgOut(nargout);
            %cd = td.channelDescriptorsByName.(name);
            
            if isempty(p.Results.as)
                as = [name sprintf('_%dback', n)];
            else
                as = p.Results.as;
            end
            
            values = td.getParamNBack(name, n);
            td = td.addParam(as, 'like', name, 'values', values);
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
            if isempty(channelDescriptors)
                names = {};
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'SpikeChannelDescriptor'), channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function td = addSpikeChannel(td, unitStr, varargin)
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addOptional('spikes', {}, @isvector);
            p.addParameter('isAligned', true, @isscalar);
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsTime', [], @isvector); % common time vector to be shared for ALL waveforms for this channel
            p.addParameter('waveformsField', sprintf('%s_waveforms', unitStr), @ischar);
            p.addParameter('sortQuality', NaN, @isscalar); % numeric scalar metric of sort quality
            p.addParameter('sortMethod', '', @ischar);
            p.addParameter('sortQualityEachTrial', [], @isvector);
            p.addParameter('blankingRegions', {}, @iscell); % nTrials x 1 cell of nIntervals x 2 matrices
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cd = SpikeChannelDescriptor.build(unitStr);
            cd.sortQuality = p.Results.sortQuality;
            cd.sortMethod = p.Results.sortMethod;
            
            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end
            spikes = cellfun(@plus, p.Results.spikes, num2cell(offsets), 'UniformOutput', false);

            channelData = {spikes};
            
            if ~isempty(p.Results.waveforms)
                cd = cd.addWaveformsField(p.Results.waveformsField, 'time', p.Results.waveformsTime);
                channelData{end+1} = makecol(p.Results.waveforms);
            end
            
            if ~isempty(p.Results.sortQualityEachTrial)
                assert(numel(p.Results.sortQualityEachTrial) == td.nTrials);
                cd = cd.addSortQualityEachTrialField();
                channelData{end+1} = makecol(p.Results.sortQualityEachTrial);
            end
            
            if ~isempty(p.Results.blankingRegions)
                blanking = p.Results.blankingRegions;
                assert(numel(blanking) == td.nTrials);
                assert(all(cellfun(@(x) isempty(x) || (ismatrix(x) && size(x, 2) == 2), blanking)), 'Blanking regions must be nTrials cellvec with nRegions x 2 matrices');
                
                cd = cd.addBlankingRegionsField();
                
                % align the blanking regions
                blanking = cellfun(@plus, blanking, num2cell(offsets), 'UniformOutput', false);
                
                % manually blank the spike times in the appropriate trials
                channelData{1} = TrialDataUtilities.SpikeData.blankRegionsEachTrial(channelData{1}, channelData{1}, blanking);
                
                channelData{end+1} = makecol(blanking);
            end
            
            td = td.addChannel(cd, channelData);
        end

        function td = setSpikeChannel(td, name, times, varargin)
            % td = setSpikeChannel(td, name, times, varargin)
            % options:
            %  isAligned
            %  waveforms
            p = inputParser();
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('blankingRegions', {}, @iscell); % nTrials x 1 cell of nIntervals x 2 matrices
            p.addParameter('sortQualityByTrial', [], @isvector); % nTrials x 1 vector of per-trial ratings
            p.KeepUnmatched = true; % we don't want to capture a fieldMask because this would necessitate thinking about the logic of whether to apply the blanking 
            p.parse(varargin{:});
            
            td.warnIfNoArgOut(nargout);
            
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            assert(isa(cd, 'SpikeChannelDescriptor'));

            if p.Results.isAligned
                % add the zero offset to the time vector for each trial
                % this is mostly for TDCA, so that alignments info is
                % preserved
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end
            
            if iscell(times)
                times = cellfun(@plus, makecol(times), num2cell(offsets), 'UniformOutput', false);
            else
                times = makecol(times) + offsets;
            end
            
            channelData = {times};
            channelFieldMask = true;
            
            % next field is waveforms if present
            if ~isempty(p.Results.waveforms)
                newWaves = p.Results.waveforms; % newWaves needed below
                channelData = cat(1, channelData, {newWaves});
                channelFieldMask(end+1) = true;
                
                if ~cd.hasWaveforms
                    error('Channel %s does not have waveforms, use addSpikeChannel to include waveforms', name);
                end
                
                nSpikesProvided = cellfun(@numel, times);
                nSpikesWave = cellfun(@(w) size(w, 1), newWaves);
                assert(all(isequaln(nSpikesProvided, nSpikesWave)), 'Number of spikes in each trial must match number of waveforms');
            
            else
                % no waveforms provided
                if cd.hasSpikeWaveforms
                    % no waveforms provided, but channel has waveforms
                    % check that number of spikes isn't changing so that
                    % the correspondence is maintained
                    nSpikesProvided = cellfun(@numel, times);
                    newWaves = td.getRawSpikeWaveforms(chName);
                    nSpikesWave = cellfun(@(w) size(w, 1), newWaves);
                    assert(all(isequaln(nSpikesProvided, nSpikesWave), 'Number of spikes in each trial cannot change unless waveforms are also provided'));
                    
                    channelData = cat(1, channelData, {{}});
                    channelFieldMask(end+1) = false;
                end
            end
            
            if ~isempty(p.Results.sortQualityByTrial)
                if cd.hasSortQualityEachTrial
                    channelData = cat(1, channelData, {p.Results.sortQualityByTrial});
                    channelFieldMask(end+1) = true;
                else
                    error('Channel %s does not have sortQualityByTrial, use addSpikeChannel to add these', name);
                end
            else
                if cd.hasSortQualityEachTrial
                    channelData = cat(1, channelData, {{}});
                    channelFieldMask(end+1) = false;
                end
            end
                
            if ~isempty(p.Results.blankingRegions)
                blanking = p.Results.blankingRegions;
                
                assert(numel(blanking) == td.nTrials);
                %assert(all(cellfun(@(x) isempty(x) || (ismatrix(x) && size(x, 2) == 2), blanking)), 'Blanking regions must be nTrials cellvec with nRegions x 2 matrices');
                
                % align the blanking regions with the zero offset to
                % unalign the blanking regions specified
                blanking = cellfun(@plus, blanking, num2cell(offsets), 'UniformOutput', false);
                
                % add field if not found
                if ~cd.hasBlankingRegions
                    cd = cd.addBlankingRegionsField();
                    td.channelDescriptorsByName.(name) = cd;
                    
                    blanking = TrialDataUtilities.SpikeData.removeOverlappingIntervals(blanking);
                else
                    % need to merge with existing blanking regions
                    blankingExisting = td.getSpikeBlankingRegions(name);
                    % these will be aligned, so we unalign them to match
                    % new blanking
                    alignedOffsets = td.getTimeOffsetsFromZeroEachTrial();
                    
                    blankingExisting = cellfun(@plus, blankingExisting, num2cell(alignedOffsets), 'UniformOutput', false);
                    
                    % combine new with existing now that they are both
                    % unaligned
                    blanking = TrialDataUtilities.SpikeData.removeOverlappingIntervals(blanking, blankingExisting);
                end
                
            elseif td.channelDescriptorsByName.(name).hasBlankingRegions
                % need to blank new times with old blanking regions
                blanking = td.getSpikeBlankingRegions(name);
                % these will be aligned, so we unalign them to match
                % new blanking
                alignedOffsets = td.getTimeOffsetsFromZeroEachTrial();

                blanking = cellfun(@plus, blanking, num2cell(alignedOffsets), 'UniformOutput', false);
            else
                blanking = {};
            end
                
            if ~isempty(blanking)
                % manually blank the spike times in the appropriate trials
                oldSpikeTimes = channelData{1};
                channelData{1} = TrialDataUtilities.SpikeData.blankRegionsEachTrial(channelData{1}, oldSpikeTimes, blanking);
                
                if cd.hasWaveforms
                    % manually blank the spike waveforms too
                    % be sure to use the oldSpikeTimes to maintain the
                    % correspondence
                    newWaves = TrialDataUtilities.SpikeData.blankRegionsEachTrial(newWaves, oldSpikeTimes, blanking);
                    
                    if ~isempty(p.Results.waveforms)
                        % waveforms provided, replace channelData{2}
                        channelData{2} = newWaves;
                    else
                        % no waveforms provided, now we need to add them to
                        % channelData so they get blanked
                        channelData = cat(1, channelData(1), {newWaves}, channelData(2:end));
                        channelFieldMask = cat(1, channelFieldMask(1), true, makecol(channelFieldMask(2:end)));
                    end
                end
                
                % and add the blanking info as the last channel data
                channelData = cat(1, channelData, {blanking});
                channelFieldMask(end+1) = true;
            else
                if cd.hasBlankingRegions
                    channelData = cat(1, channelData, {{}});
                    channelFieldMask(end+1) = false;
                end
            end

            % make the final update
            td = td.setChannelData(name, channelData, 'fieldMask', channelFieldMask, p.Unmatched);
        end
        
        function td = mergeSpikeChannels(td, chList, varargin)
            assert(iscell(chList));
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addParameter('as', chList{1}, @ischar);
            p.parse(varargin{:});
            
            % gather all spike times and waveforms
            nCh = numel(chList);
            [wavesCell, timesCell] = deal(cell(nCh, td.nTrials));
            waveTvecCell = cell(nCh, 1);
            [idxZero, nTimeWave, deltaTimeWave] = nanvec(nCh);
            blankingRegionsCell = cell(nCh, 1);
            
            for c = 1:nCh
                [wavesCell(c, :), waveTvecCell{c}, timesCell(c, :)] = ...
                    td.getRawSpikeWaveforms(chList{c});
                [~, idxZero(c)] = min(abs(waveTvecCell{c}));
                nTimeWave(c) = numel(waveTvecCell{c});
                deltaTimeWave(c) = mean(diff(waveTvecCell{c}));
                blankingRegionsCell{c} = td.getRawSpikeBlankingRegions(chList{c});
            end
            
            if max(deltaTimeWave) - min(deltaTimeWave) > median(deltaTimeWave) * 0.05
                warning('Spike channel waveforms have different sampling rates. Using median');
            end
            deltaTimeWave = median(deltaTimeWave);
            
            % figure out how to match the waveform time vectors
            zeroInd = max(idxZero);
            nRightZeroGlobal = max(nTimeWave - idxZero);
            waveTvecGlobal = -deltaTimeWave*(zeroInd-1) : deltaTimeWave : deltaTimeWave*nRightZeroGlobal;
            
            % pad the waveforms to the left or right with NaNs to match sizes
            for c = 1:nCh
                padLeft = zeroInd - idxZero(c);
                padRight = nRightZeroGlobal - (nTimeWave(c) - idxZero(c));
                if padLeft > 0 || padRight > 0
                    debug('Padding waveforms for %s to match new combined waveforms time vector\n', chList{c});          
                    for t = 1:size(wavesCell, 2)
                        nWaves = size(wavesCell{c, t}, 1);
                        wavesCell{c, t} = cat(2, zeros(nWaves, padLeft), wavesCell{c, t}, zeros(nWaves, padRight));
                    end
                end
            end
            
            % combine the data across channels
            [wave, times] = cellvec(td.nTrials);
            prog = ProgressBar(td.nTrials, 'Combining spike data');
            for iT = 1:td.nTrials
                prog.update(iT);
                wave{iT} = cat(1, wavesCell{:, iT});
                times{iT} = cat(1, timesCell{:, iT});  
            end
            prog.finish();
            
            % combine blanking regions, add spike channel will take care of
            % blanking the data in times and waves for us
            blankingRegions = TrialDataUtilities.SpikeData.removeOverlappingIntervals(blankingRegionsCell{:});
            
            % add the new channel with combined data
            td = td.dropChannels(chList);
            td = td.addSpikeChannel(p.Results.as, times, 'waveforms', wave, ...
                'waveformsTime', waveTvecGlobal, ...
                'isAligned', false, ...
                'blankingRegions', blankingRegions);
        end
        
        function td = blankSpikesInTimeIntervals(td, name, intervalCell, varargin)
            % adds a blanking region to the spiking data
            % this both removes the spikes from that period of time each
            % trial, and will inform the spike rate filtering that this
            % time interval is not observed, rather than simply has no
            % spikes
            %
            % intervalCell is nTrials x 1 cell with nIntervals(iT) x 2 matrices inside 
            
            td.warnIfNoArgOut(nargout);
            
            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            assert(isa(cd, 'SpikeChannelDescriptor'));
            
            p = inputParser();
            p.addParameter('isAligned', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            % need to request all spike times
            [rawWaves, ~, rawTimes] = td.getRawSpikeWaveforms(name);
            
            if ~iscell(intervalCell)
                % assume same for all trials
                assert(size(intervalCell, 2) == 2);
                intervalCell = repmat({intervalCell}, td.nTrials, 1);
            else
                nonEmpty = ~cellfun(@isempty, intervalCell);
                nCols = cellfun(@(x) size(x, 2), intervalCell(nonEmpty));
                assert(all(nCols == 2), 'Interval cell contents must contain matrices of intervals with 2 columns');
            end
            assert(iscell(intervalCell));
            intervalCell = makecol(intervalCell);
            
            % then we need to unalign the intervals to match the raw times
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
                intervalCell = cellfun(@plus, intervalCell, num2cell(offsets), 'UniformOutput', false);
            end
            
            % setSpikeChannel will take care of blanking the waveforms
            % be sure we specify that the data are no longer aligned
            td = td.setSpikeChannel(name, rawTimes, 'isAligned', false, ...
                'blankingRegions', intervalCell, 'waveforms', rawWaves);
        end
        
        function intervalCell = getRawSpikeBlankingRegions(td, unitName)
            cd = td.channelDescriptorsByName.(unitName);
            fld = cd.blankingRegionsField;
            if isempty(fld)
                % no blanking regions, just return blank
                intervalCell = cellvec(td.nTrials);
            else
                intervalCell = makecol({td.data.(fld)});
                intervalCell = TrialDataUtilities.SpikeData.removeOverlappingIntervals(intervalCell);
            end
        end
        
        function intervalCell = getSpikeBlankingRegions(td, unitName)
            intervalCell = td.getRawSpikeBlankingRegions(unitName);
            intervalCell = td.replaceInvalidMaskWithValue(intervalCell, []);
        end
        
        function timesCell = getRawSpikeTimes(td, unitNames)
            if ischar(unitNames)
                timesCell = {td.data.(unitNames)}';
            elseif iscellstr(unitNames)
                nUnits = numel(unitNames);
                timesCellByUnit = cell(td.nTrials, nUnits);
                for iU = 1:nUnits
                    timesCellByUnit(:, iU) = {td.data.(unitNames{iU})}';
                end
                timesCell = cellvec(td.nTrials);
                for iT = 1:td.nTrials
                    timesCell{iT} = cat(1, timesCellByUnit{iT, :});
                end
            else
                error('Unsupported unit name argument');
            end
        end

        function timesCell = getSpikeTimes(td, unitNames, varargin) 
            timesCell = td.getSpikeTimesUnaligned(unitNames);
        end
        
        function timesCell = getSpikeTimesUnaligned(td, unitNames)
            timesCell = td.getRawSpikeTimes(unitNames);
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
        end
            
        function counts = getSpikeCounts(td, unitName)
            counts = cellfun(@numel, td.getSpikeTimes(unitName));
            counts = td.replaceInvalidMaskWithValue(counts, NaN);
        end
        
        function rates = getSpikeMeanRate(td, unitName)
            counts = td.getSpikeCounts(unitName);
            durations = td.getValidDurationsForSpikeChannel(unitName);
            rates = counts ./ durations * td.timeUnitsPerSecond;
        end
        
        function [tf, blankingRegionsField] = hasBlankingRegions(td, unitName)
            if ~td.hasSpikeChannel(unitName)
                tf = false;
                return;
            end
            tf = td.channelDescriptorsByName.(unitName).hasBlankingRegions;
            blankingRegionsField = td.channelDescriptorsByName.(unitName).blankingRegionsField;
        end
        
        function tf = hasSpikeWaveforms(td, unitName)
            if ~td.hasSpikeChannel(unitName)
                tf = false;
                return;
            end
            wavefield = td.channelDescriptorsByName.(unitName).waveformsField;
            tf = ~isempty(wavefield);
        end
        
        function [wavesCell, waveTvec, timesCell] = getRawSpikeWaveforms(td, unitName)
            cd = td.channelDescriptorsByName.(unitName);
            wavefield = cd.waveformsField;
            assert(~isempty(wavefield), 'Unit %s does not have waveforms', unitName);
            wavesCell = {td.data.(wavefield)}';
            % scale to appropriate units
            wavesCell = cd.scaleWaveforms(wavesCell);
            waveTvec = td.channelDescriptorsByName.(unitName).waveformsTime;
            
            % check number of timepoints 
            waveMat = TrialDataUtilities.Data.getFirstNonEmptyCellContents(wavesCell);
            if ~isempty(waveMat)
                nSampleWave = size(waveMat, 2);
                if nSampleWave < numel(waveTvec)
                    warning('Waveforms have %d samples but waveformsTime has %d samples. Shortening waveforms to match', nSampleWave, numel(waveTvec));
                    waveTvec = waveTvec(1:nSampleWave);
                elseif nSampleWave > numel(waveTvec)
                    error('Waveforms have %d samples but waveformsTime has %d samples. Provide new waveform time vector', nSampleWave, numel(waveTvec));
                end
            end
            timesCell = td.getRawSpikeTimes(unitName);
        end
        
        function [wavesCell, waveTvec, timesCell] = getSpikeWaveforms(td, unitName)
            [wavesCell, waveTvec, timesCell] = td.getRawSpikeWaveforms(unitName);
            wavesCell = td.replaceInvalidMaskWithValue(wavesCell, []);
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
        end
    end
    
    methods % spike / continuous neural channel correspondence
        function info = listArrayElectrodesWithSpikeChannelsAsTable(td)
            % info is struct with fields array (char) and electrode
            % (numeric)
            names = td.listSpikeChannels();
            cdCell = td.getChannelDescriptorMulti(names);
            
            array = cellfun(@(cd) cd.array, cdCell, 'UniformOutput', false);
            electrode = cellfun(@(cd) cd.electrode, cdCell);
            
            t = table(array, electrode);
            [info, ~, idx] = unique(t);
            
            channels = cellfun(@(cd) cd.name, cdCell, 'UniformOutput', false);
            
            chMatch = cellvec(height(info));
            for r = 1:height(info)
                chMatch{r} = channels(idx == r);
            end
            info.channelList = chMatch;
        end
        
        function [names, units] = listSpikeChannelsOnArray(td, arrayName)
            % [names, channelDescriptors] = getChannelsOnArray(td, arrayName)
            names = td.listSpikeChannels();
            units = nanvec(numel(names));
            mask = falsevec(numel(names));
            for iC = 1:numel(names)
                cd = td.channelDescriptorsByName.(names{iC});
                mask(iC) = strcmp(cd.array, arrayName);
                units(iC) = cd.unit;
            end
                
            names = names(mask);
            units = units(mask);
        end
        
        function [names] = listContinuousNeuralChannelsOnArray(td, arrayName)
        % [names, channelDescriptors] = getContinuousNeuralChannelsOnArray(td, arrayName)
            names = td.listContinuousNeuralChannels();
            mask = falsevec(numel(names));
            for iC = 1:numel(names)
                cd = td.channelDescriptorsByName.(names{iC});
                mask(iC) = strcmp(cd.array, arrayName);
            end
                
            names = names(mask);
        end
        
        function [names, units] = listSpikeChannelsOnArrayElectrode(td, arrayName, electrodeNum, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroUnit', false, @islogical);
            p.parse(varargin{:});
            
            % [names, channelDescriptors] = getChannelsOnArrayElectrode(td, arrayName, electrodeNum)
            names = td.listSpikeChannels();
            units = nanvec(numel(names));
            mask = falsevec(numel(names));
            for iC = 1:numel(names)
                cd = td.channelDescriptorsByName.(names{iC});
                mask(iC) = any(strcmp(cd.array, arrayName)) && any(cd.electrode == electrodeNum);
                units(iC) = cd.unit;
            end
                
            if p.Results.ignoreZeroUnit
                mask = mask & units ~= 0;
            end
            names = names(mask);
            units = units(mask);
        end
        
        function [names] = listContinuousNeuralChannelsOnArrayElectrode(td, arrayName, electrodeNum)
            % [names, channelDescriptors] = getContinuousNeuralChannelsOnArrayElectrode(td, arrayName, electrodeNum)
            names = td.listContinuousNeuralChannels();
            mask = falsevec(numel(names));
            for iC = 1:numel(names)
                cd = td.channelDescriptorsByName.(names{iC});
                mask(iC) = strcmp(cd.array, arrayName) && cd.electrode == electrodeNum;
            end
                
            names = names(mask);
        end
        
        function [names, units] = listSpikeChannelsOnSameArrayElectrodeAs(td, chName, varargin)
            % [names, units] = listSpikeChannelsOnSameArrayElectrodeAs(td, chName)
            % chName is spike channel or continuous neural channel name
            cd = td.channelDescriptorsByName.(chName);
            [names, units] = td.listSpikeChannelsOnArrayElectrode(cd.array, cd.electrode, varargin{:});
        end
        
         function [names] = listContinuousNeuralChannelsOnSameArrayElectrodeAs(td, chName)
             % [names, units] = listContinuousNeuralChannelsOnSameArrayElectrodeAs(td, chName)
             % chName is spike channel or continuous neural channel name
            cd = td.channelDescriptorsByName.(chName);
            [names] = td.listContinuousNeuralChannelsOnArrayElectrode(cd.array, cd.electrode);
         end
    end

    methods % Generic add data methods
        function td = updatePostDataChange(td, fieldsAffected) %#ok<INUSD>
            % call this after any changes to td.data.(fieldsAffected)
            td.warnIfNoArgOut(nargout);
        end
        
        function td = updatePostChannelDataChange(td, chName)
            td.warnIfNoArgOut(nargout);
            fields = td.channelDescriptorsByName.(chName).dataFields;
            td = td.updatePostDataChange(fields);
        end
        
        function offsets = getTimeOffsetsFromZeroEachTrial(td)
            % when adding new data to the trial, all times are stored relative
            % to the current time zero. This will be overridden in 
            % TrialDataConditionAlign. This will be added automatically
            % to all new channel time data to match the offsets produced when
            % getting data
            offsets = zerosvec(td.nTrials);
        end
        
        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrialRaw(td)
            % overriden by TDCA to be zero-relative
            tMinByTrial = td.getEvent('TrialStart');
            tMaxByTrial = td.getEvent('TrialStop');
        end
        
        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrial(td)
            [tMinByTrial, tMaxByTrial] = td.getTimeStartStopEachTrialRaw();
            tMinByTrial(~td.valid) = NaN;
            tMaxByTrial(~td.valid) = NaN;
        end
        
        function tf = alignIncludesFullTrial(td)
            tf = true;
        end 
        
        function td = addChannels(td, cds, varargin)
            td.warnIfNoArgOut(nargout);
            
            % handle struct with .name = cd, or array of ChannelDescriptors
            if isstruct(cds)
                flds = fieldnames(cds);
                getFn = @(i) cds.(flds{i});
                nChannels = numel(flds);
            else
                getFn = @(i) cds(i);
                nChannels = numel(cds);
            end
            
            for iC = 1:nChannels
            	cd = getFn(iC);
                td = td.addChannel(cd, varargin{:});
            end
        end
        
        function td = addChannelsMissing(td, cds)
            td.warnIfNoArgOut(nargout);
            if iscell(cds)
                has = cellfun(@(cd) td.hasChannel(cd.name), cds);
                cds = cds(~has);
            elseif isstruct(cds)
                flds = fieldnames(cds);
                has = structfun(@(cd) td.hasChannel(cd.name), cds);
                cds = rmfield(cds, flds(has));
            else
                error('Must be cell array or struct of ChannelDescriptors');
            end
            if any(~has)
                td = td.addChannels(cds);
            end
        end
        
        function td = addChannel(td, cd, varargin)
            % adds a new channel described by ChannelDescriptor cd.
            % valueCell must be nDataFields x 1 cell each with nTrials x 1
            % cell within. Note that if any trials are marked invalid, this
            % data will be cleared when adding the channel.
            p = inputParser();
            p.addOptional('valueCell', {}, @(x) true);
            p.addParameter('clearIfPresent', false, @islogical);
           % p.addParameter('ignoreOverwriteChannel', false, @islogical);
            p.addParameter('updateValidOnly', true, @islogical);
            p.parse(varargin{:});
            valueCell = makecol(p.Results.valueCell);
            
            td.warnIfNoArgOut(nargout);
            
            % check for overwrite if requested
            assert(isa(cd, 'ChannelDescriptor'), 'Argument cd must be ChannelDescriptor');  
            
            % give the channelDescriptor a chance to initialize itself
            cd = cd.initialize();
            
            alreadyHasChannel = td.hasChannel(cd.name);
            if alreadyHasChannel
                warning('Overwriting existing channel with name %s', cd.name);
                td = td.dropChannels(cd.name);
                
%                 % check that existing channel matches
%                 assert(isequaln(cd, td.channelDescriptorsByName.(cd.name)), ...
%                     'ChannelDescriptor for channel %s does not match existing channel', cd.name);
            end

            td.channelDescriptorsByName.(cd.name) = cd;
            
            % touch each of the value fields to make sure they exist
            for iF = 1:cd.nFields
                if ~isfield(td.data, cd.dataFields{iF});
                    td.data(end).(cd.dataFields{iF}) = [];
                end
            end     
            
            if isempty(valueCell) && (~alreadyHasChannel || p.Results.clearIfPresent)
                td = td.clearChannelData(cd.name, 'fieldMask', ~cd.isShareableByField);
            elseif ~isempty(valueCell)
                % clear on fields where no values provided and it's not shared, 
                % set on fields where values are provided
                nonEmptyMask = ~cellfun(@isempty, valueCell);
                td = td.clearChannelData(cd.name, 'fieldMask', ~nonEmptyMask & ~cd.isShareableByField);
                td = td.setChannelData(cd.name, valueCell, 'fieldMask', nonEmptyMask, ...
                    'updateValidOnly', p.Results.updateValidOnly);
            end
        end
        
        function td = copyChannel(td, oldName, newName)
            % copy channel oldName to newName
            td.warnIfNoArgOut(nargout); 
            
            alreadyHasChannel = td.hasChannel(newName);
            if alreadyHasChannel
                warning('Overwriting existing channel with name %s', newName);
                td = td.dropChannels(newName);
            end

            % ask the channel descriptor to rename itself and store the copy
            [cd, dataFieldRenameMap] = td.channelDescriptorsByName.(oldName).rename(newName);            
            td.channelDescriptorsByName.(newName) = cd;
            
            % then copy changed channel fields
            flds = fieldnames(dataFieldRenameMap);
            for iF = 1:numel(flds)
                oldField = flds{iF};
                newField = dataFieldRenameMap.(oldField);
                td = td.renameDataField(oldField, newField, oldName, true); % rename but make a copy
            end
            
            % if analog channel in group, separate it from the group
            % this will copy the data, but it won't get rid of the old data
            % since the original channel still references that column
            if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix
                td = td.separateAnalogChannelFromGroup(newName);
            end
                  
            td = td.updatePostChannelDataChange(newName);
        end
        
        function td = renameChannel(td, oldName, newName)
            % rename channel name to newName
            % if channel name's primary data field is also name, rename
            % that field too to newName
            td.warnIfNoArgOut(nargout); 

            % update channel descriptor directly
            [cd, dataFieldRenameMap] = td.channelDescriptorsByName.(oldName).rename(newName);            
            td.channelDescriptorsByName = rmfield(td.channelDescriptorsByName, oldName);
            td.channelDescriptorsByName.(newName) = cd;
            
            % then rename channel fields
            flds = fieldnames(dataFieldRenameMap);
            for iF = 1:numel(flds)
                oldField = flds{iF};
                newField = dataFieldRenameMap.(oldField);
                td = td.renameDataField(oldField, newField, oldName);
            end
        end
        
        function td = renameDataField(td, field, newFieldName, ignoreChannelList, copy)
           % rename field number fieldIdx of channel name to newFieldName
           % if that field is shared by a channel not in ignoreChannelList, make a copy
           td.warnIfNoArgOut(nargout); 
           
           if nargin < 4
               ignoreChannelList = {};
           end
           
           if nargin < 5
               copy = false;
           end

           otherChannels = setdiff(td.getChannelsReferencingFields(field), ignoreChannelList);
           if ~isempty(otherChannels) || copy
               % being used by other channels, make a copy
               td.data = copyStructField(td.data, td.data, field, newFieldName);
           else
               % move field and delete old one
               td.data = mvfield(td.data, field, newFieldName);
           end
        end
        
        function td = copyRenameSharedChannelFields(td, names, fieldMask)
            % copy all fields which belong to this channel which are shared
            % with other channes, rename those fields, and mark the
            % updated fields inside that fields channelDescriptor. if name
            % is a cellstr of channel names, only consider fields shared
            % with channels not in that set of channels.
           td.warnIfNoArgOut(nargout); 
           
           if ischar(names)
               names = {names};
           end
           cdCell = td.getChannelDescriptorMulti(names);
           
           for iC = 1:numel(names)
               cd = cdCell{iC};
               
               if nargin < 3 || isempty(fieldMask)
                   fieldMaskThis = true(cd.nFields, 1);
               else
                   fieldMaskThis = TensorUtils.vectorIndicesToMask(fieldMask, cd.nFields);
               end
               
               isShareable = cd.isShareableByField; 
               dataFields = cd.dataFields;
               fieldsRenamed = struct();
               for iF = 1:cd.nFields
                   if ~isShareable(iF) || ~fieldMaskThis(iF), continue; end
                   if ismember(dataFields{iF}, fieldnames(fieldsRenamed))
                       % already renamed earlier, change cd to reflect that
                       newName = fieldsRenamed.(dataFields{iF});
                       cd = cd.renameDataField(iF, newName);
                   end
                       

                   fullList = td.getChannelsReferencingFields(dataFields{iF});
                   otherChannels = setdiff(fullList, names);
                   if isempty(otherChannels), continue; end

                   % being used by other channels outside of names, rename and copy
                   % use suggestion from ChannelDescriptor first, then
                   % uniquify
                   template = cd.suggestFieldName(iF);
                   template = template(1:min(numel(template), namelengthmax()-10));
                   newName = matlab.lang.makeUniqueStrings(template, fieldnames(td.data));
                   td.data = copyStructField(td.data, td.data, dataFields{iF}, newName);

                   fieldsRenamed.(dataFields{iF}) = newName;
                   cd = cd.renameDataField(iF, newName);
               end
               
               td.channelDescriptorsByName.(cd.name) = cd;
           end
           
        end
        
        function td = clearChannelData(td, name, varargin)
            p = inputParser();
            p.addParameter('fieldMask', [], @islogical);
            p.parse(varargin{:});
            
            cd = td.channelDescriptorsByName.(name);
            fieldMask = p.Results.fieldMask;
            if isempty(fieldMask)
                fieldMask = true(cd.nFields, 1);
            end
            if ~any(p.Results.fieldMask)
                return;
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
                td.data = assignIntoStructArray(td.data, dataFields{iF}, val);
%                 for iT = 1:td.nTrials
%                     td.data(iT).() = val;
%                 end
            end
            
        end
        
        function [td, trialsUpdated, fieldsUpdated] = setChannelData(td, name, valueCell, varargin)
            % note that by default, updateValidOnly is true, meaning that
            % the values on invalid trials will not be updated
            % trialsUpdated is a mask indicating which trials were changed
            % (or cleared)
            % fieldsUpdated is a mask indicating which of the data fields
            % were updated as well
            p = inputParser();
            p.addParameter('fieldMask', [], @islogical);
            p.addParameter('clearForInvalid', false, @islogical);
            p.addParameter('updateValidOnly', true, @islogical);
            p.addParameter('updateMask', [], @isvector);
            p.parse(varargin{:});
            
            updateMaskManual = p.Results.updateMask;
            if isempty(updateMaskManual)
                updateMaskManual = true(td.nTrials, 1);
            end
            assert(numel(updateMaskManual) == td.nTrials, ...
                'Size of updateMask must match nTrials');
            
            td.warnIfNoArgOut(nargout);
            
            cd = td.channelDescriptorsByName.(name);
            fieldMask = p.Results.fieldMask;
            if isempty(fieldMask)
                fieldMask = false(cd.nFields, 1);
                fieldMask(1:numel(valueCell)) = true; % assume that only the first few fields are specified
            else
                assert(numel(fieldMask) == cd.nFields);
            end
            fieldMask = makecol(fieldMask);
            
            % expand valueCell with missing data
            if numel(valueCell) < cd.nFields
                for iF = 1:cd.nFields
                    if numel(valueCell) < iF
                        if ~fieldMask(iF)
                            % fill in with empty
                            valueCell{iF} = cell(td.nTrials, 1);
                        else
                            error('valueCell missing field %d', iF);
                        end
                    end
                end
            end
            
            dataFields = cd.dataFields;
            nFields = numel(dataFields);
            
            % check that one value list was provided for each data field
            % referenced by the ChannelDescriptor
            assert(numel(valueCell) == nFields, ...
                'Channel Descriptor references %d fields but only %d field value lists provided', ...
                numel(nFields), numel(valueCell));
            
            % check for fields that don't exist and clear them so that they
            % do exist
            fieldMissing = makecol(~ismember(dataFields, fieldnames(td.data)));
            td = td.clearChannelData(cd.name, 'fieldMask', fieldMissing & fieldMask);
            
            % avoid overwriting data shared by other channels
            td = td.copyRenameSharedChannelFields(cd.name, fieldMask);
            cd = td.channelDescriptorsByName.(name);
            dataFields = cd.dataFields;      
            
            updateMask = updateMaskManual;
            if p.Results.updateValidOnly
                updateMask = updateMask & td.valid;
            end
                    
            for iF = 1:nFields
                % only touch specified fields
                if ~fieldMask(iF), continue; end
                
                if p.Results.clearForInvalid
                    % here we want the update mask to stay the same as
                    % updateMaskManual so that everything gets updated with
                    % the cleared value
                    valueCell{iF} = td.replaceInvalidMaskWithValue(valueCell{iF}, cd.missingValueByField{iF});
                end
               
%                 for iT = 1:numel(td.data)
%                     if updateMask(iT)
%                         td.data(iT).(dataFields{iF}) = valueCell{iF}{iT};
%                     end
%                 end
                if iF == 1 && isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix
                    colIdx = cd.primaryDataFieldColumnIndex;
                    % shared analog channel, just replace that column
                    prog = ProgressBar(td.nTrials, 'Updating %s (%s) data', name, dataFields{iF});
                    for iT = 1:td.nTrials
                        prog.update(iT);
                        if updateMask(iT)
                            if td.valid(iT)
                                td.data(iT).(dataFields{iF})(:, colIdx) = valueCell{iF}{iT};
                            else
                                td.data(iT).(dataFields{iF})(:, colIdx) = NaN;
                            end
                        end
                    end
                    prog.finish();
                else
                    % normal field
                    td.data = assignIntoStructArray(td.data, dataFields{iF}, valueCell{iF}(updateMask, :), updateMask);
                end
            end 
            
            % convert data, also give cd a chance to update its memory
            % storage class for the data to reflect what was passed in
            [td.data, cd] = cd.repairData(td.data);

            td.channelDescriptorsByName.(cd.name) = cd;
            td = td.updatePostDataChange(dataFields(fieldMask));
            
            fieldsUpdated = fieldMask;
            trialsUpdated = updateMask;
        end
    end

    methods % Plotting functions
        function str = getAxisLabelForChannel(td, name)
            str = td.channelDescriptorsByName.(name).getAxisLabelPrimary();
        end

        function str = getTimeAxisLabel(td)
           str = sprintf('Time (%s)',  td.timeUnitName);
        end
        
        function plotAnalogRawEachTrial(td, name, varargin) 
            p = inputParser();
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            [dataCell, timeCell] = td.getAnalogRaw(name);     

            dataCell = dataCell(td.valid);
            timeCell = timeCell(td.valid);

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
            
%             AutoAxis.replace(axh);
            
            hold(axh, 'off');
        end
    end
    
    
    % save / load wrappers for CacheCustomSaveLoad
    properties
        cacheWithSaveFast = false;
    end
    
    methods
        function tf = getUseCustomSaveLoad(td, info) %#ok<INUSD>
            tf = td.cacheWithSaveFast;
        end
        
        function token = saveCustomToLocation(td, location)
            td.saveFast(location);
            token = [];
        end
    end
    
    methods(Static)
        function data = loadCustomFromLocation(location, token) %#ok<INUSD>
            data = TrialData.loadFast(location);
        end
    end
end

