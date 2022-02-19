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

        datasetMeta = struct(); % arbitrary, user determined data, use setMetaKey and getMetaKey to access

        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = ''; % how did we access the original data?

        timeUnitName = 'ms';

        timeUnitsPerSecond = 1000;

        channelDescriptorsByName  = struct(); % struct with ChannelDescriptor for each channel, by name

        trialDescriptionExtraParams = {} % used in generating descriptions of each trial
    end

    % properties that scale with the number of trials, should be referenced
    % in all select and merge methods
    properties(SetAccess=protected)
        data = struct([]);  % standardized format nTrials x 1 struct with all trial data
        manualValid = true(0, 1);
        manualInvalidCause = {};
        temporaryValid = [];
        temporaryInvalidCause = {};
    end

    properties(Access=protected, Hidden)
        odc % TrialDataOnDemandCache instance

        nTrialsManual % used when loading just the metadata, otherwise ignored
    end

    % properties which are stored in odc, see get/set below
    properties(Dependent)
        valid
        invalidCause % cellstr of explanations

        permanentlyInvalid
    end

    % Convenience dependent properties
    properties(Dependent)
        ch % a struct which contains all channel names as fields, whose values are the channel names. This is a hack to facilitate tab completion instead of strings.
        nTrials
        nTrialsValid

        nTrialsPermanentlyInvalid

        nChannels
        channelNames % cell array of channel names

        timeUnitsPerMs
    end

    % Initializing and building
    methods
        function td = TrialData(varargin)
            p = inputParser();
            p.addOptional('buildFrom', [], @(x) isa(x, 'TrialData') || isa(x, 'TrialDataInterface'));
            %             p.addParamValue('quiet', false, @islogical); % completely silent output
            p.addParameter('suppressWarnings', false, @islogical); % don't warn about any minor issues
            p.parse(varargin{:});
            %             quiet = p.Results.quiet;

            td = td.rebuildOnDemandCache();

            if isempty(p.Results.buildFrom)
                from = ManualTrialDataInterface();
                suppressWarnings = true;
            else
                from = p.Results.buildFrom;
                suppressWarnings = p.Results.suppressWarnings;
            end

            if isa(from, 'TrialData')
                td = td.initializeFromTrialData(from, ...
                    'suppressWarnings', suppressWarnings);
            elseif isa(from, 'TrialDataInterface')
                td = td.initializeFromTrialDataInterface(from, ...
                    'suppressWarnings', suppressWarnings);
            else
                error('Unknown initializer');
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
            p.addParameter('suppressWarnings', false, @islogical);
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
            specialNames = string({specialParams.name}');
            regularChannels = makecol(tdi.getChannelDescriptors('suppressWarnings', p.Results.suppressWarnings));
            if isempty(regularChannels)
                regularNames = string([]);
            else
                regularNames = string({regularChannels.name}');
            end

            % check for reserved channel names
            overlap = intersect(specialNames, regularNames);
            if ~isempty(overlap)
                error('getChannelDescriptors() returned the following reserved channel names: %s', TrialDataUtilities.String.strjoin(overlap, ', '));
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
            [data, td.channelDescriptorsByName] = td.validateDataInternal(data, td.channelDescriptorsByName, 'suppressWarnings', p.Results.suppressWarnings);
            td.data = makecol(data);

            td.manualValid = truevec(td.nTrials);

            % check that all analog data matches the number of time samples
            % and clear trials where they don't. Warnings will be shown for
            % violating channels
            td = td.fixCheckAnalogDataMatchesTimeVectors();
            td.initialized = true;
        end

        function td = addNewTrialsFromTrialDataInterface(td, varargin)
            % this is a reduced form of initializeFromTDI that requests only new channel descriptors and new data
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addRequired('trialDataInterface', @(tdi) isa(tdi, 'TrialDataInterface'));
            p.addParameter('suppressWarnings', false, @islogical);
            p.parse(varargin{:});
            tdi = p.Results.trialDataInterface;

            % check equivalence of basic details from the TrialDataInterface
            assert(strcmp(td.trialDataInterfaceClass, class(tdi)));
            assert(strcmp(td.timeUnitName, tdi.getTimeUnitName()));
            assert(td.timeUnitsPerSecond == tdi.getTimeUnitsPerSecond());

            % request channel descriptors for both special params and regular channels
            newChannels = makecol(tdi.getNewChannelDescriptors());
            maskRemove = ismember({newChannels.name}, td.listChannels('includeNamedSubChannels', false));
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
            prog = ProgressBar(numel(trialData), 'Adding trials from TrialData instances');
            for i = 1:numel(trialData)
                prog.update(i);
                td = td.addNewTrialsRaw(trialData{i}, varargin{:});
            end
            prog.finish();
        end

        function td = addNewTrialsRaw(td, dataOrTrialData, varargin)
            td.warnIfNoArgOut(nargout);
            p = inputParser();
            p.addParameter('suppressWarnings', false, @islogical);
            p.parse(varargin{:});

            %             data = struct([]);  % standardized format nTrials x 1 struct with all trial data
            %             manualValid = true(0, 1);
            %             manualInvalidCause = {};
            %             temporaryValid = [];
            %             temporaryInvalidCause = {};

            concatData = td.data;
            if isa(dataOrTrialData, 'TrialData')
                tdNew = dataOrTrialData;
                newData = tdNew.data;
            else
                tdNew = [];
                newData = dataOrTrialData;
            end

            % validate new data against all channel descriptors (old + new)
            %debug('Validating new channel data...\n');
            newData = td.validateDataInternal(newData, td.channelDescriptorsByName, 'addMissingFields', true, ...
                'suppressWarnings', p.Results.suppressWarnings);

            % concatenate onto the old data
            td.data = TrialDataUtilities.Data.structcat(1, concatData, newData);

            % concatenate the manual valid array
            if isempty(tdNew)
                newValid = truevec(numel(newData));
                newCause = cellvec(numel(newData));
            else
                newValid = tdNew.manualValid;
                newCause = tdNew.manualInvalidCause;
                assert(numel(newValid) == numel(newData), 'manualValid must be same length as data');
            end
            td.manualValid = cat(1,  td.manualValid, newValid);
            td.manualInvalidCause = cat(1, td.manualInvalidCause, newCause);

            % concatenate the temporary valid array
            if isempty(tdNew)
                newValid = truevec(numel(newData));
                newCause = cellvec(numel(newData));
            else
                newValid = tdNew.temporaryValid;
                newCause = tdNew.temporaryInvalidCause;
                assert(numel(newValid) == numel(newData), 'manualValid must be same length as data');
            end
            td.temporaryValid = cat(1,  td.temporaryValid, newValid);
            td.temporaryInvalidCause = cat(1, td.temporaryInvalidCause, newCause);

            td = td.postAddNewTrials(); % update valid here, also to allow TDCA to update itself before update valid since TDCA's override executes first
        end

        function td = postAddNewTrials(td)
            td.warnIfNoArgOut(nargout);

            % copy and flush valid
            td.odc = td.odc.copy();
            td.odc.flush();

            td = td.invalidateValid();
        end

        % performs a validation of all channel data against the specified channelDescriptors,
        % also fixing empty values appropriately. We will also adjust
        % channelDescriptors whose data classes do not match the data in
        % data (having convert data to match that class)
        function td = validateData(td, varargin)
            % wrapper for validateDataInternal that allows this to be done
            % internally
            td.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addParameter('progress', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            progress = p.Results.progress;

            [td.data, td.channelDescriptorsByName] = td.validateDataInternal(td.data, td.channelDescriptorsByName, 'progress', progress, p.Unmatched);

            td = td.fixCheckAnalogDataMatchesTimeVectors('progress', progress);

            td = td.fixAnalogChannelGroups();
        end

        function [data, channelDescriptorsByName] = validateDataInternal(td, data, channelDescriptorsByName, varargin) %#ok<INUSL>
            p = inputParser();
            p.addParameter('progress', true, @islogical);
            p.addParameter('addMissingFields', false, @islogical); % if true, don't complain about missing channels, just add the missing fields
            p.addParameter('suppressWarnings', false, @islogical); % don't warn about any minor issues
            p.parse(varargin{:});
            suppressWarnings = p.Results.suppressWarnings;

            progress = p.Results.progress;
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
                impl = chd.getImpl();

                % check to make sure all fields were provided as expected
                [ok(iChannel), missing{iChannel}] = impl.checkData(data);
                chDescs{iChannel} = chd.describe();
                required(iChannel) = chd.required;

                if ~ok(iChannel)
                    data = impl.addMissingFields(data);
                    if ~suppressWarnings
                        if ~chd.required
                            % fill in missing values for optional channels but
                            % issue a warning
                            tcprintf('inline', '{yellow}Warning: {none}Missing optional channel {light blue}%s {none}fields {purple}%s\n', ...
                                chd.describe(), TrialDataUtilities.String.strjoin(missing{iChannel}, ', '));
                        else
                            tcprintf('inline', '{yellow}Warning: {none}Missing required channel {light blue}%s {none}fields {purple}%s\n', ...
                                chd.describe(), TrialDataUtilities.String.strjoin(missing{iChannel}, ', '));
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
                        chDescs{i}, TrialDataUtilities.String.strjoin(missing{i}, ', '));
                end
                error('Required channel data fields not provided by getChannelData');
            end

            % here we both change classes of data in memory and update
            % channelDescriptors depending on what we find
            if progress
                prog = ProgressBar(nChannels, 'Repairing and converting channel data');
            end
            for iChannel = 1:nChannels
                if progress
                    prog.update(iChannel, 'Repairing and converting %s', names{iChannel});
                end
                chd = channelDescriptorsByName.(names{iChannel});
                impl = chd.getImpl();
                if isa(chd, 'AnalogChannelDescriptor') && chd.isColumnOfSharedMatrix, continue, end
                for iF = 1:chd.nFields
                    fld = chd.dataFields{iF};

                    valueCell = {data.(fld)}';
                    [chd, valueCell] = impl.checkConvertDataAndUpdateMemoryClassToMakeCompatible(iF, valueCell);

                    if ~iscell(valueCell)
                        % nTrials x 1 vectors were num2cell'd inside
                        % assignIntoStructArray, byt nTrials x nValues
                        % param channels presented an issue
                        assert(size(valueCell, 1) == numel(data));
                        valueCell = mat2cell(valueCell, ones(numel(data), 1), size(valueCell, 2));
                    end
                    data = TrialDataUtilities.Data.assignIntoStructArray(data, fld, valueCell);
                end
                channelDescriptorsByName.(names{iChannel}) = chd;
            end

            if progress
                prog.finish();
            end
        end

        function [td, globalOkay] = fixCheckAnalogDataMatchesTimeVectors(td, varargin)
            p = inputParser();
            p.addParameter('progress', true, @islogical);
            p.parse(varargin{:});
            progress = p.Results.progress;
            
            td.warnIfNoArgOut(nargout);

            td = td.reset();

            [groupNames, ~, channelsByGroup] = td.listAnalogChannelGroups('includeTransformChannels', false);
            analogChNotInGroup = td.listAnalogChannels('includeNamedSubChannels', false, 'includeTransformChannels', false);

            globalOkay = true;

            nTotal = numel(analogChNotInGroup) + numel(groupNames);
            dataFields = cellvec(nTotal);
            timeFields = cellvec(nTotal);
            isGroup = falsevec(nTotal);
            groupIdx = nanvec(nTotal);
            [memoryClass, primaryChannelName] = stringvec(nTotal);
            for iA = 1:numel(analogChNotInGroup)
                ch = analogChNotInGroup{iA}; %#ok<*PROP>
                cd = td.channelDescriptorsByName.(ch);
                dataFields{iA} = cd.primaryDataField;
                timeFields{iA} = cd.timeField;
                isGroup(iA) = false;
                memoryClass{iA} = cd.memoryClassByField{1};
                primaryChannelName{iA} = ch;
            end

            i0 = numel(analogChNotInGroup);
            for iA = (1:numel(groupNames))
                % chList = td.listAnalogChannelsInGroup(groupNames{iA});
                cd = td.channelDescriptorsByName.(groupNames{iA});
                dataFields{i0+iA} = groupNames{iA};
                timeFields{i0+iA} = cd.timeField;
                isGroup(i0+iA) = true;
                groupIdx(i0+iA) = iA;
                memoryClass{i0+iA} = cd.memoryClassByField{1};
                primaryChannelName{iA} = groupNames{iA};
            end

            if progress
                prog = ProgressBar(nTotal, 'Checking sample count vs. times for analog channels');
            end
            for iA = 1:nTotal
                if progress
                    prog.update(iA, 'Checking sample count vs. times for %s', dataFields{iA});
                end

                dataField = dataFields{iA};
                timeField = timeFields{iA};

                okay = truevec(td.nTrials);
                [transpose_time, transpose_data] = falsevec(td.nTrials);

                for iT = 1:td.nTrials
                    sz = size(td.data(iT).(timeField));
                    if sz(2) > 1 && sz(1) == 1
                        transpose_time(iT) = true;
                        nTime = size(td.data(iT).(timeField), 2);
                    else
                        nTime = size(td.data(iT).(timeField), 1);
                    end

                    sz = size(td.data(iT).(dataField));
                    if sz(2) > 1 && sz(1) == 1 && ~isGroup(iA)
                        transpose_data(iT) = true;
                        nData = size(td.data(iT).(dataField), 2);
                    else
                        nData = size(td.data(iT).(dataField), 1);
                    end

                    okay(iT) = nTime == nData;
                end

                % do the actual transposing down here for efficiency
                if any(transpose_time)
                    tr_time = cellfun(@transpose, {td.data(transpose_time).(timeField)}, 'UniformOutput', false);
                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, tr_time, transpose_time);
                end

                if any(transpose_data)
                    tr_data = cellfun(@transpose, {td.data(transpose_data).(dataField)}, 'UniformOutput', false);
                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, tr_data, transpose_data);
                end

                if any(~okay)
                    if progress
                        prog.pause_for_output();
                    end
                    warning('%d trials have differing number of data samples in %s as timestamps in %s. Fixing by clearing data and time fields.', ...
                        nnz(~okay), dataField, timeField);
                    globalOkay = false;

                    if isGroup(iA)
                        td = td.copyRenameSharedChannelFields(channelsByGroup{groupIdx(iA)}, [false true]);
                        timeField = td.channelDescriptorsByName.(channelsByGroup{groupIdx(iA)}{1}).timeField; % update post rename
                        emptyVal = zeros(0, numel(channelsByGroup{groupIdx(iA)}), memoryClass{iA});
                    else
                        td = td.copyRenameSharedChannelFields(dataFields{iA}, [false true]);
                        timeField = td.channelDescriptorsByName.(dataFields{iA}).timeField; % update post rename

                        % for enum channels
                        if strcmp(memoryClass{iA}, 'char')
                            emptyVal = '';
                        else
                            emptyVal = zeros(0, 1, memoryClass{iA});
                        end
                    end

                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, zeros(0, 1), ~okay);
                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, dataField, emptyVal, ~okay);
                end

                timeData = {td.data.(timeField)}';
                emptyMask = cellfun(@(t) numel(t) < 2, timeData);
                timeDelta = median(cellfun(@(x) median(diff(x), 'omitnan'), timeData(~emptyMask)), 'omitnan');

                % check sorted and no duplicates
                resort = truevec(td.nTrials);
                for iT = 1:td.nTrials
                    timeThis = timeData{iT};
                    timeThis = TrialDataUtilities.Data.removeSmallTimeErrors(timeThis, timeDelta, 0);
                    resort(iT) = ~issorted(timeThis) || numel(timeThis) > numel(unique(timeThis));
                end
                if any(resort)
                    if progress
                        prog.pause_for_output();
                    end
                    warning('%d trials have duplicate or non-monotonically increasing timestamps for %s. Fixing by deleting non-montonic samples.', nnz(resort), dataField);

                    % sort all affected channels so that the shared time field can be preserved
                    [channelsSharingTime, whichField] = td.getChannelsReferencingFields(timeField);
                    if any(whichField ~= 2)
                        error('Time field shared with another channel but not as time field');
                    end

                    [timeSorted, timeSortIdx] = cellfun(@getTimeSorted, timeData(resort), 'UniformOutput', false);

                    for iF = 1:numel(channelsSharingTime)
                        cdShared = td.channelDescriptorsByName.(channelsSharingTime{iF});

                        if ~(isa(cdShared, 'AnalogChannelGroupDescriptor') || isa(cdShared, 'AnalogChannelDescriptor')) || cdShared.isTransform
                            continue;
                        end

                        dataField = cdShared.dataFieldPrimary;
                        dataSorted = cellfun(@sortData, {td.data(resort).(dataField)}', timeSortIdx, 'UniformOutput', false);
                        td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, dataField, dataSorted, resort);
                    end
                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, timeSorted, resort);
%
%
%                     if isGroup(iA)
%                         td = td.copyRenameSharedChannelFields(primaryChannelName{iA}, [false true]);
%                         timeField = td.channelDescriptorsByName.(primaryChannelName{iA}).timeField; % update post rename
%                     else
%                         td = td.copyRenameSharedChannelFields(primaryChannelName{iA}, [false true]);
%                         timeField = td.channelDescriptorsByName.(primaryChannelName{iA}).timeField; % update post rename
%                     end

                    % use the time vectors with small errors removed
%                     [timeInsert, dataInsert] = cellfun(@resortTimeDedup, timeData(resort), {td.data(resort).(dataField)}', 'UniformOutput', false);
%                     td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, timeInsert, resort);
%                     td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, dataField, dataInsert, resort);
                end

            end
            if progress
                prog.finish();
            end

            function [time, sortIdx] = getTimeSorted(time)
                time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, 0);
                [time, sortIdx] = unique(time, 'last');
            end

            function data = sortData(data, sortIdx)
                data = data(sortIdx, :, :, :, :, :, :, :, :);
            end

%             function [time, data] = resortTimeDedup(time, data)
%                 time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, 0);
%                 [time, idx] = unique(time, 'last');
%                 data = data(idx, :, :, :, :);
%             end
        end

        function td = fixAnalogChannelGroups(td)
            % first, this function sets the sample size of each analog
            % channel group correctly in the descriptor
            %
            % this function finds analog channels that are part of an
            % analog channel group but that have their own descriptor
            % it then incorporates this channels info into the
            % AnalogChannelGroupDescriptor directly.
            % This is a change made on 20181023

            td.warnIfNoArgOut(nargout);

            % first identify any channels missing size information
            groups = td.listAnalogChannelGroups();
            for iG = 1:numel(groups)
                gcd = td.channelDescriptorsByName.(groups(iG));
                if isempty(gcd.sampleSize)
                    sz = ChannelImpl.getCellElementSize({td.data.(gcd.dataFieldPrimary)}');
                    gcd.sampleSize = sz(2:end);
                    td.channelDescriptorsByName.(groups(iG)) = gcd;
                end
            end

            % first fix or incorporate each analog channel descriptor
            candidates = td.listAnalogChannels('includeNamedSubChannels', false);
            cds = td.getChannelDescriptor(candidates);
            for iA = 1:numel(candidates)
                if cds(iA).isColumnOfSharedMatrix
                    if ~ismember(cds(iA).primaryDataField, groups)
                        % add missing group with sample size info
                        sz = ChannelImpl.getCellElementSize({td.data.(cds(iA).primaryDataField)}');
                        gcd = cds(iA).buildGroupChannelDescriptor('sampleSize', sz(2:end));

                        td = td.addChannel(gcd, {}, 'ignoreDataFields', true, 'ignoreExisting', true);
                        groups = union(groups, gcd.name);
                    else
                        % set this channel's name within existing group
                        parent = cds(iA).primaryDataField;
                        gcd = td.channelDescriptorsByName.(parent);
                        if isempty(gcd.sampleSize)
                            sz = ChannelImpl.getCellElementSize({td.data.(cds(iA).dataFieldPrimary)}');
                            gcd.sampleSize = sz(2:end);
                        end
                        gcd = gcd.setSubChannelInfo(cds(iA).primaryDataFieldColumnIndex, cds(iA).name, cds(iA).unitsPrimary);
                        td.channelDescriptorsByName.(parent) = gcd;
                    end
                    td.channelDescriptorsByName = rmfield(td.channelDescriptorsByName, cds(iA).name);
                end
            end
        end

        function tf = hasMetaKey(td, key)
            tf = isstruct(td.datasetMeta) && isfield(td.datasetMeta, key);
        end

        function v = getMetaKey(td, key, default)
            if isstruct(td.datasetMeta) && isfield(td.datasetMeta, key)
                v = td.datasetMeta.(key);
            elseif nargin > 2
                v = default;
            else
                v = [];
            end
        end

        function td = setMetaKey(td, key, value)
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
        function list = expandChannelListGroups(td, list, keepGroupNamesInList)
            list = string(list);
            if nargin < 3
                keepGroupNamesInList = true;
            end
            mask = td.hasAnalogChannelGroup(list);

            groups = unique(list(mask));
            subCh = arrayfun(@(grp) td.listAnalogChannelsInGroup(grp), groups, 'UniformOutput', false);
            if ~keepGroupNamesInList
                list = list(~mask);
            end
            if ~isempty(subCh)
                list = union(list, cat(1, subCh{:}));
            end
        end

        function list = listFieldsReferencedExclusivelyByChannels(td, channels, ignoreSubChannels)
            % this function will not resolve issues arising from channel groups
            % if an analog channel group is included but not one of its sub channels
            % then the data fields of that group will be included if ignoreSubChannels == true,
            % and not included if false
            if nargin < 3
                ignoreSubChannels = false;
            end

            cdsAll = td.channelDescriptorsByName;

            if ignoreSubChannels
                channels = td.expandChannelListGroups(channels);
            end

            cdsInside = rmfield(cdsAll, setdiff(fieldnames(cdsAll), channels));
            cdsOutside = rmfield(cdsAll, intersect(fieldnames(cdsAll), channels));

            fieldsInside = cellfun(@(fld) makecol(cdsInside.(fld).dataFields), fieldnames(cdsInside), 'UniformOutput', false);
            fieldsOutside = cellfun(@(fld) makecol(cdsOutside.(fld).dataFields), fieldnames(cdsOutside), 'UniformOutput', false);

            fieldsInside = cat(1, fieldsInside{:});
            fieldsOutside = cat(1, fieldsOutside{:});

            list = setdiff(fieldsInside, fieldsOutside);
        end

        function saveFast(td, location, varargin)
            % saveFast(td, location, ['partitions', partitionStruct])
            %
            % a partition defines a set of channels / channel groups that will be saved separately from the main trial data .mat files
            % and can be independently loaded in when loadFast is called, when the partition name is passed to loadFast
            %
            % if partitionWaveforms is set true, this will create a partition named waveforms that will save the spike waveforms into a
            % separate partition from the spike data

            p = inputParser();
            p.addParameter('partitions', td.saveFastPartitionInfo, @isstruct);
            p.addParameter('partitionWaveforms', td.saveFastPartitionWaveforms, @islogical);
            p.addParameter('convertCategoricals', true, @islogical);
            p.addParameter('partialLoadParams', ["saveTag", "trialId"], @isstringlike);
            p.addParameter('separateNonPartitioned', true, @islogical); % if false, td without partitioned channels will be saved together, and only paritioned data will be split into trials
            p.addParameter('trialsPerChunk', td.saveFastTrialsPerChunk, @isscalar); % split elements into files with this many elements per file
            p.addParameter('validate', true, @islogical);
            p.addParameter('progress', true, @islogical);
            p.parse(varargin{:});

            trialsPerChunk = p.Results.trialsPerChunk;
            progress = p.Results.progress;
            
            if p.Results.validate
                td = td.validateData();
            end
            
            data = td.data;

            % save a separate meta file with these param values
            partialLoadParams = string(p.Results.partialLoadParams);
            partialLoadParams = intersect(partialLoadParams, td.listParamChannels());
            if ~isempty(partialLoadParams)
                partialLoadData = struct();
                for iF = 1:numel(partialLoadParams)
                    partialLoadData.(partialLoadParams{iF}) = td.getParamRaw(partialLoadParams{iF});
                end
            else
                partialLoadData = [];
            end

            % categoricals are super slow to load, here we convert them to
            % char strings just to save time, since they'll be fixed during
            % validateData on load
            if p.Results.convertCategoricals
                chList = td.listParamChannels('isCategorical', true);

                for iC = 1:numel(chList)
                    fld = td.channelDescriptorsByName.(chList{iC}).dataFieldPrimary;
                    vals = {data.(fld)}';
                    vals = cellfun(@char, vals, 'UniformOutput', false);
                    data = TrialDataUtilities.Data.assignIntoStructArray(data, fld, vals);
                end
            end

            keepfields = @(s, flds) rmfield(s, setdiff(fieldnames(s), flds));

            partitionMeta = struct();
            partitionFields = struct();

            % to support partitions, we need to define the fields of .data that are used within a partition,
            % that aren't needed by other channels, and strip these fields and channels off
            partitions = p.Results.partitions;
            partitionNames = fieldnames(partitions);

            % spike waveforms get handled separately. other partitions live at the channel boundary, but spike waveforms are
            % part of a spike channel
            if p.Results.partitionWaveforms
                assert(~ismember('waveforms', partitionNames), 'Partition named waveforms reserved for partitionWaveforms');

                spikeCh = [td.listSpikeChannels('includeArraySubChannels', false); td.listExplicitSpikeArrays()];
                if ~isempty(spikeCh)
                    wavefields = cellvec(numel(spikeCh));
                    cdWithWaveforms = struct();
                    hasWaveforms = falsevec(numel(spikeCh));
                    for iU = 1:numel(spikeCh)
                        ch = spikeCh{iU};
                        hasWaveforms(iU) = td.channelDescriptorsByName.(ch).hasWaveforms;

                        if hasWaveforms(iU)
                            wavefields{iU} = td.channelDescriptorsByName.(ch).waveformsField;
                            % hold onto the original channel descriptor to store in the partition meta
                            cdWithWaveforms.(ch) = td.channelDescriptorsByName.(ch);
                            % and strip the waveforms from the descriptor saved with the core trial data
                            % so that it will be overwritten if the waveforms partition is loaded
                            td.channelDescriptorsByName.(ch) = td.channelDescriptorsByName.(ch).removeWaveformsField();
                        end
                    end

                    if any(hasWaveforms)
                        partitionMeta.waveforms.channelDescriptorsByName = cdWithWaveforms;
                        partitionFields.waveforms = wavefields(hasWaveforms);
                    end
                end
            end

            channelDescriptorsByName = td.channelDescriptorsByName;

            % figure out field subset and channel descriptors for each user-specified partition
            for p = 1:numel(partitionNames)
                pname = partitionNames{p};

                chThisPartition = partitions.(pname);
                if ischar(chThisPartition)
                    chThisPartition = {chThisPartition};
                end

                % check that fields aren't analog channel groups
                mask = falsevec(numel(chThisPartition));
                for iC = 1:numel(chThisPartition)
                    mask(iC) = isfield(td.channelDescriptorsByName, chThisPartition{iC});
                    if mask(iC)
                        cd = td.channelDescriptorsByName.(chThisPartition{iC});
                        if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix && ~ismember(cd.groupName, chThisPartition)
                            error('%s is an analog channel in group %s but group %s is not included in the partition', cd.name, cd.groupName, cd.groupName);
                        end
                    end
                end
                chThisPartition = chThisPartition(mask);

                % include the group sub channel names
                chThisPartition = td.expandChannelListGroups(chThisPartition, true);
                fieldsThisPartition = td.listFieldsReferencedExclusivelyByChannels(chThisPartition);

                if isempty(fieldsThisPartition)
                    % no need to save empty partitions
                    continue;
                end
                
                % store the channel descriptors in partition meta
                partitionMeta.(pname) = struct('channelDescriptorsByName', keepfields(channelDescriptorsByName, chThisPartition));
                partitionFields.(pname) = fieldsThisPartition;

                % strip these fields off channelDescriptorsByName, leave them in data though since saveArray will do this
                % we don't want to strip these off of td here, because then
                % fields that overlap between two partitions will be
                % stripped off by the second partition, rather than left in
                % the core td object, which is what we want.
                channelDescriptorsByName = rmfield(channelDescriptorsByName, intersect(fieldnames(channelDescriptorsByName), chThisPartition));
            end

            td.data = 'saved separately, load with loadFast';
            td.channelDescriptorsByName = channelDescriptorsByName;
            td.odc = []; %#ok<MCHV2>

            mkdirRecursive(location);
            TrialDataUtilities.Save.savefast(fullfile(location, 'td.mat'), 'td');

            % save elements of data
            msg = sprintf('Saving TrialData to %s', location);
            TrialDataUtilities.Data.SaveArrayIndividualized.saveArray(location, data, 'message', msg, 'progress', progress, ...
                'elementsPerChunk', trialsPerChunk, ...
                'partitionFieldLists', partitionFields, 'partitionMeta', partitionMeta, 'partialLoadData', partialLoadData);
        end
    end

    % Loading from disk
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

        function td = loadFast(location, varargin)
            p = inputParser();
            p.addParameter('maxTrials', Inf, @isscalar);
            p.addParameter('partitions', {}, @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('loadAllPartitions', [], @(x) isempty(x) || islogical(x));
            p.addParameter('ignoreMissingPartitions', false, @islogical);
            p.addParameter('objectName', "TrialData", @isstringlike);
            p.addParameter('validate', true, @islogical);
            p.addParameter('progress', true, @islogical);
            p.KeepUnmatched = true; % unmatched fields will match parameter values that were provided to load fast
            p.parse(varargin{:});
            
            % strip extension
            %             [path, name, ext] = fileparts(location);
            %             location = fullfile(path, name);

            % by default, load all
            partitions = string(p.Results.partitions);
            loadAllPartitions = p.Results.loadAllPartitions;
            if isempty(loadAllPartitions)
                loadAllPartitions = isempty(partitions);
            end
            progress = p.Results.progress;
            obj_name = string(p.Results.objectName);

            if exist(location, 'file') == 2
                ld = load(location);
                if isfield(ld, 'td')
                    td = ld.td;
                elseif numel(fieldnames(ld)) == 1
                    fld = fieldnames(ld);
                    td = ld.(fld{1});
                else
                    error('File %s has no td variable and more than one variable was found');
                end
                td = td.rebuildOnDemandCache();

            elseif exist(location, 'dir')
                loaded = load(fullfile(location, 'td.mat'));
                td = loaded.td;

                % load elements of data
                msg = sprintf('Loading %s from %s', obj_name, location);
                [data, partitionMeta, partialLoadMask] = TrialDataUtilities.Data.SaveArrayIndividualized.loadArray(location, ...
                    'message', msg, 'progress', progress, 'maxElements', p.Results.maxTrials, ...
                    'partitions', partitions, 'loadAllPartitions', loadAllPartitions, ...
                    'ignoreMissingPartitions', p.Results.ignoreMissingPartitions, 'partialLoadSpec', p.Unmatched);
                td.data = TensorUtils.inflateMaskedTensor(makecol(data), 1, partialLoadMask); % we'll select the partially loaded trials later, but for now needs to be numel(td.data) == nTrials (full)

                % add channel descriptors in partitionMeta, overwriting existing channels if overlapping (used for partitionWaveforms)
                partitions = fieldnames(partitionMeta);
                for iF = 1:numel(partitions)
                    td.channelDescriptorsByName = structMerge(td.channelDescriptorsByName, partitionMeta.(partitions{iF}).channelDescriptorsByName);
                end

                if ~all(partialLoadMask)
                    td = td.selectTrials(partialLoadMask);
                end

                td = td.rebuildOnDemandCache();

                if p.Results.validate
                    td = td.validateData('progress', progress);
                end
            else
                error('Directory %s not found. Did you save with saveFast?', location);
            end

            function s = structMerge(s, s2)
                flds = fieldnames(s2);
                for f = 1:numel(flds)
                    s.(flds{f}) = s2.(flds{f});
                end
            end
        end

        function td = loadFastMetaOnly(location, varargin)
            p = inputParser();
            p.addParameter('partitions', {}, @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('loadAllPartitions', [], @(x) isempty(x) || islogical(x));
            p.addParameter('ignoreMissingPartitions', false, @islogical);
            p.parse(varargin{:});
            
            % returns TrialData without the .data field, which will provide
            % access to metadata and channel info but not the data

            if exist(location, 'file') == 2
                td = TrialData.loadFast(location);
            elseif exist(location, 'dir')
                loaded = load(fullfile(location, 'td.mat'));
                td = loaded.td;
                td.nTrialsManual = TrialData.loadFastTrialCount(location);
                
                % load and merge in any partition meta
                partitionMeta = TrialDataUtilities.Data.SaveArrayIndividualized.loadPartitionMeta(location, ...
                    'partitions', p.Results.partitions, 'loadAllPartitions', p.Results.loadAllPartitions, ...
                    'ignoreMissingPartitions', p.Results.ignoreMissingPartitions);
                partitions = fieldnames(partitionMeta);
                for iF = 1:numel(partitions)
                    td.channelDescriptorsByName = structMerge(td.channelDescriptorsByName, partitionMeta.(partitions{iF}).channelDescriptorsByName);
                end
                
            else
                error('Directory %s not found. Did you save with saveFast?', location);
            end
            
            function s = structMerge(s, s2)
                flds = fieldnames(s2);
                for f = 1:numel(flds)
                    s.(flds{f}) = s2.(flds{f});
                end
            end
        end

        function nTrials = loadFastTrialCount(location)
            if exist(location, 'file') == 2
                td = TrialData.loadFast(location);
                nTrials = td.nTrials;
            elseif exist(location, 'dir')
                nTrials = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(location);
            else
                error('Location %s not found', location);
            end
        end

        function list = loadFastListPartitions(location)
            list = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(location);
        end

        function params = loadFastListPartialLoadParams(location)
            params = TrialDataUtilities.Data.SaveArrayIndividualized.listPartialLoadFields(location);
        end

        function tf = loadFastIsValidLocation(location)
            if exist(location, 'file') == 2
                info = whos('-file', location);
                names = {info.name};
                tf = numel(names) == 1 || ismember('td', names);
            else
                tf = exist(location, 'dir') && ...
                     exist(fullfile(location, 'td.mat'), 'file') && ...
                     TrialDataUtilities.Data.SaveArrayIndividualized.isValidLocation(location);
            end
        end

        function saveFastLinkPartitionFromOtherLocation(locationNameRef, locationNameDest, varargin)
            TrialDataUtilities.Data.SaveArrayIndividualized.linkPartitionFromOtherLocation(locationNameRef, locationNameDest, varargin{:});
        end
    end

    % simple build from scratch type methods
    methods(Static)
        % copy everything over from the TrialDataInterface
        function td = buildEmptyWithTrialDurations(durations, varargin)
            td = TrialData.buildEmptyWithTrialStartTrialEnd(0*durations, durations, varargin{:});
        end

        function td = buildEmptyWithTrialStartTrialEnd(trialStart, trialEnd, varargin)
            p = inputParser();
            p.addParameter('timeUnitName', 'ms', @ischar);
            p.parse(varargin{:});

            tdi = ManualTrialDataInterface('timeUnitName', p.Results.timeUnitName, ...
                'TrialStart', makecol(trialStart), 'TrialEnd', makecol(trialEnd));
            td = TrialData(tdi, 'suppressWarnings', true);
        end

        function td = buildForAnalogChannelGroupTensor(groupName, data, varargin)
            % data is nTrials x T x C, time is T x 1
            p = inputParser();
            p.addOptional('time', (0:size(data, 2)-1)', @(x) isempty(x) || isnumeric(x) && isvector(x));
            p.addParameter('timeUnitName', 'ms', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            time = p.Results.time;

            % args are name, tensor (nTrials x nTime x nChannels)
            assert(isvector(time) && isnumeric(data));
            nTrials = size(data, 1);
            nTime = size(data, 2);
            if isempty(time)
                time = 0:nTime-1;
            else
                assert(numel(time) == nTime);
            end

            % add extra time so that padding takes the full window when
            % needed
            pad = (time(2) - time(1)) / 2;

            start = repmat(min(time) - pad, nTrials, 1);
            stop = repmat(max(time) + pad, nTrials, 1);

            td = TrialData.buildEmptyWithTrialStartTrialEnd(start, stop, 'timeUnitName', p.Results.timeUnitName);
            td = td.addAnalogChannelGroup(groupName, data, time, p.Unmatched);
        end

    end

    % general utils
    methods(Static)
        function delta = computeDeltaFromTimes(time)
            if ~iscell(time)
                delta = nanmedian(diff(time));
            else
                % median of medians is faster and close enough
                emptyMask = cellfun(@(t) numel(t) < 2, time);
                delta = nanmedian(cellfun(@(x) nanmedian(diff(x)), time(~emptyMask)));
            end
        end

        % general utility to send plots to the correct axis
        function [axh, unmatched] = getRequestedPlotAxis(varargin)
            if isa(varargin{1}, 'TrialData') % used to be non-static method
                varargin = varargin(2:end);
            end

            p = inputParser();
            p.addParameter('figh', [], @(x) isempty(x) || ishandle(x));
            p.addParameter('axh', [], @(x) isempty(x) || ishandle(x));
            p.KeepUnmatched = false;
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

            aa = AutoAxis.recoverForAxis(axh);
            if ~isempty(aa)
                aa.uninstall();
                delete(aa);
            end
            unmatched = p.Unmatched;
        end

        function tdCombined = mergeTrialsFromMultipleTrialData(varargin)
            if isempty(varargin)
                error('Please provide at least 1 argument');
            end
            tdCombined = varargin{1};
            assert(isa(tdCombined, 'TrialData'), 'Individual arguments must be TrialData instances');
            if numel(varargin) > 1
                % check that channels  in common are all equivalent
                TrialData.assertCommonChannelsEquivalentMultipleTrialData(varargin{:});

                % pass as a cell array of trial data
                tdCombined = tdCombined.addTrialsFromTrialData(varargin(2:end));
            end
        end

        function assertCommonChannelsEquivalentMultipleTrialData(varargin)
            % check that for every cchannel in common between varargin{:}
            % that the channel descriptors are equal

            if numel(varargin) < 2
                return;
            end

            chListEach = cellfun(@(td) td.listChannels('includeNamedSubChannels', false), varargin, 'UniformOutput', false);
            chListAll = chListEach{1};
            for i = 2:numel(varargin)
                chListAll = union(chListAll, chListEach{i});
            end

            okay = truevec(numel(chListAll));
            for iC = 1:numel(chListAll)
                ch = chListAll{iC};
                % check that channel descriptors are okay for all
                % channels held in common between any
                hasChMask = cellfun(@(list) ismember(ch, list), chListEach);
                if nnz(hasChMask) > 1
                    cdCell = cellfun(@(td) td.channelDescriptorsByName.(ch), varargin(hasChMask), 'UniformOutput', false);
                    okay(iC) = isequaln(cdCell{:});
                end
            end

            if ~all(okay)
                error('TrialData instances differ on ChannelDescriptor for channel(s): %s', TrialDataUtilities.String.strjoin(chListAll(~okay), ', '));
            end
        end
    end

    % printing and display utilities
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
            % split param and event channels by display group
            [paramDisplayGroups, paramChannelsByGroup, paramChannelsNotInGroup] = td.listParamChannelDisplayGroups();
            [eventDisplayGroups, eventChannelsByGroup, eventChannelsNotInGroup] = td.listEventChannelDisplayGroups();

            % parse analog channels into grouped and non grouped
            analogCh = td.listAnalogChannels('includeDerivedChannelTypes', false, 'includeNamedSubChannels', false);
            continuousCh = td.listContinuousNeuralChannels('includeGroupSubChannels', false);
            [contGroups, ~, contGroupChannels] = td.listContinuousNeuralChannelGroups();
            [groups, ~, groupChannels] = td.listAnalogChannelGroups('includeDerivedChannelTypes', false);
            [imgChannels, ~, imgGroupChannels] = td.listImageChannels();


            tcprintf('inline', '{yellow}Param: {none}%s\n', TrialDataUtilities.String.strjoin(paramChannelsNotInGroup, ', '));
            printDisplayGroups(paramDisplayGroups, paramChannelsByGroup);
            tcprintf('inline', '{yellow}Event: {none}%s\n', TrialDataUtilities.String.strjoin(eventChannelsNotInGroup, ', '));
            printDisplayGroups(eventDisplayGroups, eventChannelsByGroup);
            tcprintf('inline', '{yellow}Analog: {none}%s\n', TrialDataUtilities.String.strjoin(analogCh, ', '));
            printAnalogGroups(groups, groupChannels);

            str = '{yellow}Spike: {none}';

            % display spike channel arrays
            arrays = td.listExplicitSpikeArrays();
            for iA = 1:numel(arrays)
                cd = td.channelDescriptorsByName.(arrays{iA});
                if cd.hasWaveforms
                    wavesStr = '+w';
                else
                    wavesStr = '';
                end
                str = [str, sprintf('{bright blue}%s{none} (%d on %d electrodes){bright blue}%s', arrays{iA}, cd.nChannels, cd.nElectrodes, wavesStr)]; %#ok<AGROW>
                if iA < numel(arrays)
                    str = [str, ', ']; %#ok<AGROW>
                end
            end
            if numel(arrays) > 0
                str = [str, ', '];
            else
                str = [str, ' '];
            end
            tcprintf('inline', str);
            str = '';

            % display spike channels indicating waveforms
            spikeCh = td.listSpikeChannels('includeArraySubChannels', false);
            hasWaves = td.hasSpikeWaveforms(spikeCh);
            for iS = 1:numel(spikeCh)
                if ~hasWaves(iS)
                    str = [str, spikeCh{iS}]; %#ok<AGROW>
                else
                    str = [str, spikeCh{iS}, '{bright blue}+w{none}']; %#ok<AGROW>
                end
                if iS < numel(spikeCh)
                    str = [str, ', ']; %#ok<AGROW>
                end
            end
            tcprintf('inline', str);

            tcprintf('inline', '\n{yellow}Continuous Neural: {none}%s\n', TrialDataUtilities.String.strjoin(continuousCh, ', '));
            printAnalogGroups(contGroups, contGroupChannels, false);

            tcprintf('inline', '{yellow}Image: {none}\n');
            printAnalogGroups(imgChannels, imgGroupChannels);

            function printDisplayGroups(groups, groupChannels)
                for iG = 1:numel(groups)
                    sz = numel(groupChannels{iG});
                    tcprintf('inline', '  {bright blue}%s{none} (%s): %s\n', groups{iG}, ...
                        TrialDataUtilities.String.strjoin(string(sz), ','), ...
                        TrialDataUtilities.String.strjoin(groupChannels{iG}, ', '));
                end
            end

            function printAnalogGroups(groups, groupChannels, showNamedSubChannels)
                if nargin < 3
                    showNamedSubChannels = true;
                end
                if ~isempty(groupChannels) && showNamedSubChannels
                    groupsNamedMask = ~cellfun(@isempty, groupChannels);
                else
                    groupsNamedMask = false(size(groups));
                end

                for iG = 1:numel(groups)
                    if groupsNamedMask(iG)
                        sz = td.getAnalogChannelGroupSize(groups{iG});
                        tcprintf('inline', '  {bright blue}%s{none} (%s): %s\n', groups{iG}, ...
                            TrialDataUtilities.String.strjoin(string(sz), ','), ...
                            TrialDataUtilities.String.strjoin(groupChannels{iG}, ', '));
                    end
                end

                if any(~groupsNamedMask)
                    gstr = '  ';
                    for iG = 1:numel(groups)
                        if ~groupsNamedMask(iG)
                            if any(~groupsNamedMask(iG+1:end))
                                commaStr = ', ';
                            else
                                commaStr = [' ' newline];
                            end
                            sz = td.getAnalogChannelGroupSize(groups{iG});
                            gstr = [gstr, sprintf('{bright blue}%s {none}(%s)%s', groups{iG}, ...
                                TrialDataUtilities.String.strjoin(string(sz), ','), commaStr)]; %#ok<AGROW>
                        end
                    end
                    tcprintf('inline', gstr);
                end
            end
        end

        function disp(td)
            td.printDescriptionShort();
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end

        function tbl = summarizeInvalidCauseAsTable(td)
            temporary = ~td.valid & ~td.permanentlyInvalid;
            cause = td.invalidCause;

            tbl = table(cause, temporary);
            [tbl, ~, idx] = unique(tbl);
            nTrials = histcounts(idx, 1:(size(tbl, 1)+1))';
            tbl = addvars(tbl, nTrials, 'Before', 'cause');
            tbl = sortrows(tbl, 'nTrials', 'descend');
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
            % support nTrialsManual only for loadMetaOnly
            if isstruct(td.data) || isempty(td.nTrialsManual)
                nTrials = numel(td.data);
            else
                nTrials = td.nTrialsManual;
            end
        end

        function nTrials = get.nTrialsPermanentlyInvalid(td)
            nTrials = nnz(~td.getManualValid);
        end

        function v = get.permanentlyInvalid(td)
            v = ~td.getManualValid;
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

        function n = get.timeUnitsPerMs(td)
            n = td.timeUnitsPerSecond / 1000;
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

        function td = markTrialsPermanentlyInvalid(td, mask, reason, varargin)
            p = inputParser();
            p.addParameter('validOnly', true, @islogical);  % if false, mark trial invalid even if its not currently valid
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            mask = makecol(TensorUtils.vectorIndicesToMask(mask, td.nTrials));
            if p.Results.validOnly
                mask = mask & td.valid;
            else
                % don't override reason for trials already permanently
                % invalid
                mask = mask & ~td.permanentlyInvalid;
            end

            td.manualValid(mask) = false;
            if nargin > 2
                if isempty(td.manualInvalidCause)
                    td.manualInvalidCause = cellvec(td.nTrials);
                end
                reason = cellstr(reason);
                if numel(reason) == 1
                    td.manualInvalidCause(mask) = reason;
                else
                    td.manualInvalidCause(mask) = reason(mask);
                end
            end
            td = td.invalidateValid();
        end

        function td = markValidTrialsPermanentlyInvalid(td, reason)
            td.warnIfNoArgOut(nargout);
            if nargin < 2
                reason = '';
            end
            td = td.markTrialsPermanentlyInvalid(td.valid, reason, 'validOnly', true);
        end

        function td = markInvalidTrialsPermanentlyInvalid(td, reason)
            td.warnIfNoArgOut(nargout);
            if nargin < 2
                reason = '';
            end
            td = td.markTrialsPermanentlyInvalid(~td.valid, reason, 'validOnly', false);
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

        function td = markTrialsTemporarilyInvalid(td, mask, reason, varargin)
            p = inputParser();
            p.addParameter('validOnly', true, @islogical);  % if false, mark trial invalid even if its not currently valid
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            mask = makecol(TensorUtils.vectorIndicesToMask(mask, td.nTrials));
            if p.Results.validOnly
                mask = mask & td.valid;
            end

            tempValid = td.getTemporaryValid();
            tempValid(mask) = false;
            if nargin > 2
                assert(ischar(reason));
                td.temporaryInvalidCause(mask) = {reason};
            end
            td.temporaryValid = makecol(tempValid);
            td = td.invalidateValid();
        end

        function td = markValidTrialsTemporarilyInvalid(td, reason)
            td.warnIfNoArgOut(nargout);
            if nargin < 2
                reason = '';
            end
            td = td.markTrialsTemporarilyInvalid(td.valid, reason, 'validOnly', true);
        end

        function td = withTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            mask = TensorUtils.vectorIndicesToMask(mask, td.nTrials);
            td = td.markTrialsTemporarilyInvalid(~mask);
        end

%         function td = setTrialsTemporarilyInvalid(td, mask)
%             td.warnIfNoArgOut(nargout);
%
%             td.temporaryValid = TensorUtils.vectorIndicesToMask(mask, td.nTrials);
%             td = td.invalideValid();
%         end

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

        function td = resetHard(td)
            td.warnIfNoArgOut(nargout);
            % don't touch .manualValid. this is not consistent with what reset means
            % for TrialDataConditionAlign

            td = td.clearTrialsTemporarilyInvalid();
            td = td.setManualValidTo(true(td.nTrials, 1));
        end
        
        function saveTags = listSaveTags(td)
            saveTags = td.getParamUnique('saveTag');
        end

        function td = selectTrialsFromSaveTag(td, saveTags)
            td.warnIfNoArgOut(nargout);

            mask = ismember(td.getParam('saveTag'), saveTags);
            td = td.selectTrials(mask);
        end

        function td = withTrialsFromSaveTag(td, saveTags)
            td.warnIfNoArgOut(nargout);

            mask = ismember(td.getParam('saveTag'), saveTags);
            td = td.withTrials(mask);
        end
    end

    % Internal code warnings
    methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect', ...
                    class(obj));
            end
        end
    end

    % Utility methods related to valid trials
    methods(Access=protected) 
        function [valid, cause] = getManualValid(td)
            if isempty(td.manualValid)
                valid = truevec(td.nTrials);
                cause = cellstrvec(td.nTrials);
            else
                valid = makecol(td.manualValid);

                % generate the list of causes
                if isempty(td.manualInvalidCause)
                    cause = cellstrvec(td.nTrials);
                else
                    cause = td.manualInvalidCause;
                end

                emptyMask = cellfun(@isempty, cause);
                cause(~valid & emptyMask) = {'marked invalid manually'};
                cause(valid) = {''};
            end
        end

        function [valid, cause] = getTemporaryValid(td)
            if isempty(td.temporaryValid)
                valid = truevec(td.nTrials);
                cause = cellstrvec(td.nTrials);
            else
                valid = makecol(td.temporaryValid);

                % generate the list of causes, prefixed with (temporary)
                if isempty(td.temporaryInvalidCause)
                    cause = cellstrvec(td.nTrials);
                else
                    cause = td.temporaryInvalidCause;
                end

                for iT = 1:td.nTrials
                    if ~valid(iT)
                        if isempty(cause{iT})
                            cause{iT} = 'marked invalid temporarily';
                        end
                    end
                end
                cause(valid) = {''};
            end
        end

        function vals = replaceInvalidMaskWithValue(td, vals, value)
            vals = td.replaceMaskedValuesWithValue(vals, value, ~td.valid);
        end

        function vals = replaceInvalidOrEmptyWithValue(td, vals, value)
            vals = makecol(vals);
            empty = cellfun(@isempty, vals);
            vals = td.replaceMaskedValuesWithValue(vals, value, ~td.valid | empty);
        end

        function vals = replaceMaskedValuesWithValue(td, vals, value, mask) %#ok<INUSL>
            % (valid, :) notation is to allow vals to be high dimensional
            assert(isequal(size(vals, 1), numel(mask)));
            if iscell(vals)
                [vals{mask, :}] = deal(value);
            else
                if isempty(value)
                    value = NaN;
                end
                sel = TensorUtils.maskByDimCellSelectAlongDimension(size(vals), 1, mask);
                vals(sel{:}) = value;
            end
        end

        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
        
        function td = dropChannelFields(td, fieldsRemove)
            % for the removed data channels' fields, figure out which ones
            % aren't referenced by any other channels
            otherChannelsReferencingFields = td.getChannelsReferencingFields(fieldsRemove);
            maskRemove = cellfun(@isempty, otherChannelsReferencingFields);
            fieldsRemove = fieldsRemove(maskRemove);
            td.data = rmfield(td.data, fieldsRemove);
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
                td.manualInvalidCause = cellvec(td.nTrials);
            end
            if isempty(td.temporaryValid)
                td.temporaryValid = truevec(td.nTrials);
                td.temporaryInvalidCause = cellvec(td.nTrials);
            end
        end

        function tf = hasChannel(td, names, varargin)
            p = inputParser();
            p.addParameter('includeNamedSubChannels', true, @islogical);
            p.addParameter('includeIndexedSubChannels', true, @islogical);
            p.addParameter('channelDescriptorClass', 'ChannelDescriptor', @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('groupChannelDescriptorClass', 'ChannelDescriptor', @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('includeDerivedChannelTypes', true, @islogical);

            p.parse(varargin{:});

            names = string(names);
            channelList = td.listChannels('includeNamedSubChannels', p.Results.includeNamedSubChannels, ...
                'channelDescriptorClass', p.Results.channelDescriptorClass, ...
                'groupChannelDescriptorClass', p.Results.groupChannelDescriptorClass, ...
                'includeDerivedChannelTypes', p.Results.includeDerivedChannelTypes);
            tf = ismember(names, channelList);

            if any(~tf) && p.Results.includeIndexedSubChannels
                groupNames = td.listChannelGroups('channelDescriptorClass', p.Results.groupChannelDescriptorClass, ...
                    'includeDerivedChannelTypes', p.Results.includeDerivedChannelTypes);
                for iN = 1:numel(names)
                    if contains(names(iN), '(')
                        [parent, cidx] = ChannelDescriptor.parseIndexedChannelName(names(iN));
                        if ~isnan(cidx) && ismember(parent, groupNames)
                            nCh = td.channelDescriptorsByName.(parent).nChannels;
                            tf(iN) = TrialDataUtilities.Data.indexInRange(cidx, nCh);
                        end
                    end
                end
            end
        end

        function assertHasChannel(td, names, varargin)
            p = inputParser();
            p.addParameter('includeNamedSubChannels', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            names = string(names);
            tf = td.hasChannel(names, p.Results, p.Unmatched);
            if p.Results.includeNamedSubChannels
                assert(all(tf), 'TrialData does not have channel(s) %s', strjoin(names(~tf), ','));
            else
                assert(all(tf), ...
                    'TrialData does not have channel(s) %s or else they are sub-channels, which are not supported for this operation', strjoin(names(~tf), ','));
            end
        end


        function [parent, chidx] = findContainerChannelForSubChannel(td, names)
            [subNames, parentNames, subChIdx] = td.listNamedSubChannels();
            [parent, chidx] = arrayfun(@findSingle, string(names));

            function [parent, chidx] = findSingle(name)
                parent = "";
                chidx = nan;
                if contains(name, '(')
                    % treat as indexed channel
                    [p, c] = ChannelDescriptor.parseIndexedChannelName(name);
                    if isfield(td.channelDescriptorsByName, p)
                        cd = td.channelDescriptorsByName.(p);
                        if cd.hasSubChannel(name)
                            parent = string(p);
                            chidx = c;
                        end
                    end
                else
                    [tf, idx] = ismember(name, subNames);
                    if tf
                        parent = parentNames(idx);
                        chidx = subChIdx(idx);
                    else
                        parent = "";
                        chidx = NaN;
                    end
                end
            end
        end

        function cd = getChannelDescriptor(td, names)
            if isempty(names)
                cd = [];
                return;
            end
            names = string(names);
            subNames = [];

            for iN = 1:numel(names)
                name = names(iN);
                if isfield(td.channelDescriptorsByName, name)
                    cd(iN) = td.channelDescriptorsByName.(name); %#ok<AGROW>
                elseif contains(name, '(')
                    [ch, chidx] = ChannelDescriptor.parseIndexedChannelName(name);
                    if isfield(td.channelDescriptorsByName, ch)
                        gcd = td.channelDescriptorsByName.(ch);
                        cd(iN) = gcd.buildSubChannelDescriptor(chidx); %#ok<AGROW>
                    else
                        error('Could not find channel %s', name);
                    end
                else
                    if isempty(subNames)
                        % compute these once on demand
                        [subNames, parentNames, chidx] = td.listNamedSubChannels();
                    end
                    [tf, idx] = ismember(name, subNames);
                    if ~tf
                        error('Could not find channel or sub-channel %s', name);
                    end
                    parentCd = td.channelDescriptorsByName.(parentNames(idx));
                    cd(iN) = parentCd.buildSubChannelDescriptor(chidx(idx)); %#ok<AGROW>
                end
            end

            cd = makecol(cd);
        end

        % This Should be disabled! It is mostly a hack for fixing data
        % issues post-hoc
%         function td = setChannelDescriptor(td, name, cd)
%             td.warnIfNoArgOut(nargout);
%             td.assertHasChannel(name);
%             assert(isa(cd, 'ChannelDescriptor'));
%             td.channelDescriptorsByName.(name) = cd;
%         end

        function type = getChannelType(td, name)
            type = td.getChannelDescriptor(name).getType();
        end

        function units = getChannelUnitsPrimary(td, names)
            % return a string describing the units of a given channel
            names = string(names);

            function units = getUnitsSingle(name)
%                 if td.hasSpikeChannel(name)
%                     units = "spikes/sec";
%                 else
                    units = string(td.getChannelDescriptor(name).unitsPrimary);
                    if strcmp(units, 'enum') % MatUdp did this at one point
                        units = "";
                    end
%                 end
            end

            units = arrayfun(@(n) getUnitsSingle(n), names);
            uniq_units = unique(units);
            if numel(uniq_units) == 1
                units = uniq_units(1);
            else
                units = "mixed";
            end
        end

        function td = setChannelUnitsPrimary(td, name, units)
            td.warnIfNoArgOut(nargout);
            td.channelDescriptorsByName.(name) = td.getChannelDescriptor(name).setPrimaryUnits(units);
        end

        function tf = isChannelVectorizable(td, name)
            cd = td.getChannelDescriptor(name);
            tf = cd.isVectorizableByField(1);
        end

        function tf = isChannelNumericScalar(td, name)
            cd = td.getChannelDescriptor(name);
            tf = cd.isNumericScalarByField(1);
        end

        function tf = isChannelCategorical(td, name)
            cd = td.channelDescriptorsByName.(name);
            tf = ~cd.collectAsCellByField(1) && strcmp(cd.accessClassByField{1}, 'categorical');
        end
        
        function td= setChannelMeta(td, name, meta)
            td.warnIfNoArgOut(nargout);
            assert(isstruct(meta) && isscalar(meta));
            
            td.channelDescriptorsByName.(name).meta = meta;
        end

        function td = setChannelMetaKey(td, name, key, value)
            td.warnIfNoArgOut(nargout);
            assert(ischar(key), 'Key must be string');

            if ~isstruct(td.channelDescriptorsByName.(name).meta)
                td.channelDescriptorsByName.(name).meta = struct(key, value);
            else
                td.channelDescriptorsByName.(name).meta.(key) = value;
            end
        end

        function meta = getChannelMeta(td, name)
            meta = td.channelDescriptorsByName.(name).meta;
            if isempty(meta)
                meta = struct();
            end
        end

        function value = getChannelMetaKey(td, name, key)
            if isstruct(td.channelDescriptorsByName.(name).meta) && isfield(td.channelDescriptorsByName.(name).meta, key)
                value = td.channelDescriptorsByName.(name).meta.(key);
            else
                value = [];
            end
        end

        function [names, channelDescriptors] = listChannels(td, varargin)
            % this lists all channels (including groups and optionally their sub
            % channels) matching a certain channel descriptor class (and
            % optionally their subclasses)
            p = inputParser();
            p.addParameter('includeNamedSubChannels', true, @islogical);
            p.addParameter('includeOnlyNamedSubChannels', false, @islogical);
            p.addParameter('channelDescriptorClass', 'ChannelDescriptor', @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('groupChannelDescriptorClass', string([]), @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.parse(varargin{:});

            channelDescriptors = td.getChannelDescriptorArray();
            if isempty(channelDescriptors)
                names = string([]);
                return;
            end

            cdClassList = string(p.Results.channelDescriptorClass);
            if p.Results.includeDerivedChannelTypes
                if numel(cdClassList) > 1
                    isaMulti = @(obj) any(arrayfun(@(cls) isa(obj, cls), cdClassList));
                    mask = arrayfun(isaMulti, channelDescriptors);
                else
                    mask = arrayfun(@(cd) isa(cd, p.Results.channelDescriptorClass), channelDescriptors);
                end
            else
                isaMutliNonSubClass = @(obj) any(arrayfun(@(cls) strcmp(class(obj), cls), cdClassList));
                mask = arrayfun(isaMutliNonSubClass, channelDescriptors);
            end
            channelDescriptors = channelDescriptors(mask);
            names = string({channelDescriptors.name})';

            if p.Results.includeNamedSubChannels || p.Results.includeOnlyNamedSubChannels
                if nargout > 1
                    groupChannels = td.listChannelGroups('channelDescriptorClass', p.Results.groupChannelDescriptorClass);
                    [namesSub, ~, ~, cdsSub] = td.listNamedSubChannels(groupChannels);
                    if p.Results.includeOnlyNamedSubChannels
                        names = namesSub;
                        channelDescriptors = cdsSub;
                    else
                        names = [names; namesSub];
                        channelDescriptors = [channelDescriptors; cdsSub];
                    end
                else
                    namesSub = td.listNamedSubChannels();
                    if p.Results.includeOnlyNamedSubChannels
                        names = namesSub;
                    else
                        names = [names; namesSub];
                    end
                end
            end
        end

        function [groups, channelDescriptors, channelsByGroup] = listChannelGroups(td, varargin)
            % this lists all channel groups
            p = inputParser();
            p.addParameter('channelDescriptorClass', string([]), @(x) ischar(x) || isstring(x));
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.parse(varargin{:});

            if isempty(p.Results.channelDescriptorClass)
                cdc = ["AnalogChannelGroupDescriptor", "SpikeArrayChannelDescriptor"];
            else
                cdc = p.Results.channelDescriptorClass;
            end
            [groups, channelDescriptors] = td.listChannels('channelDescriptorClass', cdc, ...
                'includeDerivedChannelTypes', p.Results.includeDerivedChannelTypes, ...
                'includeNamedSubChannels', false);
            channelsByGroup = arrayfun(@(g) td.channelDescriptorsByName.(g).listNamedSubChannels(), groups, 'UniformOutput', false);
        end

        function [subNames, parentNames, subChIdx, channelDescriptors] = listNamedSubChannels(td, groupNames, varargin)

            if nargin < 2
                groupNames = td.listChannelGroups();
            end

            if isempty(groupNames)
                [subNames, parentNames] = deal(strings(0, 1));
                subChIdx = nan(0, 1);
                channelDescriptors = [];
            else
                [subNames, subChIdx] = arrayfun(@(g) td.channelDescriptorsByName.(g).listNamedSubChannels(), groupNames, 'UniformOutput', false);
                nSubEach = cellfun(@numel, subNames);
                subNames = cat(1, subNames{:});
                subChIdx = cat(1, subChIdx{:});

                parentNames = arrayfun(@(g, n) repmat(g, n, 1), groupNames, nSubEach, 'UniformOutput', false);
                parentNames = cat(1, parentNames{:});

                if nargout > 3 % this is time consuming so only do if requested
                    channelDescriptors = td.getChannelDescriptor(subNames);
                end
            end
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

        function [names, channelDescriptors] = listSpecialChannels(td, varargin)
            [names, channelDescriptors] = td.listChannels(varargin{:}, 'includeNamedSubChannels', false);
            mask = arrayfun(@(cd) cd.special, channelDescriptors);
            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

        function [names, channelDescriptors] = listNonSpecialChannels(td, varargin)
            [names, channelDescriptors] = td.listChannels(varargin{:});
            mask = arrayfun(@(cd) ~cd.special, channelDescriptors);
            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

        function durations = getDurationsRaw(td)
            starts = td.getEventRawFirst('TrialStart');
            ends = td.getEventRawFirst('TrialEnd');
            durations = ends - starts;
        end

        function durations = getValidDurations(td)
            % get the time window for each trial
            starts = td.getEventFirst('TrialStart');
            ends = td.getEventLast('TrialEnd');
            durations = ends-starts;
            durations(durations < 0) = 0;
            durations = td.replaceInvalidMaskWithValue(durations, NaN);
        end

        function [durations, containsBlanked] = getValidDurationsForSpikeChannel(td, unitName, varargin)
            p = inputParser();
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            % similar to getValidDurations, except factors in the blanking
            % region
            durations = td.getValidDurations();

            %blankIntervals = td.getSpikeBlankingRegions(unitName, 'combine', p.Results.combine);
            blankIntervals = [];
            if ~isempty(blankIntervals) 
                % blankIntervals are guaranteed to be non-overlapping and lie
                % within the alignment window, so we can just add up the
                % individual intervals
                totalFn = @(mat) sum(mat(:, 2) - mat(:, 1), 1);
                nUnits = size(blankIntervals, 2);
                blankDurations = zeros(td.nTrials, nUnits);
            
                for iT = 1:td.nTrials
                    for iU = 1:nUnits
                        if ~isempty(blankIntervals{iT, iU})
                            blankDurations(iT, iU) = totalFn(blankIntervals{iT, iU});
                        end
                    end
                end
                containsBlanked = blankDurations > 0;
                durations = durations - blankDurations;
                if any(durations < 0)
                    error('Internal issue with blanking region durations');
                end
            else
                containsBlanked = false(size(durations));
            end
        end

        function cds = getChannelDescriptorArray(td, channels)
            % this sort operation is what puts channel names in order
            % during display, don't remove
            if nargin < 2
                channels = sort(fieldnames(td.channelDescriptorsByName));
            end
            if isempty(channels)
                cds = [];
            else
                for iF = 1:length(channels)
                    cds(iF) = td.channelDescriptorsByName.(channels{iF}); %#ok<AGROW>
                end
            end
            cds = makecol(cds);
        end

        function [namesByField, fieldIdxEach] = getChannelsReferencingFields(td, fields)
            % channels may reference common data fields in .data to prevent
            % data duplication. This function returns the list of channels
            % which reference a particular data field. fields is a cellstr
            % and nameCell is a cell of cellstr
            %
            % Note: sub channels are not included. This is important for
            % resolving conflicts

            if ischar(fields) || (isstring(fields) && isscalar(fields))
                wasCell = false;
                fields = {char(fields)};
            else
                wasCell = true;
            end
            channels = fieldnames(td.channelDescriptorsByName);
            fieldsByChannel = structfun(@(cd) cd.dataFields, ...
                td.channelDescriptorsByName, 'UniformOutput', false);

            [namesByField, fieldIdxEach] = cellvec(numel(fields));
            for iF = 1:numel(fields)
                % don't need the assert since this might be called before
                % the field is written
                if ~isfield(td.data, fields{iF}), continue; end
                %                 assert(isfield(td.data, fields{iF}), 'TrialData does not have data field %s', fields{iF});
                [mask, whichField] = structfun(@(chFields) ismember(string(fields{iF}), string(chFields)), fieldsByChannel);
                namesByField{iF} = channels(mask);
                fieldIdxEach{iF} = whichField(mask);
            end
            if ~wasCell
                namesByField = namesByField{1};
                fieldIdxEach = fieldIdxEach{1};
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
            removeNames = setdiff(td.listChannels('includeNamedSubChannels', false), names);
            td = td.dropChannels(removeNames);
        end

        function td = dropChannel(td, name)
            td.warnIfNoArgOut(nargout);
            td = td.dropChannels(name);
        end

        function td = dropAnalogChannelGroup(td, groupName) %#ok<INUSD>
            error('use dropChannels instead');
%             td.warnIfNoArgOut(nargout);
% 
%             groupName = string(groupName);
%             groupName = groupName(td.hasAnalogChannelGroup(groupName));
% 
%             for iG = 1:numel(groupName)
%                 % first process entire groups to be removed
%                 cd = td.channelDescriptorsByName.(groupName{iG});
%                 td.channelDescriptorsByName = rmfield(td.channelDescriptorsByName, groupName{iG});
% 
%                 % for the removed data channels' fields, figure out which ones
%                 % aren't referenced by any other channels
%                 fieldsRemove = cd.dataFields;
%                 otherChannelsReferencingFields = td.getChannelsReferencingFields(fieldsRemove);
%                 maskRemove = cellfun(@isempty, otherChannelsReferencingFields);
%                 fieldsRemove = fieldsRemove(maskRemove);
% 
%                 td.data = rmfield(td.data, fieldsRemove);
%             end
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
            names = intersect(names, td.listChannels('includeNamedSubChannels', true));
            if isempty(names)
                return;
            end

            names = string(names);

            % first hold onto the to-be-removed channel descriptors
            cds = cellvec(numel(names));
            inGroup = falsevec(numel(names));
            groupPartialSignalsRemove = struct();
            for i = 1:numel(names)
                cds{i} = td.getChannelDescriptor(names{i});
                if isa(cds{i}, 'AnalogChannelDescriptor') && cds{i}.isColumnOfSharedMatrix
                    groupName = cds{i}.primaryDataField;
                    inGroup(i) = true;
                    if ~ismember(groupName, names) % if not already included as a whole group
                        if isfield(groupPartialSignalsRemove, groupName)
                            groupPartialSignalsRemove.(groupName) = union(groupPartialSignalsRemove.(groupName), names(i));
                        else
                            groupPartialSignalsRemove.(groupName) = names(i);
                        end
                    end
                elseif isa(cds{i}, 'SpikeChannelDescriptor') && cds{i}.isColumnOfArray
                    error('Not yet implemented');
                end
%                 elseif isa(cds{i}, 'AnalogChannelGroupDescriptor')
%                     isGroup(i) = true;
            end

            % remove whole groups first - this doesn't do anything special anymore
%             if any(isGroup)
%                 td = td.dropAnalogChannelGroup(names(isGroup));
%             end

            % remove the channel descriptors
            cdFieldsRemove = intersect(names(~inGroup), fieldnames(td.channelDescriptorsByName));
            if ~isempty(cdFieldsRemove)
                td.channelDescriptorsByName = rmfield(td.channelDescriptorsByName, cdFieldsRemove);
            end

            % for the removed data channels' fields, figure out which ones
            % aren't referenced by any other channels
            fieldsRemoveByChannel = cellfun(@(cd) makecol(cd.dataFields), ...
                cds, 'UniformOutput', false);
            fieldsRemove = unique(cat(1, fieldsRemoveByChannel{:}));

            otherChannelsReferencingFields = td.getChannelsReferencingFields(fieldsRemove);
            maskRemove = cellfun(@isempty, otherChannelsReferencingFields);
            fieldsRemove = fieldsRemove(maskRemove);

            if ~isempty(fieldsRemove)
                td.data = rmfield(td.data, fieldsRemove);
            end

            % remove signals from groups by trimming analog data columns
            groups = fieldnames(groupPartialSignalsRemove);
            for iG = 1:numel(groups)
                group = groups{iG};
                removeSubNames = groupPartialSignalsRemove.(group);

                td = td.removeColumnsFromAnalogChannelGroup(group, removeSubNames);
            end

            td = td.postDataChange(fieldsRemove);
        end

        function td = dropChannelsMatchingWildcard(td, search)
            list = td.listChannelsMatchingWildcard(search);
            td = td.dropChannels(list);
        end

        function td = dropChannelsMatchingRegex(td, search)
            list = td.listChannelsMatchingRegex(search);
            td = td.dropChannels(list);
        end

        function td = trimAllChannelsToTrialStartEnd(td)
            td.warnIfNoArgOut(nargout);
            td = td.trimAllChannelsRaw(); % defaults to TrialStart / TrialEnd
        end
        
        function td = trimAllChannelsRaw(td, varargin)
            % Timepoints that lie outside of TrialStart (or startTimes and TrialStop (or endTimes) will
            % never be accessible via getTimes since they will be filtered
            % out by the AlignInfo. Optionally, re-zero all trials to
            % zeroTimes; ignored if left empty;
            p = inputParser();
            p.addOptional('startTimes', {}, @(x) isempty(x) || isvector(x));
            p.addOptional('stopTimes', {}, @(x) isempty(x) || isvector(x));
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            % in TDCA this is the equivalent of getAlignedTimesCell, but we
            % don't have that infrastructure in TD, so we just do it here
            % directly
            trialStart = cellfun(@(x) x(1), td.getEventRaw('TrialStart'));
            if isempty(p.Results.startTimes)
                startTimes = trialStart;
            else
                startTimes = p.Results.startTimes;
            end

            trialEnd = cellfun(@(x) x(1), td.getEventRaw('TrialEnd'));
            if isempty(p.Results.stopTimes)
                stopTimes = trialEnd;
            else
                stopTimes = p.Results.stopTimes;
            end

            startTimes(~td.valid) = NaN;
            stopTimes(~td.valid) = NaN;

            debug('Trimming event channels\n');
            ev = setdiff(td.listEventChannels(), {'TrialStart', 'TrialEnd'});
            td = td.trimEventRaw(ev, startTimes, stopTimes);

            debug('Trimming spike channels\n');
            sp = cat(1, td.listSpikeChannels('includeArraySubChannels', false), td.listExplicitSpikeArrays());
            td = td.trimSpikeChannelRaw(sp, startTimes, stopTimes);

            debug('Trimming analog channels\n');
            an = td.listAnalogChannelsNonGrouped();
            td = td.trimAnalogChannelRaw(an, startTimes, stopTimes);

            debug('Trimming analog channel groups\n');
            ag = td.listAnalogChannelGroups();
            td = td.trimAnalogChannelGroupRaw(ag, startTimes, stopTimes);

            % move TrialStart and TrialEnd events inward to match
            trialStart = max(trialStart, startTimes);
            td = td.setEvent('TrialStart', trialStart, 'isAligned', false);

            trialEnd = min(trialEnd, stopTimes);
            td = td.setEvent('TrialEnd', trialEnd, 'isAligned', false);
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
            p.addOptional('times', {}, @(x) isempty(x) || ischar(x) || iscell(x) || isvector(x)); % char or time cell
            p.addParameter('timeField', '', @ischar);
            p.addParameter('timeDelta', [], @(x) ismempty(x) || isscalar(x)); % specify this to ensure the channel is evenly sampled
            p.addParameter('units', '', @ischar);
            p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're identical
            p.addParameter('isAligned', true, @islogical);
            p.addParameter('clearForInvalid', false, @islogical);
            p.addParameter('scaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('scaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('dataInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is
            p.addParameter('displayGroup', '', @isstringlike);
            p.parse(varargin{:});
            times = p.Results.times;
            values = p.Results.values;
            units = p.Results.units;
            isAligned = p.Results.isAligned;

            td.warnIfNoArgOut(nargout);

            if isempty(times)
                assert(~isempty(p.Results.timeField))
                timeField = p.Results.timeField;
                times = [];
            elseif ischar(times)
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
                        isAligned = false;

                    elseif isfield(td.data, timeField)
                        % use directly specified time field in .data
                        times = {td.data.(timeField)};
                        isAligned = false;
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
%                 else
%                     % check that we're not overwriting another channel's
%                     % time field
%                     if isfield(td.data, timeField)
%                         otherChannels = setdiff(td.getChannelsReferencingFields(timeField), name);
%                         if ~isempty(otherChannels)
%                             error('Analog channel time field %s conflicts with existing channel(s) %s. If you meant to reference their timeField, specify ''timeField'' parameter and leave times argument empty', ...
%                                 name, TrialDataUtilities.String.strjoin(otherChannels, ','));
%                         end
%                     end
                end

                % also need to trim data coming in, in case it was
                % based off pre-trimmed data
                [times, values] = td.trimIncomingAnalogChannelData(times, values, 'isAligned', isAligned);
            end

            % build a channel descriptor for the data
            if p.Results.isContinuousNeural
                cd = ContinuousNeuralChannelDescriptor.buildVectorAnalogFromValues(name, timeField, units, td.timeUnitName, values, times);
            else
                cd = AnalogChannelDescriptor.buildVectorAnalogFromValues(name, timeField, units, td.timeUnitName, values, times);
            end

            cd.scaleFromLims = p.Results.scaleFromLims;
            cd.scaleToLims = p.Results.scaleToLims;
            cd.displayGroup = p.Results.displayGroup;
            
            td = td.addChannel(cd);

            if ~isempty(values)
                td = td.setAnalog(name, values, times, ...
                    'clearForInvalid', p.Results.clearForInvalid, 'isAligned', isAligned, ...
                    'keepScaling', true, 'dataInMemoryScale', p.Results.dataInMemoryScale);
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

        function td = addAnalogTransform(td, name, fromChannelNames, transformFn, varargin)
            % td = td.addAnalog(td, groupName, values, times)
            % values may be cell of matrices
            % corresponding to each trial.
            % times may be cell of vectors or single vector (for the matrix
            % values). Alternatively, times may be blank, and then the
            % parameter 'timeField' can specify which time field in which
            % to find the times for this channel
            p = inputParser();
            p.addParameter('units', '', @ischar);
%             p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're subclassed
%             p.addParameter('isImage', false, @islogical); % shortcut for making image channels since they're subclassed
            p.parse(varargin{:});

            chList = td.getChannelsReferencingFields(name);
            if ~isempty(chList)
                error('Channels %s already reference field %s', TrialDataUtilities.String.strjoin(chList));
            end
            if ~iscell(fromChannelNames)
                fromChannelNames = {fromChannelNames};
            end
            td.assertHasChannel(fromChannelNames);
            transformGroupDescriptors = td.getChannelDescriptor(fromChannelNames);
            cd = AnalogChannelDescriptor.buildTransformAnalogChannel(name, transformGroupDescriptors, transformFn, p.Results);
            td = td.addChannel(cd, {}, 'ignoreDataFields', true);
        end

        function td = setAnalog(td, name, values, varargin)
            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);

            cd = td.channelDescriptorsByName.(name);
            assert(~cd.isTransform, 'Can not set tranform channels');

            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) ||  isnumeric(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', true, @islogical); % if false, make the channel unscaled and convert the data to access class
            p.addParameter('dataInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is

            % these pass thru to setChannelData but are useful for error
            % checking here
            p.addParameter('updateValidOnly', true, @islogical);

            % for internal use mostly, used to indicate that the data provided is for the full trial,
            % not just the aligned window. this is for when times are not
            % provided, then when this flag is true, we know that the times
            % are not being updated (since we're not losing times outside
            % the alignment window). Used currently by
            % setAnalogWithinAlignWindow
            p.addParameter('timesMatchFullTrialRaw', false, @islogical);

            p.addParameter('separateChannelFromGroupIfNeeded', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            times = makecol(p.Results.times);

            % check the values and convert to nTrials cellvec
            if ismatrix(values) && (isnumeric(values) || islogical(values))
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

            if p.Results.updateValidOnly
                % don't worry about data counts not matching for invalid
                % trials if we're not updating them
                validMask = td.valid;
                times(~validMask) = {[]};
                values(~validMask) = {[]};
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
                if p.Results.separateChannelFromGroupIfNeeded
                    debug('Separating channel %s from column of shared matrix %s\n', name, cd.primaryDataField);
                    td = td.separateAnalogChannelFromGroup(name, 'separateTimeField', true); % no separate time field, no copy data
                    cd = td.channelDescriptorsByName.(name);
                else
                    error('Setting analog channel while ''keepScaling'' is false or when specifying new sample times requires this channel to be separated from its analog channel group. Use separateAnalogChannelFromGroup if you want to do this. Or use setAnalogChannelGroup to set all channels in the group at once.');
                end
            end

            if ~p.Results.keepScaling
                % data being passed in is now in original units
                % so change scaling factors
                cd = cd.withNoScaling();

                if p.Results.dataInMemoryScale
                    impl = cd.getImpl();
                    values = impl.convertMemoryDataCellToAccess(1, values);
                end
            else
                if ~p.Results.dataInMemoryScale
                    % take new data back into scaled values
                    impl = cd.getImpl();
                    values = impl.convertAccessDataCellToMemory(1, values);
                end
            end

            % update the channel descriptor accordingly
            td.channelDescriptorsByName.(name) = cd;

            % setChannelData will call repairData which will update
            % memoryDataClassByField{1} to reflect the type of values
            if updateTimes
                td = td.setChannelData(name, {values, times}, 'updateValidOnly', p.Results.updateValidOnly, p.Unmatched);
            else
                td = td.setChannelData(name, {values}, 'updateValidOnly', p.Results.updateValidOnly, p.Unmatched);
            end
        end

        function td = addOrUpdateAnalog(td, name, data, times, varargin)
            % set time samples of channel group if it exists where mask is true.
            % By default mask is non-empty cells in times
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', ~cellfun(@isempty, data), @isvector);
            %             p.addOptional('times', [], @(x) iscell(x) ||  ismatrix(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', false, @islogical);

            p.addParameter('dataInMemoryScale', false, @islogical);

            p.addParameter('timeField', '', @ischar);
            p.addParameter('units', '', @ischar);
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);

            p.parse(varargin{:});

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;

            assert(iscell(data) && numel(data) == td.nTrials, 'Data must be nTrials cell vector');

            td.warnIfNoArgOut(nargout);
            if td.hasAnalogChannel(name)
                td = td.setAnalog(name, data, times, 'updateMask', mask, ...
                    'isAligned', p.Results.isAligned, 'dataInMemoryScale', p.Results.dataInMemoryScale, ...
                    'separateChannelFromGroupIfNeeded', true);
            else
                % clear masked out cells
                [times{~mask}] = deal([]);
                [data{~mask}] = deal([]);
                td = td.addAnalog(name, data, times, ...
                    TrialDataUtilities.Data.keepfields(p.Results, {'timeField', 'units', 'isContinuousNeural', 'isAligned', ...
                    'scaleFromLims', 'scaleToLims', 'dataInMemoryScale'}));
            end
        end

        function td = addOrUpdateContinuousNeural(td, name, data, times, varargin)
            td.warnIfNoArgOut(nargout);

            td = td.addOrUpdateAnalog(name, data, times, varargin{:}, 'isContinuousNeural', true);
        end

        function td = scaleAnalogTimeField(td, names, multiplyBy, varargin)
            % pass the full list of
            td.warnIfNoArgOut(nargout);
            assert(isscalar(multiplyBy) && isnumeric(multiplyBy));

            td = td.copyRenameSharedChannelFields(names, 2);
            timeField = td.channelDescriptorsByName.(names{1}).timeField;

            prog = ProgressBar(td.nTrials, 'Scaling time field %s', timeField);
            for iT = 1:td.nTrials
                prog.update(iT);
                td.data(iT).(timeField) = multiplyBy * td.data(iT).(timeField);
            end
            prog.finish();

            td = td.postDataChange(timeField);
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

        function tf = hasAnalogChannel(td, name, allowSubChannel)
            if nargin < 3
                allowSubChannel = true; %
            end

            tf = td.hasChannel(name, 'includeNamedSubChannels', allowSubChannel, ...
                'channelDescriptorClass', "AnalogChannelDescriptor", 'groupChannelDescriptor', 'AnalogChannelGroupDescriptor');
        end

        function assertHasAnalogChannel(td, name)
            assert(td.hasAnalogChannel(name), 'No analog channel %s found', name);
        end

        function timeField = getAnalogChannelTimeField(td, name)
            td.assertHasChannel(name);
            timeField = td.channelDescriptorsByName.(name).timeField;
        end

        function [tf, timeField] = checkAnalogChannelsShareTimeField(td, names)
            names = string(names);
            timeFields = arrayfun(@(name) td.getAnalogChannelTimeField(name), names);
            tf = numel(unique(timeFields)) == 1;
            if tf
                timeField = timeFields(1);
            else
                timeField = "";
            end
        end

        function [names, channelDescriptors] = listAnalogChannels(td, varargin)
            % lists analog channels that have a named entry in the
            % channelDescriptors struct. this includes lone channels as
            % well as channels in an AnalogChannelGroup that have their own
            % AnalogChannelGroup in the table
            p = inputParser();
            p.addParameter('includeNamedSubChannels', true, @islogical);
            p.addParameter('channelDescriptorClass', 'AnalogChannelDescriptor', @ischar);
            p.addParameter('groupChannelDescriptorClass', 'AnalogChannelGroupDescriptor', @ischar);
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.addParameter('includeTransformChannels', true, @islogical);
            p.parse(varargin{:});

            [names, channelDescriptors] = td.listChannels(rmfield(p.Results, 'includeTransformChannels'));

            mask = true(numel(names), 1);
            if ~p.Results.includeTransformChannels
                mask = mask & arrayfun(@(cd) ~cd.isTransform, channelDescriptors);
            end

            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

%         function [subNames, parentNames, subChIdx, channelDescriptors]  = listAnalogNamedSubChannels(td)
%             groups = td.listAnalogChannelGroups();
%             [subNames, parentNames, subChIdx, channelDescriptors] = td.listNamedSubChannels(groups);
%         end

        function [names, cds] = listAnalogChannelsNonGrouped(td, varargin)
            [names, cds] = td.listAnalogChannels('includeNamedSubChannels', false, varargin{:});
        end

        function deltaMs = getAnalogTimeDeltaMs(td, name)
            delta = td.getAnalogTimeDelta(name);
            deltaMs = delta * td.timeUnitsPerSecond / 1000;
        end

        function delta = getAnalogTimeDelta(td, name)
            % compute the median delta betwen successive samples of an
            % analog channel(s), returns the minimum timeDelta across all channels

            name = string(name);

            delta = nanvec(numel(name));
            for i = 1:numel(name)
                cd = td.getChannelDescriptor(name{i});
                if (isa(cd, 'AnalogChannelDescriptor') || isa(cd, 'AnalogChannelGroupDescriptor')) && cd.isUniform
                    % trust the channel descriptor if it is a uniform channel
                    delta(i) = cd.timeDelta;
                else
                    time = td.getAnalogTime(name{i});
                    delta(i) = TrialData.computeDeltaFromTimes(time);
                end
            end

            % pick good sampling rate for all channels, since we'll be
            % resampling anyway, all we need to do is respect the nyquist
            % rates
            delta = nanmin(delta);
        end

        function fsHz = getAnalogSamplingRateHz(td, name)
            fsHz = td.timeUnitsPerSecond / td.getAnalogTimeDelta(name);
        end

        function timeField = getAnalogTimeField(td, name)
            cd = td.getChannelDescriptor(name);
            if isa(cd, 'AnalogChannelDescriptor') || isa(cd, 'AnalogChannelGroupDescriptor')
                timeField = string(cd.timeField);
            else
                error('%s is not an analog channel or analog channel group', name);
            end
        end

        function time = getAnalogTimeRaw(td, name)
            if td.hasAnalogChannel(name)
                timeField = td.getAnalogTimeField(name);
            elseif td.hasAnalogChannelGroup(name)
                timeField = td.getAnalogTimeField(name);
            elseif isfield(td.data, name)
                timeField = name;
            else
                error('%s is not an analog channel or time field', name);
            end

            time = {td.data.(timeField)}';
            time = cellfun(@makecol, time, 'UniformOutput', false);
        end

        function time = getAnalogTimeRawMinusOffsetFromZero(td, name)
            time = td.getAnalogTimeRaw(name);
            offsets = td.getTimeOffsetsFromZeroEachTrial();

            time = cellfun(@minus, time, num2cell(offsets), 'UniformOutput', false);
        end

        function [ch, idx, cdSubChannel] = parseIndexedAnalogChannelName(td, name)
            % if name is a normal analog channel, do nothing
            % if it is a name like analogGroup(5), make a sub channel
            % channelDescriptor for analogGroup referring to column 5

            tokens = regexp(name, '(?<ch>[\w_]+)\((?<idx>\d+)\)', 'names');
            if isempty(tokens)
                ch = name;
                idx = NaN;
                if td.hasChannel(name)
                    cdSubChannel = td.channelDescriptorsByName.(name);
                else
                    cdSubChannel = [];
                end
            else
                ch = tokens.ch;
                idx = str2double(tokens.idx);
                if td.hasAnalogChannelGroup(ch)
                    chName = sprintf('%s_%d', ch, idx);
                    cdSubChannel = td.channelDescriptorsByName.(ch).buildIndividualSubChannel(chName, idx);
                else
                    cdSubChannel = [];
                end
            end
        end

        function ch = generateIndexedAnalogChannelName(~, groupName, idx)
            ch = sprintf('%s(%d)', groupName, idx);
        end

        function [data, time] = getAnalogRaw(td, name, varargin)
            p = inputParser();
%             p.addParameter('sort', false, @islogical);
            p.addParameter('applyScaling', true, @islogical);

            % this functionality is not used by TDCA, which implements its own logic
            % but useful here for low level manipulations of analog
            % channels. Primarily for incorporating into channel groups
            p.addParameter('ensureUniformSampling', false, @islogical);
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @ischar); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @ischar);
            p.addParameter('sort', false, @islogical);
            p.KeepUnmatched = false;
            p.parse(varargin{:});

            cd = td.getChannelDescriptor(name);
            if isa(cd, 'AnalogChannelGroupDescriptor')
                [data, time] = td.getAnalogChannelGroup(name, varargin{:});
                return;
            end
            
            impl = cd.getImpl();
            assert(isa(cd, 'AnalogChannelDescriptor'), 'Channel %s is not analog', name);

            if cd.isTransform
                % now we do a transform of the data if the group is a virtual
                % transformation of another group. this function will handle
                % the slicing since it might be efficient not to compute the
                % elements that are not needed
                impl = cd.getImpl();
                [data, time] = impl.computeTransformDataRaw(td, 'applyScaling', p.Results.applyScaling, ...
                    'sort', p.Results.sort);
            else
                if cd.isColumnOfSharedMatrix
                    data = arrayfun(@(t) t.(cd.dataFields{1})(:, cd.primaryDataFieldColumnIndex), ...
                         td.data, 'UniformOutput', false, 'ErrorHandler', @(varargin) []);
                else
                    data = {td.data.(cd.dataFields{1})}';
                end
                time = {td.data.(cd.dataFields{2})}';

                for i = 1:td.nTrials
                    if numel(data{i}) == numel(time{i}) - 1
                        time{i} = makecol(time{i}(1:end-1));
                    else
                        time{i} = makecol(time{i});
                    end
                    data{i} = makecol(data{i});
                end

                if p.Results.applyScaling
                    % do scaling and convert to double
                    data = impl.convertDataCellOnAccess(1, data);
                end

%                 if p.Results.sort % this really shouldn't be necessary if data have been validated
%                     for iT = 1:td.nTrials
%                         [time{iT}, idx] = sort(time{iT}, 'ascend');
%                         data{iT} = data{iT}(idx);
%                     end
%                 end
            end

            % this functionality is not used by TDCA, which implements its own logic
            if p.Results.ensureUniformSampling || ~isempty(p.Results.timeDelta)
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = td.getAnalogTimeDelta(name);
                end
                [tMin, tMax] = td.getTimeStartStopEachTrialRaw();
                time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, p.Results.timeReference);

                [data, time] = TrialDataUtilities.Data.resampleDataCellInTime(data, time, 'timeDelta', timeDelta, ...
                    'timeReference', p.Results.timeReference, 'binAlignmentMode', p.Results.binAlignmentMode, ...
                    'resampleMethod', p.Results.resampleMethod, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'tMinExcludingPadding', tMin, 'tMaxExcludingPadding', tMax);
            end
        end

        function time = getAnalogTime(td, name)
            time = td.getAnalogTimeRaw(name);
        end

        % same as raw, except empty out invalid trials
        function [data, time] = getAnalog(td, name, varargin)
            [data, time] = td.getAnalogRaw(name, varargin{:});
            data = td.replaceInvalidOrEmptyWithValue(data, nan(0, 1));
            time = td.replaceInvalidOrEmptyWithValue(time, nan(0, 1));
        end

        function [dataVec, timeVec] = getAnalogSample(td, name, varargin)
            [dataCell, timeCell] = td.getAnalog(name);
            dataVec = cellfun(@(v) v(1), dataCell(td.valid), 'ErrorHandler', @(varargin) NaN);
            dataVec = TensorUtils.inflateMaskedTensor(dataVec, 1, td.valid, NaN);
            if nargout > 1
                timeVec = cellfun(@(v) v(1), timeCell(td.valid), 'ErrorHandler', @(varargin) NaN);
                timeVec = TensorUtils.inflateMaskedTensor(timeVec, 1, td.valid, NaN);
            end
        end

        function [dataUnif, timeUnif, delta] = getAnalogUniformlySampled(td, name, varargin)
%             [dataUnif, timeUnif, delta] = td.getAnalogRawUniformlySampled(name, varargin{:});
%             dataUnif = td.replaceInvalidMaskWithValue(dataUnif, []);
%             timeUnif = td.replaceInvalidMaskWithValue(timeUnif, []);


              [dataUnif, timeUnif] = td.getAnalog(name, varargin{:}, 'ensureUniformSampling', true);
              delta = TrialData.computeDeltaFromTimes(timeUnif);
        end

        function [dataCell, timeCell] = getAnalogMulti(td, name, varargin)
            % [data, time] = getAnalogMulti(td, chNames, varargin)
            % data and time are cell(nTrials, nChannels)
            p = inputParser();
            p.addParameter('raw', false, @islogical);
            p.addParameter('applyScaling', true, @islogical);
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
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalogRaw(name{c}, 'applyScaling', p.Results.applyScaling);
                else
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalog(name{c}, 'applyScaling', p.Results.applyScaling, p.Unmatched);
                end
            end
            prog.finish();
        end

        function [dataCell, timeCell] = getAnalogMultiCommonTime(td, names, varargin)
            % dataCell is cell(nTrials, 1) containing nTime x nChannels mat
            % timeCell is cell(nTrials, 1) containing nTime x 1 vectors
            % All data vectors will be interpolated to a
            % common time vector independently on each trial. Use
            % getAnalogMultiAsMatrix to register to a common time vector
            % across all trials.

            p = inputParser;
            p.addParameter('raw', false, @islogical);
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('interpolateMethod', 'linear', @ischar); % see interp1 for details
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp

            p.KeepUnmatched = true;
            p.parse(varargin{:});

            names = string(names);
            [inGroup, groupName] = td.checkAnalogChannelsInSameGroup(names);

            if inGroup
                % fetch the data as a group and rearrange the signals to
                % match the requested order

                % this is a much faster way of fetching the data whole, and
                % getAnalogChannelGroup will do the resampling
                dataCell = cellvec(td.nTrials);

                args = {groupName, 'timeDelta', p.Results.timeDelta, ...
                    'includeEdgeBins', true, 'timeReference', p.Results.timeReference, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'binAlignmentMode', p.Results.binAlignmentMode, 'resampleMethod', p.Results.resampleMethod, p.Unmatched};
                if p.Results.raw
                    [matCell, timeCell] = td.getAnalogChannelGroupRaw(args{:});
                else
                    [matCell, timeCell] = td.getAnalogChannelGroup(args{:});
                end

                % then go and grab the correct columns
                colIdx = td.getAnalogChannelColumnIdxInGroup(names);
                prog = ProgressBar(td.nTrials, 'Slicing columns from analog channel group');
                for iT = 1:td.nTrials
                    prog.update(iT)
                    dataCell{iT} = matCell{iT}(:, colIdx);
                end
                prog.finish();

            else
                % pick common sampling rate up front
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    % use common sampling rate and upsample
                    timeDelta = td.getAnalogTimeDelta(names); % will choose smallest sampling interval
                    if isempty(timeDelta)
                        timeDelta = td.alignInfoActive.minTimeDelta;
                    end
                end

                [dataCellRaw, timeCellRaw] = td.getAnalogMulti(names, 'raw', p.Results.raw, 'timeDelta', timeDelta, ...
                    'timeReference', p.Results.timeReference, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'binAlignmentMode', p.Results.binAlignmentMode, 'resampleMethod', p.Results.resampleMethod, p.Unmatched);

                % interpolate to common per-trial time vector in quick insertion mode
                % mat is nTrials x nTime
                [dataCell, timeCell] = deal(cell(td.nTrials, 1));

                for iT = 1:td.nTrials
                    if ~td.valid(iT), continue; end
                    % dataCell{iT} will be 1 x T x nChannels and needs to
                    % be squeezed along dim 1
                    [dataCell{iT}, timeCell{iT}]  = TrialDataUtilities.Data.embedTimeseriesInMatrix(...
                        dataCellRaw(iT, :), timeCellRaw(iT, :), ...
                        'assumeUniformSampling', true);  % since getAnalog already ensured uniform sampling
                    assert(size(dataCell{iT}, 1) == 1);
                    dataCell{iT} = permute(dataCell{iT}, [2 3 4 1]);
                end
            end
        end
    end

    methods % Analog Group and multi-channel analog methods
        function td = addAnalogChannelGroup(td, groupName, varargin)
            % td = td.addAnalog(td, groupName, chNames, values, times)
            % values may be cell of matrices
            % corresponding to each trial.
            % times may be cell of vectors or single vector (for the matrix
            % values). Alternatively, times may be blank, and then the
            % parameter 'timeField' can specify which time field in which
            % to find the times for this channel
            p = inputParser();
            p.addOptional('values', {}, @(x) iscell(x) || isnumeric(x) || islogical(x));
            p.addOptional('times', {}, @(x) isempty(x) || ischar(x) || iscell(x) || isvector(x)); % char or time cell
            p.addParameter('isAligned', true, @islogical);
            p.addParameter('clearForInvalid', false, @islogical);
            p.addParameter('scaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('scaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('dataInMemoryScale', false, @islogical); % treat data as in memory scaling and class (don't reverse the scaling)
            p.addParameter('raw', false, @islogical);

            % provide this manually, or provide the group of fields below
            p.addParameter('channelDescriptor', [], @(x) isempty(x) || isa(x, 'ChannelDescriptor'));

            % if channelDescriptor is provided, these fields are not used:
            p.addParameter('subChannelNames', [], @(x) isstring(x) || iscellstr(x));
            p.addParameter('subChannelUnits', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.addParameter('sampleSize', [], @(x) isempty(x) || isvector(x));
            p.addParameter('timeField', '', @ischar);
            p.addParameter('units', '', @isstringlike);
            p.addParameter('unitsByChannel', {}, @(x) isempty(x) || iscellstr(x) || isstring(x));
            p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're subclassed
            p.addParameter('continuousNeuralElectrodes', [], @(x) true);
            p.addParameter('isImage', false, @islogical); % shortcut for making image channels since they're subclassed
            p.addParameter('displayGroup', '', @ischar);
            p.parse(varargin{:});

            times = p.Results.times;
            values = p.Results.values;
            units = p.Results.units;
            isAligned = p.Results.isAligned;

            % create AnalogChannelDescriptors for each column?
            subChannelNames = string(p.Results.subChannelNames);
            subChannelUnits = string(p.Results.subChannelUnits);
            nameIndividualChannels = ~isempty(subChannelNames);

            td.warnIfNoArgOut(nargout);

            if ischar(times)
                % allow specifying char in times as well as timeField
                timeField = times;
                times = [];
            else
                timeField = p.Results.timeField;
            end

            if td.hasChannel(groupName)
                warning('Overwriting existing analog channel group %s', groupName);
                td = td.dropChannel(groupName);
            end

            existingSub = td.hasChannel(subChannelNames);
            if any(existingSub)
                warning('Overwriting existing analog channels %s', TrialDataUtilities.String.strjoin(subChannelNames(existingSub)));
                td = td.dropChannels(subChannelNames(existingSub));
            end

            % times can either be a field/channel name, or it can be raw
            % time values
            if isempty(times)
                assignTimesIntoData = false;
                % reference another time field or another analog channels
                % time field
                if ~isempty(timeField)
                    if td.hasChannel(timeField)
                        % treat timeField as analog channel name
                        % share that existing channel's time field
                        cd = td.channelDescriptorsByName.(timeField);
                        assert(isa(cd, 'AnalogChannelDescriptor') || isa(cd, 'AnalogChannelGroupDescriptor'), ...
                            'Channel %s is not an analog channel or channel group', timeField);
                        timeField = cd.timeField;
                        times = {td.data.(timeField)};
                        isAligned = false;

                    elseif isfield(td.data, timeField)
                        % use directly specified time field in .data
                        times = {td.data.(timeField)};
                        isAligned = false;
                    else
                        % timeField not found
                        error('%s is not a channel or data field name', timeField);
                    end

                else
                    error('Time vector or cell of vectors must be passed in when not referencing an existing timeField');
                end
            else
                % need to trim data coming in to TrialStart:TrialStop to
                % avoid conflicts with other channels
                [times, values] = td.trimIncomingAnalogChannelData(times, values, 'isAligned', isAligned);

                % times were provided
                % we'll set the field value for the new times field

                if isempty(timeField)
                    assignTimesIntoData = true;

                    % generate unique time field name
                    timeField = matlab.lang.makeUniqueStrings(sprintf('%s_time', groupName), fieldnames(td.data));
                else
                    if td.hasChannel(timeField)
                        % treat timeField as analog channel name
                        % share that existing channel's time field
                        cd = td.channelDescriptorsByName.(timeField);
                        assert(isa(cd, 'AnalogChannelDescriptor'), ...
                            'Channel %s is not an analog channel', timeField);
                        timeField = cd.timeField;
                    end

                    if isfield(td.data, timeField)
                        % use the other channel's time field to avoid duplicated timestamps
                        % but check that the time data is the same
                        assignTimesIntoData = false;

                        % need to trim the associated channels so that
                        % everything matches
                        td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeField);

                        % then check for equality
                        if isAligned
                            existingTimeData = td.getAnalogTimeRawMinusOffsetFromZero(timeField);
                        else
                            existingTimeData = td.getAnalogTimeRaw(timeField);
                        end

                        eqByTrial = cellfun(@(x, y) (isempty(x) && isempty(y)) || isequal(makecol(x), makecol(y)), times, existingTimeData);
                        assert(all(eqByTrial), 'Time passed in differ from existing times in .times on some trials. Note that incoming data is trimmed.');
                    else
                        % new, specified time field
                        assignTimesIntoData = true;
                    end
                end
            end

            if isfield(td.data, groupName)
                error('Field %s already exists in data structure', groupName);
            end

            sampleSize = p.Results.sampleSize;
            if isempty(sampleSize)
                sampleSize = ChannelImpl.getCellElementSize(values);
                sampleSize = sampleSize(2:end);
            end

            % assign an empty field
            if assignTimesIntoData
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, []);
            end
            td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, groupName, []);

            if isempty(p.Results.channelDescriptor)
                % loop over channels and create AnalogChannelDescriptors
                opts.scaleFromLims = p.Results.scaleFromLims;
                opts.scaleToLims = p.Results.scaleToLims;
                opts.dataClass = ChannelImpl.getCellElementClass(values);
                opts.timeClass = ChannelImpl.getCellElementClass(times);
                opts.sampleSize = sampleSize;

                % create a channel descriptor for the group as a whole
                if p.Results.isContinuousNeural
                    cdGroup = ContinuousNeuralChannelGroupDescriptor.buildAnalogGroup(groupName, timeField, ...
                        p.Results.continuousNeuralElectrodes, units, td.timeUnitName, opts);
                elseif p.Results.isImage
                    cdGroup = ImageChannelDescriptor.buildImage(groupName, timeField, units, td.timeUnitName, opts);
                else
                    cdGroup = AnalogChannelGroupDescriptor.buildAnalogGroup(groupName, timeField, units, td.timeUnitName, opts);
                end

                if nameIndividualChannels
                    cdGroup = cdGroup.setSubChannelInfo(1:numel(subChannelNames), subChannelNames, subChannelUnits);
                end
            else
                cdGroup = p.Results.channelDescriptor;
            end
            cdGroup.displayGroup = p.Results.displayGroup;
            td = td.addChannel(cdGroup, {}, 'ignoreDataFields', true);

            % set the data and time field in td.data
            if ~assignTimesIntoData
                times = []; % times weren't specified but copied from another channel, so there's no need to write into them in setAnalogChannelGroup
            end
            td = td.setAnalogChannelGroup(groupName, values, times, 'isAligned', isAligned, ...
                'keepScaling', true, 'dataInMemoryScale', p.Results.dataInMemoryScale, 'raw', p.Results.raw);
        end

        function td = addOrUpdateAnalogChannelGroup(td, groupName, data, times, varargin)
            % set time samples of channel group if it exists where mask is true.
            % By default mask is non-empty cells in times
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', ~cellfun(@isempty, data), @isvector);
            %             p.addOptional('times', [], @(x) iscell(x) ||  ismatrix(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', false, @islogical);

            p.addParameter('dataInMemoryScale', false, @islogical);

            p.addParameter('timeField', '', @ischar);
            p.addParameter('units', '', @ischar);
            p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're identical
            p.addParameter('isImage', false, @islogical);
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);

            p.KeepUnmatched = true;

            p.parse(varargin{:});
            td.warnIfNoArgOut(nargout);

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;

            td.warnIfNoArgOut(nargout);
            if td.hasAnalogChannelGroup(groupName)
                td = td.setAnalogChannelGroup(groupName, data, times, 'updateMask', mask, ...
                    'isAligned', p.Results.isAligned, 'dataInMemoryScale', p.Results.dataInMemoryScale);
            else
                % clear masked out cells
                [times{~mask}] = deal([]);
                [data{~mask}] = deal([]);
                td = td.addAnalogChannelGroup(groupName, data, times, ...
                    TrialDataUtilities.Data.keepfields(p.Results, {'timeField', 'units', 'isImage', 'isContinuousNeural', 'continuousNeuralElectrodes', 'isAligned', ...
                    'scaleFromLims', 'scaleToLims', 'dataInMemoryScale'}), p.Unmatched);
            end
        end

        function td = addAnalogChannelGroupTransform(td, groupName, fromGroupNames, transformFn, outputSize, varargin)
            % td = td.addAnalog(td, groupName, chNames, values, times)
            % values may be cell of matrices
            % corresponding to each trial.
            % times may be cell of vectors or single vector (for the matrix
            % values). Alternatively, times may be blank, and then the
            % parameter 'timeField' can specify which time field in which
            % to find the times for this channel
            p = inputParser();
            p.addParameter('units', '', @ischar);
%             p.addParameter('isContinuousNeural', false, @islogical); % shortcut for making LFP channels since they're subclassed
%             p.addParameter('isImage', false, @islogical); % shortcut for making image channels since they're subclassed
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            chList = td.getChannelsReferencingFields(groupName);
            if ~isempty(chList)
                error('Channels %s already reference field %s', TrialDataUtilities.String.strjoin(chList));
            end

            if ~iscell(fromGroupNames)
                fromGroupNames = {fromGroupNames};
            end
            td.assertHasChannel(fromGroupNames);
            transformGroupDescriptor = td.getChannelDescriptor(fromGroupNames);
            cd = AnalogChannelGroupDescriptor.buildTransformAnalogGroup(groupName, transformGroupDescriptor, transformFn, outputSize, p.Results);
            td = td.addChannel(cd, {}, 'ignoreDataFields', true);
        end

        function td = addAnalogChannelGroupAffineTransform(td, groupName, fromGroupNames, WinByOut, bOut, varargin)
            % creates a virtual transform analog channel group that
            % performs an affine transform of another group, as:
            % out = outputFn(in * W + b)

            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addParameter('outputFn', @(x) x, @(fn) isa(fn, 'function_handle')); % pointwise nonlinearity applied after
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [transformFn, szOut] = TrialDataUtilities.Transform.buildAnalogAffineTransformFn(WinByOut, bOut, 'outputFn', p.Results.outputFn);

            td = td.addAnalogChannelGroupTransform(groupName, fromGroupNames, transformFn, szOut, p.Unmatched);
        end

        function td = addAnalogChannelGroupSliceTransform(td, groupName, fromGroupNames, sliceArgs, varargin)
            % creates a virtual transform analog channel group that
            % performs an affine transform of another group, as:
            % out = outputFn(in * W + b)

            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addParameter('outputFn', @(x) x, @(fn) isa(fn, 'function_handle')); % pointwise nonlinearity applied after
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [transformFn, szOut] = TrialDataUtilities.Transform.buildAnalogSliceTransformFn(sliceArgs);

            td = td.addAnalogChannelGroupTransform(groupName, fromGroupNames, transformFn, szOut, p.Unmatched);
        end

        function td = setAnalogChannelGroupSubChannelNames(td, groupName, chNames, chIdx)
            td.warnIfNoArgOut(nargout);

            assert(td.hasAnalogChannelGroup(groupName));
            cdGroup = td.channelDescriptorsByName.(groupName);
            if ischar(chNames)
                chNames = {chNames};
            end
            if nargin < 4
                chIdx = 1:numel(chNames);
            end

            td.channelDescriptorsByName.(groupName) = cdGroup.setSubChannelInfo(chIdx, chNames);
            
            % drop the existing named channels
%             [oldNamedChannels, hasName] = td.listAnalogChannelsInGroupByColumn(groupName, chIdx);
%             if any(hasName)
%                 td = td.dropChannels(oldNamedChannels(hasName));
%             end
% 
%             for iCh = 1:numel(chNames)
%                 % build a channel descriptor for the data
%                 cd = cdGroup.buildIndividualSubChannel(chNames{iCh}, chIdx(iCh));
%                 td = td.addChannel(cd, {}, 'ignoreDataFields', true);
%             end
        end

        function [groupNames, channelDescriptors, channelsByGroup] = listAnalogChannelGroups(td, varargin)
            p = inputParser();
            p.addParameter('channelDescriptorClass', 'AnalogChannelGroupDescriptor', @ischar);
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.addParameter('includeTransformChannels', true, @islogical);
            p.parse(varargin{:});

            [groupNames, channelDescriptors, channelsByGroup] = td.listChannelGroups(rmfield(p.Results, 'includeTransformChannels'));

            mask = true(numel(groupNames), 1);
            if ~p.Results.includeTransformChannels
                mask = mask & ~arrayfun(@(cd) cd.isTransform, channelDescriptors);
            end

            groupNames = groupNames(mask);
            channelDescriptors = channelDescriptors(mask);
            channelsByGroup = channelsByGroup(mask);
        end

        function sz = getAnalogChannelGroupSize(td, groupName)
            % we ignore this call since it will fail anyway and it wastes a
            % lot of time doing lokup
            %td.assertHasAnalogChannelGroup(groupName);
            cd = td.channelDescriptorsByName.(groupName);
            if isa(cd, 'AnalogChannelGroup')
                sz = 1;
            else
                sz = cd.getSampleSize(td);
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

            cdGroup = td.getChannelDescriptor(groupName);

            chList = td.listAnalogChannelsInGroup(groupName);
            cdList = td.getChannelDescriptor(chList);

            if ~cdGroup.hasScaling
                % nothing to do
                return;
            end

            cdGroupNew = cdGroup.withNoScaling();

            % manually replace the .data field with the new values
            % we convert memory->access using the old cd, then
            % access->memory with the new cd.
            dataField = groupName;
            prog = ProgressBar(td.nTrials, 'Converting %s data to unscaled form', groupName);
            for t = 1:td.nTrials
                prog.update(t);
                unscaled = cdGroup.convertDataSingleOnAccess(1, td.data(t).(dataField));
                td.data(t).(dataField) = cdGroupNew.convertAccessDataSingleToMemory(1, unscaled);
            end
            prog.finish();

            % replace the group channel descriptor
            td.channelDescriptorsByName.(groupName) = cdGroupNew;

            % and replace all the sub channel descriptors
            if ~isempty(cdCell)
                newCdCell = arrayfun(@(cd) cd.withNoScaling(), cdList, 'UniformOutput', false);
                for i = 1:numel(newCdCell)
                    td.channelDescriptorsByName.(newCdCell{i}.name) = newCdCell{i};
                end
            end
        end

        function [tf, groupName] = checkAnalogChannelsInSameGroup(td, names)
            groupNames = arrayfun(@(name) td.getAnalogChannelGroupName(name), string(names));
            tf = numel(unique(groupNames)) == 1 && ~isempty(groupNames{1});
            if tf
                groupName = groupNames(1);
            else
                groupName = "";
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

            cd = td.getChannelDescriptor(groupName);
            [names, colIndex] = cd.listNamedSubChannels();

%             chList = td.listAnalogChannels();
%             cdCell = cellvec(numel(chList));
%             mask = falsevec(numel(chList));
%             for iC = 1:numel(chList)
%                 cd = td.channelDescriptorsByName.(chList{iC});
%                 if cd.isColumnOfSharedMatrix && strcmp(cd.primaryDataField, groupName)
%                     mask(iC) = true;
%                     cdCell{iC} = cd;
%                 end
%
%             end
%             names = chList(mask);
%             cdCell = cdCell(mask);
%
%             % order them by their column index
%             colIndex = cellfun(@(cd) cd.primaryDataFieldColumnIndex, cdCell);
%             [colIndex, idx] = sort(colIndex, 'ascend');
%             names = names(idx);
        end

        function [namesByColumn, hasName] = listAnalogChannelsInGroupByColumn(td, groupName, colIdx, varargin)
            p = inputParser();
            p.addParameter('channelNamesBySource', true, @islogical);
            p.parse(varargin{:});

            % namesByColumn will be either the name or '' if not named (or
            % groupName(###) if channelNamesBySource is true
            if nargin < 3
                sz = td.getAnalogChannelGroupSize(groupName);
                colIdx = reshape(1:prod(sz), TensorUtils.expandScalarSize(sz));
            end

            namesByColumn = cell(size(colIdx));
            hasName = false(size(colIdx));
            [names, namedColIndex] = td.listAnalogChannelsInGroup(groupName);
            for c = 1:numel(colIdx)
                [tf, which] = ismember(colIdx(c), namedColIndex);
                if tf
                    namesByColumn{c} = names{which};
                elseif p.Results.channelNamesBySource
                    namesByColumn{c} = td.generateIndexedAnalogChannelName(groupName, c);
                else
                    namesByColumn{c} = '';
                end
                hasName(c) = tf;
            end
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

            cdList = td.getChannelDescriptor(names);
            colIdx = arrayfun(@(cd) cd.primaryDataFieldColumnIndex, cdList);
        end

        function [tf, groupName] = isAnalogChannelInGroup(td, name, groupName)
            % [tf, groupName] = isAnalogChannelInGroup(td, name) % in any group?
            % tf = isAnalogChannelInGroup(td, groupName) % in specific group?

            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            if ~isa(cd, 'AnalogChannelDescriptor')
                tf = false;
                groupName = '';
                return;
            end

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
            name = string(name);
            td.assertHasChannel(name);
            cd = td.getChannelDescriptor(name);
            groupName = strings(numel(name), 1);
            for iN = 1:numel(name)
                if ~cd(iN).isColumnOfSharedMatrix
                    groupName(iN) = "";
                else
                    groupName(iN) = string(cd(iN).primaryDataField);
                end
            end
        end

        function td = separateAnalogChannelFromGroup(td, name, varargin)
            p = inputParser();
            p.addParameter('separateTimeField', false, @islogical);
            p.addParameter('filterColumnFromGroup', false, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            names = string(name);
            
            groupNames = td.getAnalogChannelGroupName(names);
            
            for iN = 1:numel(names)
                name = names(iN);
                td.assertHasChannel(name);
                groupName = groupNames(iN);

                if isempty(groupName)
                    % not in group
                    return;
                end

                colIdx = td.getAnalogChannelColumnIdxInGroup(name);

                group_cd = td.channelDescriptorsByName.(groupName);
                cd = group_cd.buildIndividualSubChannel(colIdx);

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

    %             % rename channels 
    %             group_cd = group_cd.setSubChannelInfo(colIdx, chNames);
    %             td.channelDescriptorsByName.(groupName) = group_cd;

                % unset this sub channel name
                td = td.setAnalogChannelGroupSubChannelNames(groupName, "", colIdx);

            end
            if p.Results.filterColumnFromGroup
                groupNames = unique(groupNames);
                for iN = 1:numel(groupNames)
                    td = td.removeUnnamedColumnsAnalogChannelGroup(groupNames(iN));
                end
            end
        end

        function td = removeUnnamedColumnsAnalogChannelGroup(td, groupName, varargin)
            td.warnIfNoArgOut(nargout);
            [~, colIdx] = td.listAnalogChannelsInGroup(groupName);
            td = td.filterColumnsAnalogChannelGroup(groupName, colIdx, varargin{:});
        end

        function td = removeColumnsFromAnalogChannelGroup(td, groupName, removeColIdx, varargin)
            td.warnIfNoArgOut(nargout);

            cd = td.getChannelDescriptor(groupName);
            if isnumeric(removeColIdx)
                colIdx = setdiff(1:cd.nChannels, removeColIdx);
            else
                removeSubNames = string(removeColIdx);
                colIdx = ~ismember(cd.subChannelNames, removeSubNames);
            end
            td = td.filterColumnsAnalogChannelGroup(groupName, colIdx, varargin{:});
        end

        function td = filterColumnsAnalogChannelGroup(td, groupName, colIdx, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addParameter('newSampleSize', [], @(x) isempty(x) || isvector(x));
            p.addParameter('progress', false, @islogical);
            p.parse(varargin{:});
            progress = p.Results.progress;
            
            cd = td.getChannelDescriptor(groupName);
            
            if isempty(colIdx)
                debug('Dropping analog channel group %s', groupName);
                td = td.dropChannel(groupName);
                return;
            end

            if progress, prog = ProgressBar(td.nTrials, 'Filtering columns of analog channel group %s', groupName); end
            for t = 1:td.nTrials
                if progress, prog.update(t); end
                if ~isempty(td.data(t).(groupName))
                    td.data(t).(groupName) = td.data(t).(groupName)(:, colIdx);
                end
            end
            if progress, prog.finish(); end

            % and update the channel descriptor
            cd = cd.filterSubChannels(colIdx, p.Results.newSampleSize);
            td.channelDescriptorsByName.(groupName) = cd;

            td = td.postDataChange(groupName);
        end

        function td = separateAnalogChannelsIntoSeparateGroup(td, names, newGroupName, varargin)
            p = inputParser();
            p.addParameter('filterColumnFromGroup', false, @islogical);
            p.parse(varargin{:});
            td.warnIfNoArgOut(nargout);
            assert(td.checkAnalogChannelsInSameGroup(names), 'All channels must be part of a group already');

            [chAll, ~, groupName] = td.listAnalogChannelsInGroupWith(names{1});
            if ~strcmp(groupName, newGroupName)
                % name is changing
                conflicts = td.getChannelsReferencingFields(newGroupName);
                assert(isempty(conflicts), 'Channel(s) %s already reference field %s. This field name cannot be used to store analog channel group.', TrialDataUtilities.String.strjoin(conflicts, ', '), newGroupName);

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
                if ~isempty(td.data(iT).(groupName))
                    td.data(iT).(newGroupName) = td.data(iT).(groupName)(:, colIdx);
                else
                    td.data(iT).(newGroupName) = zeros(0, numel(colIdx));
                end
            end
            prog.finish();


            % add a group descriptor for the created group
            cdGroup = td.getChannelDescriptor(groupName);
            newCdGroup = cdGroup.rename(newGroupName, false); % leave time field alone
            td = td.addChannel(newCdGroup, 'ignoreDataFields', true, 'ignoreExisting', true);

            % remove existing names from the old group
            td = td.setAnalogChannelGroupSubChannelNames(groupName, repmat("", numel(names), 1), colIdx);
            
            % update the channel descriptors to point to the new group
%             for c = 1:numel(names)
%                 cd = td.channelDescriptorsByName.(names{c});
%                 cd.primaryDataField = newGroupName;
%                 cd.primaryDataFieldColumnIndex = c;
%                 td.channelDescriptorsByName.(names{c}) = cd;
%             end

            % clear out the old channel
            if ~strcmp(groupName, newGroupName)
                if isempty(setdiff(chAll, names))
                    % no need for old field
                    td.data = rmfield(td.data, groupName);
                elseif p.Results.filterColumnFromGroup
                    % clean out the unused columns from the old field
                    td = td.removeUnnamedColumnsAnalogChannelGroup(groupName);
                end
            end

            td = td.postDataChange(newGroupName);
        end

        function td = renameAnalogChannelGroup(td, groupName, newGroupName)
            td.warnIfNoArgOut(nargout);
            chNames = td.listAnalogChannelsInGroup(groupName);
            td = td.separateAnalogChannelsIntoSeparateGroup(chNames, newGroupName);
        end

        function td = sortChannelsInAnalogChannelGroup(td, groupName, sort_order)
            td.warnIfNoArgOut(nargout);
            td = td.filterColumnsAnalogChannelGroup(groupName, sort_order);
        end

        function tf = hasAnalogChannelGroup(td, groupName)
            % this could be made faster
            groupName = string(groupName);
            tf = ismember(groupName, td.listAnalogChannelGroups());
        end

        function tf = hasAnalogChannelOrGroup(td, groupName)
            tf = td.hasAnalogChannel(groupName) | td.hasAnalogChannelGroup(groupName);
        end

        function assertHasAnalogChannelGroup(td, groupName)
            assert(all(td.hasAnalogChannelGroup(groupName)), 'No analog channel group %s found', TrialDataUtilities.String.strjoin(groupName));
        end

        function assertHasAnalogChannelOrGroup(td, groupName)
            assert(all(td.hasAnalogChannel(groupName) | td.hasAnalogChannelGroup(groupName)), ...
                'Analog channel or channel group %s found', TrialDataUtilities.String.strjoin(groupName));
        end

        function timeFields = getAnalogChannelGroupTimeField(td, groupNames)
            groupNames = string(groupNames);
            timeFields = stringvec(numel(groupNames));
            for iG = 1:numel(groupNames)
                if td.hasAnalogChannel(groupNames(iG))
                    timeFields(iG) = td.getAnalogTimeField(groupNames(iG));
                else
                    td.assertHasAnalogChannelGroup(groupNames(iG));
                    cd = td.channelDescriptorsByName.(groupNames(iG));
                    timeFields(iG) = string(cd.timeField);
                end
            end
        end

        function time = getAnalogChannelGroupTime(td, groupName)
            time = td.getAnalogChannelGroupTimeRaw(groupName);
            time = td.replaceInvalidMaskWithValue(time, []);
        end

        function time = getAnalogChannelGroupTimeRaw(td, groupName)
            td.assertHasAnalogChannelGroup(groupName);
            cd = td.channelDescriptorsByName.(groupName);
            timeField = cd.timeField;
            time = {td.data.(timeField)}';
        end

        function emptySlice = getAnalogChannelGroupEmptySlice(td, groupName, varargin)
            sz = td.getAnalogChannelGroupSize(groupName);
            emptySlice = nan([0 sz]);
        end

        function [data, time] = getAnalogChannelGroupRaw(td, groupName, varargin)
            p = inputParser();
            p.addParameter('sort', false, @islogical);

            p.addParameter('applyScaling', true, @islogical);

            % the order of operations here is:
            % scale the data if requested
            % do a transform of the data if this group is a virtual transform of another group
            % slice args get processed first to select into the samples
            % then the weightedCombination is used
            % then the averaging is performed
            p.addParameter('slice', {}, @(x) true); % this is used to index specifically into each sample
            p.addParameter('sliceSelectSubChannels', {}, @(x) true);
            p.addParameter('linearCombinationWeights', [], @(x) true); % alternatively, take a weighted combination over samples in the slice, size should be [size of analog channel, number of weighted combinations]

			% these apply to the weighted combination operation
			p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            p.addParameter('keepNaNIfAllNaNs', false, @islogical); % when replaceNaNWithZero is true, keep the result as NaN if every entry being combined is NaN
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical); % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis this enables nanmean like computations

            p.addParameter('averageOverSlice', false, @islogical); % average within each slice
            p.parse(varargin{:});

            td.assertHasAnalogChannelGroup(groupName);
            cd = td.channelDescriptorsByName.(groupName);

            sliceArgs = p.Results.slice;
            if isempty(sliceArgs) && ~isempty(p.Results.sliceSelectSubChannels)
                % lookup channels by name
                sliceSelectSubChannels = string(p.Results.sliceSelectSubChannels);
                [allNamedChannels, channelColIdx] = td.listAnalogChannelsInGroup(groupName);
                [tf, ind] = ismember(sliceSelectSubChannels, allNamedChannels);
                assert(all(tf), 'Some specified sliceSelectSubChannels not found in group %s', groupName);
                sliceArgs = {channelColIdx(ind)};
            elseif ~isempty(sliceArgs) && ~iscell(sliceArgs)
                sliceArgs = {sliceArgs};
            end
            sampleSize = cd.getSampleSize(td);
            emptySlice = nan([0 sampleSize]);

            if cd.isTransform
                % now we do a transform of the data if the group is a virtual
                % transformation of another group. this function will handle
                % the slicing since it might be efficient not to compute the
                % elements that are not needed
                impl = cd.getImpl();
                [data, time] = impl.computeTransformDataRaw(td, 'applyScaling', p.Results.applyScaling, ...
                    'slice', sliceArgs, 'sort', p.Results.sort);

            else
                % fetch the data directly
                cd = td.channelDescriptorsByName.(groupName);
                impl = cd.getImpl();
                timeField = cd.timeField;
                dataField = cd.dataFieldPrimary;

                data = {td.data.(dataField)}';
                time = {td.data.(timeField)}';

                mask_empty = true(numel(data), 1);
                for i = 1:numel(data)
                    time{i} = makecol(time{i});
                    if isempty(data{i}), time{i} = nan(0, 1); continue; end

                    % this situation can happen if a channel sharing the time
                    % vector with this channel is set and this trial was
                    % invalid when it was set.
                    if isempty(time{i}), data{i} = emptySlice; continue; end
                    assert(size(data{i}, 1) == numel(time{i}), 'Number of timepoints in data on trial %d does not match time', i);
                    assert(size(data{i}, 2) == prod(sampleSize), 'Number of channels on trial %d does not match channel count', i);
                    
                    mask_empty(i) = false;
                end

                % do scaling and convert to double
                if p.Results.applyScaling
                    data = impl.convertDataCellOnAccess(1, data);
                    emptySlice = impl.convertDataCellOnAccess(1, emptySlice);
                end

                if ~isempty(sliceArgs)
                    % take a slice through the data
                    for i = 1:numel(data)
                        if ~isempty(data{i})
                            data{i} = data{i}(:, sliceArgs{:});
                        end
                    end
                    emptySlice = emptySlice(:, sliceArgs{:});
                    data(mask_empty) = {emptySlice};
                end
                
                if p.Results.sort
                    for iT = 1:td.nTrials
                        [time{iT}, idx] = sort(time{iT}, 'ascend');
                        data{iT} = data{iT}(idx);
                    end
                end
            end

            if ~isempty(p.Results.linearCombinationWeights)
                % linearCombinationWeights  is size(slice) by nCombinations
                weights = p.Results.linearCombinationWeights;
                szWeights = TensorUtils.sizeNDims(weights, numel(sampleSize) + 1);
                assert(TensorUtils.compareSizeVectors(szWeights(1:numel(sampleSize)), sampleSize), 'Size of linearCombination weights must be [%s nCombinations]', vec2str(sampleSize));
                nCombinations = szWeights(numel(sampleSize) + 1);

                weightsNewByOld = reshape(weights, [prod(sampleSize) nCombinations])';
                emptySlice = nan(0, nCombinations);
                prog = ProgressBar(numel(data), 'Computing weighted combinations of analog channel group');
                for i = 1:numel(data)
                    prog.update(i);
                    if isempty(data{i})
                        data{i} = emptySlice;
                    else
                        nSamplesThis = size(data{i}, 1);

                        thisFlat = reshape(data{i}, [nSamplesThis, prod(sliceSize)]);
                        data{i} = TensorUtils.linearCombinationAlongDimension(thisFlat, 2, weightsNewByOld, ...
                            'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                            'keepNaNIfAllNaNs', p.Results.keepNaNIfAllNaNs, ...
                            'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN);
                    end
                end
                prog.finish();
            end

            if p.Results.averageOverSlice
                % average the whole slice down to a single timeseries
                for i = 1:numel(data)
%                     if isempty(data{i})
%                         data{i} = mean(empty;
%                     else
                        data{i} = nanmean(data{i}(:, :), 2);
%                     end
                end
            end
        end

        function [dataCell, timeCell] = getAnalogChannelGroupMulti(td, groupNames, varargin)
            % [dataCell, timeCell] = getAnalogChannelGroupMulti(td, chNames, varargin)
            % data and time are cell(nTrials, nChannels)
            % groupNames may include analog channels too
            p = inputParser();
            p.addParameter('raw', false, @islogical);

            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if ischar(groupNames)
                 groupNames = {groupNames};
            end

            % build nTrials x nChannels cell of data/time vectors
            C = numel(groupNames);
            [dataCell, timeCell] = deal(cell(td.nTrials, C));
            for c = 1:C
                gname = groupNames{c};
                if td.hasAnalogChannelGroup(gname)
                    if p.Results.raw
                        [dataCell(:, c), timeCell(:, c)] = td.getAnalogChannelGroupRaw(gname);
                    else
                        [dataCell(:, c), timeCell(:, c)] = td.getAnalogChannelGroup(gname, p.Unmatched);
                    end
                elseif td.hasAnalogChannel(gname)
                    if p.Results.raw
                        [dataCell(:, c), timeCell(:, c)] = td.getAnalogRaw(gname);
                    else
                        [dataCell(:, c), timeCell(:, c)] = td.getAnalog(gname, p.Unmatched);
                    end
                end
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
                cd = td.getChannelDeescriptor(groupName);
                tvec = makecol(td.data(trialInd).(cd.timeField));
                assert(size(dataMat, 1) == numel(tvec), 'Number of timepoints in data on trial %d does not match time', trialInd);
            end

            % do scaling and convert to double
            if p.Results.applyScaling
                dataMat = cd.convertDataSingleOnAccess(1, dataMat);
            end
        end

        function delta = getAnalogChannelGroupTimeDelta(td, groupName)
            % compute the median delta betwen successive samples of an
            % analog channel(s), returns the minimum timeDelta across all channels
            % given in units of .timeUnits

            if ischar(groupName)
                groupName = {groupName};
            end

            delta = nanvec(numel(groupName));
            for i = 1:numel(groupName)
                time = td.getAnalogChannelGroupTime(groupName{i});
                % median of medians is faster and close enough
                delta(i) = nanmedian(cellfun(@(x) nanmedian(diff(x)), time));
            end

            delta = nanmin(delta);
        end

        function fsHz = getAnalogChannelGroupSamplingRateHz(td, name)
            fsHz = td.timeUnitsPerSecond / td.getAnalogChannelGroupTimeDelta(name);
        end

        function deltaMs = getAnalogChannelGroupTimeDeltaMs(td, name)
            delta = td.getAnalogChannelGroupTimeDelta(name);
            deltaMs = delta * td.timeUnitsPerSecond / 1000;
        end

%         function [dataUnif, timeUnif, delta] = getAnalogChannelGroupRawUniformlySampled(td, groupName, varargin)
%             p = inputParser();
%             p.addParameter('method', 'linear', @ischar);
%             p.addParameter('delta', [], @(x) isscalar(x) || isempty(x)); % in ms
%             p.parse(varargin{:});
%
%             delta = p.Results.delta;
%             if isempty(delta)
%                 delta = td.getAnalogChannelGroupTimeDelta(groupName);
%             end
%             [data, time] = td.getAnalogChannelGroupRaw(groupName);  %#ok<*PROP>
%
%             [dataUnif, timeUnif] = cellvec(td.nTrials);
%
%             % get first and last valid sample for all trials
%             % very important to use this in case data has NaN samples
%             [tmins, tmaxs] = TrialDataUtilities.Data.getValidTimeExtents(time, data);
%             for iT = 1:td.nTrials
%                 if isempty(time{iT}) || isempty(data{iT})
%                     continue;
%                 end
%                 timeUnif{iT} = (tmins(iT):delta:tmaxs(iT))';
%                 if nnz(~isnan(data{iT})) > 5
%                     dataUnif{iT} = interp1(time{iT}, data{iT}, timeUnif{iT}, p.Results.method);
%                 else
%                     dataUnif{iT} = nan(size(timeUnif{iT}));
%                 end
%             end
%         end

        function [dataUnif, timeUnif, delta] = getAnalogChannelGroupUniformlySampled(td, groupName, varargin)
            % WARNING: This will be aligned if td is a trialdataconditionalign instance
            
%             [dataUnif, timeUnif, delta] = td.getAnalogChannelGroupRawUniformlySampled(groupName, varargin{:});
%             dataUnif = td.replaceInvalidMaskWithValue(dataUnif, []);
%             timeUnif = td.replaceInvalidMaskWithValue(timeUnif, []);

            [dataUnif, timeUnif] = td.getAnalogChannelGroup(groupName, varargin{:}, 'ensureUniformSampling', true);
            delta = TrialData.computeDeltaFromTimes(timeUnif);
        end

        function [data, time] = getAnalogChannelGroup(td, groupName, varargin)
            [data, time] = td.getAnalogChannelGroupRaw(groupName, varargin{:});

            % now handled within getAnalogChannelGroupRaw so that slice works properly
%             emptySlice = [];
%             for i = 1:numel(data)
%                 sz = size(data{i});
%                 if sz(2) > 1
%                     emptySlice = nan([0 sz(2:end)]);
%                 end
%             end
% 
%             data = td.replaceInvalidOrEmptyWithValue(data, emptySlice);
%             time = td.replaceInvalidOrEmptyWithValue(time, nan(0, 1));
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
            cd = td.channelDescriptorsByName.(groupName);
            assert(~cd.isTransform, 'Cannot modify transform groups');

            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) ||  ismatrix(x) && ~ischar(x));
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('keepScaling', true, @islogical); % if false, drop the scaling of the channel in memory and convert everything to access class
            p.addParameter('dataInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is
            p.addParameter('updateMask', true(td.nTrials, 1), @isvector);
            p.addParameter('units', '', @ischar);
            p.addParameter('raw', false, @islogical);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            times = makecol(p.Results.times);

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.updateMask), td.nTrials);
            if ~p.Results.raw % don't update invalid trials
                mask = mask & td.valid;
            end
            groupName = char(groupName);

            subChList = td.listAnalogChannelsInGroup(groupName);
            nCh = numel(subChList);

            % check the values and convert to nTrials cellvec
            if isnumeric(values) || islogical(values)
                % values must be nTrials x nTimes x nChannels
                assert(size(values, 1) == td.nTrials, 'Values as matrix must be nTrials along dimension 1');

                % convert to cellvec
                tensor = values;
                values = cellvec(td.nTrials);
                for r = 1:td.nTrials
                    values{r} = TensorUtils.squeezeDims(tensor(r, :, :), 1);
                end
                clear tensor;

            elseif iscell(values)
                assert(numel(values) == td.nTrials, 'Values as cell must have numel == nTrials');
                nCol = cellfun(@(x) size(x, 2), values);
                isEmpty = cellfun(@isempty, values);
                if nCh > 0 % if named channels
                    assert(all(nCol(~isEmpty) == nCh), 'All elements of values must have same number of columns as named channels (%d)', nCh);
                end
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
                if p.Results.dataInMemoryScale
                    % we need to un-scale this memory-scale data back to
                    % access scale
                    values = td.channelDescriptorsByName.(groupName).convertMemoryDataCellToAccess(1, values);
                end
                % data being passed in is now in original units
                % so change scaling factors of the channel and the data
                td = td.convertAnalogChannelGroupToNoScaling(groupName);

                % and convert back to memory anyway in case the data class
                % has changed, although this shouldn't do any scaling
                % since the channel descriptor has changed
                values = td.channelDescriptorsByName.(groupName).convertAccessDataCellToMemory(1, values);
            else
                if ~p.Results.dataInMemoryScale
                    % take new data back into scaled values to match the
                    % existing
                    cd = td.channelDescriptorsByName.(groupName);
                    impl = cd.getImpl();
                    values = impl.convertAccessDataCellToMemory(1, values);
                end
            end

            if updateTimes
                % avoid overwriting data shared by other channels,
                % including time field
                td = td.copyRenameSharedChannelFields(groupName, 2);
                timeField = td.channelDescriptorsByName.(groupName).timeField;
            end

%             prog = ProgressBar(td.nTrials, 'Writing analog channel group data on valid trials');
%             for t = 1:td.nTrials
%                 prog.update(t);
%                 if ~td.valid(t), continue; end
%
%                 td.data(t).(groupName) = values{t};
%                 if updateTimes
%                     td.data(t).(timeField) = times{t};
%                 end
%             end

            td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, groupName, values, mask);
            if updateTimes
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, times, mask);
            end

            if ~isempty(p.Results.units)
                td = setChannelUnitsPrimary(td, groupName, p.Results.units);
            end

            if updateTimes
                td = td.postDataChange({groupName, timeField});
            else
                td = td.postDataChange({groupName});
            end

%             prog.finish();
        end

        function td = trimAnalogChannelGroupToTrialStartEnd(td, names)
            td.warnIfNoArgOut(nargout);
            td = td.trimAnalogChannelGroupRaw(names); % defaults to trial start / end
        end

        function td = trimAnalogChannelGroupRaw(td, groupNames, varargin)
            td.warnIfNoArgOut(nargout);

            groupNames = string(groupNames);
            timeFields = td.getAnalogChannelGroupTimeField(groupNames);

            if ~isempty(timeFields)
                timeFields = unique(timeFields);
                prog = ProgressBar(numel(timeFields), 'Trimming analog channels');
                for i = 1:numel(timeFields)
                    prog.update(i, 'Trimming analog channels with time field %s', timeFields{i});
                    td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeFields{i}, varargin{:});
                end
                prog.finish();
            end
        end

        function td = trimAnalogChannelToTrialStartEnd(td, names)
            td.warnIfNoArgOut(nargout);
            td = td.trimAnalogChannelRaw(names); % defaults to trial start / end
        end

        function td = trimAnalogChannelRaw(td, names, varargin)
            td.warnIfNoArgOut(nargout);

            names = string(names);
            timeFields = unique(arrayfun(@(name) td.getAnalogTimeField(name), names));

            if ~isempty(timeFields)
                timeFields = unique(timeFields);
                prog = ProgressBar(numel(timeFields), 'Trimming analog channels');
                for i = 1:numel(timeFields)
                    prog.update(i, 'Trimming analog channels sharing time field %s', timeFields{i});
                    td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeFields{i}, varargin{:});
                end
                prog.finish();
            end
        end

        function td = trimAnalogChannelTimeFieldAndReferencingChannelsRaw(td, timeField, varargin)
            % delete time points outside of a certain start stop interval
            % defaults to TrialStart:TrialStop
            p = inputParser();
            p.addOptional('startTimes', {}, @isvector);
            p.addOptional('stopTimes', {}, @isvector);
            p.addParameter('clearInvalidTrials', true, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            times = {td.data.(timeField)}';
            notEmpty = ~cellfun(@isempty, times);

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

            if p.Results.clearInvalidTrials
                startTimes(~td.valid) = NaN;
                stopTimes(~td.valid) = NaN;
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
            cdList = td.getChannelDescriptor(chList);

            % separate groups from non groups
            maskIndiv = arrayfun(@(cd) isa(cd, 'AnalogChannelDescriptor'), cdList);
            chListIndiv = chList(maskIndiv);
            cdIndiv = cdList(maskIndiv);
            inGroup = arrayfun(@(cd) cd.isColumnOfSharedMatrix, cdIndiv);
            groupList = unique(arrayfun(@(cd) cd.primaryDataField, cdIndiv(inGroup), 'UniformOutput', false));
            indivChList = chListIndiv(~inGroup);

            maskGroup = arrayfun(@(cd) isa(cd, 'AnalogChannelGroupDescriptor'), cdList);
            groupList = union(groupList, chList(maskGroup));
            cdOther = cdList(~maskIndiv & ~maskGroup);
            if ~isempty(cdOther)
                warning('Trimming analog time field %s referenced by other channels %s', ...
                    timeField, TrialDataUtilities.String.strjoin({cdOther.name}));
            end

            prog = ProgressBar(numel(indivChList) + numel(groupList), 'Stripping analog channels sharing time field %s', timeField);

            % overwrite time field
            times(notEmpty) = cellfun(@(t, m) t(m), times(notEmpty), timesMask(notEmpty), 'UniformOutput', false);
            td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, timeField, times);

            for iC = 1:numel(indivChList)
                prog.increment();
                samples = {td.data.(indivChList{iC})}';
                samplesNotEmpty = ~cellfun(@isempty, samples) & notEmpty;
                samples(samplesNotEmpty) = cellfun(@(v, m) v(m), samples(samplesNotEmpty), timesMask(samplesNotEmpty), 'UniformOutput', false);
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, indivChList{iC}, samples);
            end
            for iG = 1:numel(groupList)
                prog.increment();
                samples = {td.data.(groupList{iG})}';
                samplesNotEmpty = ~cellfun(@isempty, samples) & notEmpty;
                samples(samplesNotEmpty) = cellfun(@(v, m) v(m, :), samples(samplesNotEmpty), timesMask(samplesNotEmpty), 'UniformOutput', false);
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, groupList{iG}, samples);
            end
            prog.finish();

            td = td.postDataChange(cat(1, indivChList, groupList));
        end

        function [times, values] = trimIncomingAnalogChannelData(td, times, values, varargin)
            % delete time points outside of a certain start stop interval
            % defaults to TrialStart:TrialStop
            p = inputParser();
            p.addParameter('startTimes', {}, @isvector);
            p.addParameter('stopTimes', {}, @isvector);
            p.addParameter('clearInvalidTrials', true, @islogical);
            p.addParameter('isAligned', true, @islogical);
            p.addParameter('zeroTimes', {}, @isvector);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            if isvector(times) && ~iscell(times)
                times = repmat({times}, td.nTrials, 1);
            end
            if isnumeric(values) || islogical(values)
                assert(size(values, 1) == td.nTrials, 'Values must be nTrials along dimension 1');

                valNum = values;
                values = cellvec(td.nTrials);
                for r = 1:td.nTrials
                    values{r} = TensorUtils.squeezeDims(valNum(r, :, :), 1);
                end
            end
            notEmpty = ~cellfun(@isempty, times);

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
            if isempty(p.Results.zeroTimes)
                zeroTimes = td.getTimeOffsetsFromZeroEachTrial();
            else
                zeroTimes = p.Results.zeroTimes;
            end

            if p.Results.isAligned
                offsets = zeroTimes;
            else
                offsets = zeros(size(zeroTimes));
            end

            if p.Results.clearInvalidTrials
                startTimes(~td.valid) = NaN;
                stopTimes(~td.valid) = NaN;
                offsets(~td.valid) = NaN;
            end

            timesMask = cellvec(td.nTrials);
            for iT = 1:td.nTrials
                timesMask{iT} = times{iT} + offsets(iT) >= startTimes(iT) & times{iT} + offsets(iT) <= stopTimes(iT);
            end

            maskNeedsModification = ~cellfun(@all, timesMask);
            if ~any(maskNeedsModification)
                return;
            end

            % overwrite time field
            times(notEmpty) = cellfun(@(t, m) t(m), times(notEmpty), timesMask(notEmpty), 'UniformOutput', false);

            samplesNotEmpty = ~cellfun(@isempty, values) & notEmpty;
            values(samplesNotEmpty) = cellfun(@(v, m) v(m, :, :, :, :, :, :, :), values(samplesNotEmpty), timesMask(samplesNotEmpty), 'UniformOutput', false);
        end

        function td = incorporateAnalogChannelsIntoGroup(td, subNames, groupName, varargin)
            p = inputParser();
            p.addParameter('channelDescriptor', [], @(x) isempty(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('sampleSize', numel(subNames), @isvector);

            p.KeepUnmatched = true;
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            cdsSub = td.getChannelDescriptor(subNames);
            cde = cdsSub(1);
            [data, time] = td.getAnalogMultiCommonTime(subNames, 'raw', true);
            td = td.addAnalogChannelGroup(groupName, data, time, 'raw', true, ...
                'isAligned', false, 'scaleFromLims', cde.scaleFromLims, 'scaleToLims', cde.scaleToLims, ...
                p.Unmatched);

            td = td.dropChannels(subNames);
        end
    end

    methods % Continuous Neural Channels
        function td = addContinuousNeuralChannelGroup(td, groupName, electrodes, varargin)
            % see addAnalogChannelGroup, same signature
            td.warnIfNoArgOut(nargout);
            td = td.addAnalogChannelGroup(groupName, 'isContinuousNeural', true, ...
                'continuousNeuralElectrodes', electrodes, varargin{:});
        end

        function td = addOrUpdateContinuousNeuralChannelGroup(td, groupName, electrodes, data, times, varargin)
            td.warnIfNoArgOut(nargout);

            if td.hasAnalogChannelGroup(groupName)
                cd = td.getChannelDescriptor(groupName);
                assert(isa(cd, 'ContinuousNeuralChannelGroupDescriptor'));
                assert(isequal(makecol(electrodes), cd.electrodes));
            end

            td = td.addOrUpdateAnalogChannelGroup(groupName, data, times, 'isContinuousNeural', true, 'continuousNeuralElectrodes', electrodes, varargin{:});
        end

        function [names, channelDescriptors] = listContinuousNeuralChannels(td, varargin)
            p = inputParser();
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.addParameter('includeNonGroupSubChannels', true, @islogical); % include those not in explicit arrays
            p.addParameter('includeGroupSubChannels', true, @islogical); % include those in explicit arrays
            p.addParameter('type', '', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.addParameter('array', '', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.addParameter('electrode', [], @(x) isvector(x));
            p.parse(varargin{:});

            % this should grab array sub channels if includeArraySubChannels
            [names, channelDescriptors] = td.listChannels('includeNamedSubChannels', p.Results.includeGroupSubChannels, ....
                'channelDescriptorClass', ["ContinuousNeuralChannelGroupDescriptor", "ContinuousNeuralChannelDescriptor"], 'includeDerivedChannelTypes', true);
            mask = arrayfun(@(cd) isa(cd, 'ContinuousNeuralChannelDescriptor'), channelDescriptors);
            names = names(mask);
            channelDescriptors = channelDescriptors(mask);

            mask = true(numel(names), 1);
            if ~p.Results.includeNonGroupSubChannels
                mask = mask & arrayfun(@(cd) cd.isColumnOfSharedMatrix, channelDescriptors);
            end

            if ~isempty(p.Results.type)
                type = arrayfun(@(cd) string(cd.type), channelDescriptors);
                mask = mask & ismember(type, string(p.Results.type));
            end
            if ~isempty(p.Results.array)
                array = arrayfun(@(cd) string(cd.array), channelDescriptors);
                mask = mask & ismember(array, string(p.Results.array));
            end
            if ~isempty(p.Results.electrode)
                electrode = arrayfun(@(cd) cd.electrode, channelDescriptors);
                mask = mask & ismember(electrode, p.Results.electrode);
            end

            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

        function [names, channelDescriptors] = listContinuousNeuralChannelsOnArray(td, arrayName, varargin)
            % [names, channelDescriptors] = getContinuousNeuralChannelsOnArray(td, arrayName)
            [names, channelDescriptors] = td.listContinuousNeuralChannels('array', arrayName, varargin{:});
        end

        function [names, channelDescriptors] = listContinuousNeuralChannelsOnArrayElectrode(td, arrayName, electrodeNum, varargin)
            [names, channelDescriptors] = td.listContinuousNeuralChannels(td, 'array', arrayName, 'electrode', electrodeNum, varargin{:});
        end

        function [names] = listContinuousNeuralChannelsOnSameArrayElectrodeAs(td, chName)
            % [names, units] = listContinuousNeuralChannelsOnSameArrayElectrodeAs(td, chName)
            % chName is spike channel or continuous neural channel name
            cd = td.channelDescriptorsByName.(chName);
            [names] = td.listContinuousNeuralChannelsOnArrayElectrode(cd.array, cd.electrode);
        end

        function [groupNames, channelDescriptors, channelsByGroup] = listContinuousNeuralChannelGroups(td, varargin)
            p = inputParser();
            p.addParameter('channelDescriptorClass', 'ContinuousNeuralChannelGroupDescriptor', @ischar);
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.parse(varargin{:});

            [groupNames, channelDescriptors, channelsByGroup] = td.listChannelGroups(p.Results);
        end

        function [types, arrays] = listContinuousNeuralTypeArrays(td)
            [t1, a1] = td.listExplicitContinuousNeuralArrays();
            [t2, a2] = td.listImplicitContinuousNeuralArrays();
            types = cat(1, t1, t2);
            arrays = cat(1, a1, a2);
        end

        function [types, arrays] = listImplicitContinuousNeuralArrays(td)
            [~, cds] = td.listContinuousNeuralChannels('includeGroupSubChannels', false);
            types = arrayfun(@(cd) string(cd.type), cds);
            arrays = arrayfun(@(cd) string(cd.array), cds);
            rows = unique([types, arrays], 'rows');
            types = rows(:, 1);
            arrays = rows(:, 2);
        end

        function [types, arrays, names] = listExplicitContinuousNeuralArrays(td)
            [~, cds] = td.listContinuousNeuralChannelGroups();
            types = arrayfun(@(cd) string(cd.type), cds);
            arrays = arrayfun(@(cd) string(cd.array), cds);
            [rows, idx] = sortrows([types, arrays]);
            types = rows(:, 1);
            arrays = rows(:, 2);
            names = arrayfun(@(cd) string(cd.name), cds(idx));
        end

        function td = incorporateContinuousNeuralChannelsIntoGroups(td)
            td.warnIfNoArgOut(nargout);
            [types, arrays] = td.listImplicitContinuousNeuralArrays();

            [explicitTypes, explicitArrays] = td.listExplicitContinuousNeuralArrays();
            implicitKey = arrayfun(@strcat, types, arrays);
            explicitKey = arrayfun(@strcat, explicitTypes, explicitArrays);

            bad = ismember(implicitKey, explicitKey);
            assert(~any(bad), 'Failed. Explicit spike type_array(s) %s already exists', strjoin(implicitKey(bad), ','));

            prog = ProgressBar(numel(arrays), 'Incorporating continuous neural channels');
            for iA = 1:numel(arrays)
                prog.update(iA, 'Incorporating continuous neural array %s_%s channels', types{iA}, arrays{iA});
                % find elec/unit combos in this array
                [subNames, cdsSub] = td.listContinuousNeuralChannels('type', types(iA), 'array', arrays{iA});
                electrodes = arrayfun(@(cd) cd.electrode, cdsSub);
                units = cdsSub(1).unitsPrimary;
                groupName = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArray(types{iA}, arrays{iA});
                td = td.incorporateAnalogChannelsIntoGroup(subNames, groupName, 'isContinuousNeural', true, ...
                    'continuousNeuralElectrodes', electrodes, 'units', units);
            end
            prog.finish();
        end

        function td = sortChannelsInContinuousNeuralChannelGroup(td, group)
            td.warnIfNoArgOut(nargout);
            cd = td.getChannelDescriptor(group);
            [electrodes, order] = sort(cd.electrodes);
            td = td.filterColumnsAnalogChannelGroup(group, order);

            cd = td.getChannelDescriptor(group); % request again just in case it has been modified
            cd.electrodes = electrodes;
            td = td.setChannelDescriptor(group, cd);
        end

        function td = remapContinuousNeuralChannelGroupElectrodes(td, arrayChName, new_electrodes, varargin)
            p = inputParser();
            p.addParameter('resort', false, @islogical);
            p.parse(varargin{:});
            % works for spike arrays and continuous neural channel groups
            td.warnIfNoArgOut(nargout);

            cd = td.getChannelDescriptor(arrayChName);
            assert(isa(cd, 'ContinuousNeuralChannelGroupDescriptor'));

            nElectrodes = numel(cd.electrodes);
            assert(nElectrodes == numel(new_electrodes));
            cd.electrodes = new_electrodes;
            td = td.setChannelDescriptor(arrayChName, cd);
            if p.Results.resort
                td = td.sortChannelsInContinuousNeuralChannelGroup(arrayChName);
            end
        end

        function td = setContinuousNeuralChannelArray(td, contCh, array)
            td.warnIfNoArgOut(nargout);
            contCh = string(contCh);
            td.assertHasChannel(contCh, false); % non-sub channels only for rename

            for iCh = 1:numel(contCh)
                cd = td.getChannelDescriptor(contCh{iCh});
                newName = cd.getNameWithUpdatedArray(array);
                td = td.renameChannel(contCh{iCh}, newName);
            end
        end

        function td = renameContinuousNeuralChannelArray(td, arrayCurrent, arrayNew)
            td.warnIfNoArgOut(nargout);
            if td.hasChannel(arrayCurrent, false)
            end
            contChList = td.listContinuousNeuralChannelsOnArray(arrayCurrent);
            td = td.setContinuousNeuralChannelArray(contChList, arrayNew);
        end

        function td = setContinuousNeuralChannelType(td, contCh, type)
            td.warnIfNoArgOut(nargout);
            contCh = TrialDataUtilities.Data.wrapCell(contCh);

            for iCh = 1:numel(contCh)
                cd = td.getChannelDescriptor(contCh{iCh});
                newName = cd.getNameWithUpdatedType(type);
                td = td.renameChannel(contCh{iCh}, newName);
            end
        end

    end

    methods % Image channel methods - mostly defer to AnalogChannelGroupMethods
        function [groupNames, channelDescriptors, channelsByGroup] = listImageChannels(td, varargin)
            p = inputParser();
            p.addParameter('channelDescriptorClass', 'ImageChannelDescriptor', @ischar);
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.parse(varargin{:});

            [groupNames, channelDescriptors, channelsByGroup] = td.listChannelGroups(p.Results);
        end

        function td = addImage(td, groupName, varargin)
            % see , same signature
            td.warnIfNoArgOut(nargout);
            td = td.addAnalogChannelGroup(groupName, varargin{:}, 'isImage', true);
        end

        function td = addOrUpdateImage(td, groupName, data, times, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addOrUpdateAnalogChannelGroup(groupName, data, times, varargin{:}, 'isImage', true);
        end

    end

    methods % Event channel methods
        function td = addEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('name', @TrialDataUtilities.String.isstringlike);
            p.addRequired('times', @(x) isempty(x) || isvector(x) || ismatrix(x));
            p.addParameter('isAligned', true, @islogical);
            p.addParameter('useExistingDataField', false, @islogical);
            p.addParameter('color', [], @(x) true);
            %p.addParamValue('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.addParameter('displayGroup', '', @ischar);

            p.parse(name, times, varargin{:});
            name = char(name);
            %cd = p.Results.channelDescriptor;

            if ~p.Results.useExistingDataField
                if isempty(times)
                    times = cellvec(td.nTrials);
                end
                if isscalar(times)
                    times = repmat(times, td.nTrials, 1);
                end

                if ismatrix(times) && ~isvector(times)
                    % assume each row is one trial and convert into cell
                    % ignoring nans
                    assert(size(times, 1) == td.nTrials);
                    times = TensorUtils.splitAlongDimension(times, 1, 1);
                    times = cellfun(@(x) makecol(removenan(x)), times, 'UniformOutput', false);
                end

                assert(numel(times) == td.nTrials, 'Times must be vector with length %d', td.nTrials);
            else
                times = {td.data.(name)}';
            end
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

            cd.color = p.Results.color;
            cd.displayGroup = p.Results.displayGroup;

            if p.Results.useExistingDataField
                td = td.addChannel(cd, {}, 'ignoreDataFields', true, 'ignoreExisting', true);
            else
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
        end

        function td = addOrUpdateEvent(td, name, times, varargin)
            % set values of channel name if it exists where mask is true.
            % By default mask is non-nan values of mask
            % otherwise create channel

            td.warnIfNoArgOut(nargout);
            if td.hasEventChannel(name)
                td = td.setEvent(name, times);
            else
                td = td.addEvent(name, times);
            end
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

        function [td, eventName] = addEventFromAnalogTimes(td, name)
            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelOrGroup(name);

            timeField = td.getAnalogTimeField(name);
            if ~td.hasEventChannel(timeField)
                % add an event channel in situ with the same name as the
                % time field
                td = td.addEvent(timeField, {}, 'useExistingDataField', true);
            end
            eventName = timeField;
        end

        function tf = hasEventChannel(td, name)
            tf = td.hasChannel(name, 'includeNamedSubChannels', false, ...
                'channelDescriptorClass', "EventChannelDescriptor");
        end

        function td = setEventColor(td, name, color)
            td.warnIfNoArgOut(nargout);
            if ~td.hasEventChannel(name)
                error('Event %s not found', name);
            end
            td.channelDescriptorsByName.(name).color = color;
        end

        function [names, channelDescriptors] = listEventChannels(td, varargin)
            p = inputParser();
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.parse(varargin{:});

            [names, channelDescriptors] = td.listChannels('includeNamedSubChannels', false, ...
                'channelDescriptorClass', 'EventChannelDescriptor', ...
                'includeDerivedChannelTypes', p.Results.includeDerivedChannelTypes);

%             mask = true(numel(names), 1);

%             names = names(mask);
%             channelDescriptors = channelDescriptors(mask);
        end

        function [displayGroups, channelsByGroup, channelsNotInGroup] = listEventChannelDisplayGroups(td)
            [names, channelDescriptors] = td.listEventChannels();
            displayGroupsByCh = TrialDataUtilities.String.stringarrayfun(@(cd) cd.displayGroup, channelDescriptors);

            [displayGroups, ~, whichGroup] = unique(displayGroupsByCh);
            mask = arrayfun(@(grp) grp ~= "" || ismissing(grp), displayGroups);
            channelsByGroup = arrayfun(@(gidx) names(whichGroup==gidx), (1:numel(displayGroups))', 'UniformOutput', false);

            channelsNotInGroup = cat(1, channelsByGroup{~mask});
            displayGroups = displayGroups(mask);
            channelsByGroup = channelsByGroup(mask);

            [displayGroups, isort] = sort(displayGroups);
            channelsByGroup = channelsByGroup(isort);
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
                if td.channelDescriptorsByName.(ch).isVectorizableByField(1)
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
            if td.hasEventChannel(name) || td.hasSpikeChannel(name)
                field = name;
                dataFieldIdx = 1;
            elseif td.hasAnalogChannelOrGroup(name)
                field = td.getAnalogTimeField(name);
                dataFieldIdx = 2;
            else
                error('Unknown event %s', name);
            end
            timesCell = {td.data.(field)}';

            % remove nans from the event list
            timesCell = cellfun(@(x) x(~isnan(x)), timesCell, 'UniformOutput', false);

            % convert to access class
            cd = td.channelDescriptorsByName.(name);
            impl = cd.getImpl();
            timesCell = impl.convertDataCellOnAccess(dataFieldIdx, timesCell);
        end

        function times = getEventRawFirst(td, name)
            timesCell = td.getEventRaw(name);

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
            timesCell = td.getEventRaw(name);

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
        
        function counts = getEventCountRaw(td, name)
            counts = cellfun(@numel, td.getEventRaw(name));
        end

        function tfList = getEventOccurred(td, name)
            tfList = ~cellfun(@isempty, td.getEvent(name));
        end
        
        function tfList = getEventOccurredRaw(td, name)
            tfList = ~cellfun(@isempty, td.getEventRaw(name));
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
        
        function times = getEventNthRaw(td, name, n)
            if strcmp(n, 'end')
                times = td.getEventLastRaw(name);
                return;
            end

            times = cellfun(@getNth, td.getEventRaw(name));

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
        
        function times = getEventFirstRaw(td, name)
            times = getEventNthRaw(td, name, 1);
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
        
        function times = getEventLastRaw(td, name)
            times = cellfun(@getLast, td.getEventRaw(name));

            function t = getLast(times)
                if ~isempty(times)
                    t = times(end);
                else
                    t = NaN;
                end
            end
        end

        function [td, fieldsUpdated] = setEvent(td, name, times, varargin)
            p = inputParser();
            p.addOptional('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('deferPostDataChange', false, @islogical);
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
                if ismatrix(times) && ~isvector(times)
                    % assume each row is one trial and convert into cell
                    % ignoring nans
                    assert(size(times, 1) == td.nTrials);
                    times = TensorUtils.splitAlongDimension(times, 1, 1);
                    times = cellfun(@(x) makecol(removenan(x)), times, 'UniformOutput', false);
                else
                    times = arrayfun(@removenan, times, 'UniformOutput', false);
                end
            end
            times = cellfun(@(x) makecol(sort(x)), times, 'UniformOutput', false);

            [td, ~, fieldsUpdated] = td.setChannelData(name, {times}, 'deferPostDataChange', p.Results.deferPostDataChange);
        end

        function td = addOrSetEvent(td, name, times, varargin)
            td.warnIfNoArgOut(nargout);
            if td.hasEventChannel(name)
                td = td.setEvent(name, times, varargin{:});
            else
                td = td.addEvent(name, times, varargin{:});
            end
        end

        function td = trimEventRaw(td, names, varargin)
            % delete time points outside of a certain start stop interval
            % defaults to TrialStart:TrialStop
            p = inputParser();
            p.addOptional('startTimes', {}, @isvector);
            p.addOptional('stopTimes', {}, @isvector);
            p.addParameter('clearInvalidTrials', true, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            if ischar(names)
                names = {names};
            end

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
            if p.Results.clearInvalidTrials
                startTimes(~td.valid) = NaN;
                stopTimes(~td.valid) = NaN;
            end

            prog = ProgressBar(numel(names), 'Trimming event channels');
            fieldsUpdated = cellvec(numel(names));
            for i = 1:numel(names)
                name = names{i};
                prog.update(i, 'Trimming event %s', name);
                assert(td.hasEventChannel(name), 'No event channel named %s', name);
                times = td.getEventRaw(name);
                %                 times = {td.data.(name)}';

                timesMask = cellvec(td.nTrials);
                for iT = 1:td.nTrials
                    timesMask{iT} = times{iT} >= startTimes(iT) & times{iT} <= stopTimes(iT);
                end

                maskNeedsModification = ~cellfun(@all, timesMask);
                if ~any(maskNeedsModification)
                    continue;
                end

                notEmpty = ~cellfun(@isempty, times);
                times(notEmpty) = cellfun(@(t, m) t(m), times(notEmpty), timesMask(notEmpty), 'UniformOutput', false);

                [td, fieldsUpdated{i}] = td.setEvent(name, times, ...
                    'isAligned', false, 'deferPostDataChange', true);
                %                 td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, name, times);
            end
            % defer post data change to the end so that align info need
            % only be applied once
            td = td.postDataChange(cat(1, fieldsUpdated{:}));
            prog.finish();
        end

        function td = trimEventToTrialStartEnd(td, names)
            % Timepoints that lie outside of TrialStart and TrialStop will
            % never be accessible via getTimes since they will be filtered
            % out by the AlignInfo

            td.warnIfNoArgOut(nargout);
            % default is TrialStart and TrialEnd, so just pass it along
            td = td.trimEventRaw(names);
        end

    end

    methods % Param channel methods
        function td = addParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('values', '', @(x) true);
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'ChannelDescriptor'));
            p.addParameter('like', '', @ischar);
            p.addParameter('units', '', @ischar);
            p.addParameter('displayGroup', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            %name = p.Results.name;
            name = char(name);
            values = p.Results.values;

            if ~isempty(values)
                % expand scalar values to be nTrials x 1
                if ischar(values)
                    values = repmat({values}, td.nTrials, 1);
                elseif numel(values) == 1
                    values = repmat(values, td.nTrials, 1);
                else
                    values = makecol(values);
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
            cd.displayGroup = p.Results.displayGroup;

            if ~ismember('units', p.UsingDefaults)
                cd = cd.setPrimaryUnits(p.Results.units);
            end

            td = td.addChannel(cd, {values}, p.Unmatched);
        end

        function td = addScalarParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.addParameter('units', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            values = p.Results.values;
            cd = ParamChannelDescriptor.buildScalarParam(name, class(values), p.Results.units);
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end

        function td = addVectorParamAccessAsMatrix(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addOptional('values', {}, @(x) isempty(x) || ismatrix(x) || iscell(x));
            p.addParameter('units', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            values = p.Results.values;
            if iscell(values)
                dataClass = ChannelImpl.getCellElementClass(values);
            else
                dataClass = class(values);
            end
            cd = ParamChannelDescriptor.buildVectorParamAccessAsMatrix(name, dataClass, p.Results.units);

            if isnumeric(values)
                values = TensorUtils.splitAlongDimension(values, 1);
            end
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end

        function td = addOrUpdateScalarParam(td, name, vals, varargin)
            % set values of channel name if it exists where mask is true.
            % By default mask is non-nan values of mask
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', ~isnan(vals), @isvector);
            p.addParameter('units', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;
            td.warnIfNoArgOut(nargout);
            if td.hasParamChannel(name)
                td = td.setParam(name, vals, 'updateMask', mask);
            else
                vals(~mask) = NaN;
                td = td.addScalarParam(name, vals, 'units', p.Results.units);
            end
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

        function td = addOrUpdateStringParam(td, name, vals, varargin)
            % set values of channel name if it exists where mask is true.
            % By default mask is non-empty values
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', ~cellfun(@isempty, vals), @isvector);
            p.parse(varargin{:});

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;
            td.warnIfNoArgOut(nargout);
            if td.hasParamChannel(name)
                td = td.setParam(name, vals, 'updateMask', mask);
            else
                [vals{~mask}] = deal([]);
                td = td.addStringParam(name, vals);
            end
        end

        function td = addBooleanParam(td, name, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addOptional('values', {}, @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            name = char(name);
            cd = ParamChannelDescriptor.buildBooleanParam(name);
            values = p.Results.values;
            td = td.addParam(name, values, 'channelDescriptor', cd, ...
                p.Unmatched);
        end

        function td = addOrUpdateBooleanParam(td, name, vals, varargin)
            % set values of channel name if it exists where mask is true.
            % By default mask is non-empty values
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', vals, @isvector);
            p.parse(varargin{:});

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;

            td.warnIfNoArgOut(nargout);
            if td.hasParamChannel(name)
                td = td.setParam(name, vals, 'updateMask', mask);
            else
                vals(~mask) = false;
                td = td.addBooleanParam(name, vals);
            end
        end

        function tf = hasParamChannel(td, name)
            tf = td.hasChannel(name, 'includeNamedSubChannels', false, ...
                'channelDescriptorClass', "ParamChannelDescriptor");
        end

        function [names, channelDescriptors] = listParamChannels(td, varargin)
            p = inputParser();
            p.addParameter('includeDerivedChannelTypes', true, @islogical);
            p.addParameter('isScalar', false, @islogical);
            p.addParameter('isString', false, @islogical);
            p.addParameter('isCategorical', false, @islogical);
            p.parse(varargin{:});

            [names, channelDescriptors] = td.listChannels('includeNamedSubChannels', false, ...
                'channelDescriptorClass', 'ParamChannelDescriptor', ...
                'includeDerivedChannelTypes', p.Results.includeDerivedChannelTypes);

            mask = true(numel(names), 1);
            if p.Results.isScalar
                mask = mask & arrayfun(@(cd) cd.isVectorizableByField(1), channelDescriptors);
            end
            if p.Results.isString
                mask = mask & arrayfun(@(cd) cd.isStringByField(1), channelDescriptors);
            end
            if p.Results.isCategorical
                mask = mask & arrayfun(@(cd) cd.isCategoricalByField(1), channelDescriptors);
            end

            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

        function [names, channelDescriptors] = listScalarParamChannels(td, varargin)
            [names, channelDescriptors] = td.listParamChannels('isScalar', true, varargin{:});
        end

        function [names, channelDescriptors]  = listStringParamChannels(td, varargin)
            [names, channelDescriptors] = td.listParamChannels('isString', true, varargin{:});
        end

        function [displayGroups, channelsByGroup, channelsNotInGroup] = listParamChannelDisplayGroups(td)
            [names, channelDescriptors] = td.listParamChannels();
            displayGroupsByCh = TrialDataUtilities.String.stringarrayfun(@(cd) cd.displayGroup, channelDescriptors);

            [displayGroups, ~, whichGroup] = unique(displayGroupsByCh);
            mask = arrayfun(@(grp) grp ~= "", displayGroups);
            channelsByGroup = arrayfun(@(gidx) names(whichGroup==gidx), (1:numel(displayGroups))', 'UniformOutput', false);

            channelsNotInGroup = cat(1, channelsByGroup{~mask});
            displayGroups = displayGroups(mask);
            channelsByGroup = channelsByGroup(mask);

            [displayGroups, isort] = sort(displayGroups);
            channelsByGroup = channelsByGroup(isort);
        end

        function [ch, idx] = parseIndexedParamChannelName(td, name) %#ok<INUSL>
            % if name is a normal param channel, do nothing
            % if it is a name like param(5), return param, 5

            tokens = regexp(name, '(?<ch>[\w_]+)\((?<idx>\d+)\)', 'names');
            if isempty(tokens)
                ch = name;
                idx = NaN;
            else
                ch = tokens.ch;
                idx = str2double(tokens.idx);
            end
        end

        function cd = getParamChannelDescriptor(td, name)
            chName = td.parseIndexedParamChannelName(name);
            cd = td.channelDescriptorsByName.(chName);
        end

        function values = getParamRaw(td, name)
            % grab the raw value for a parameter without considering
            % validity. if name is an analog channel, grab the first sample
            % of that channel

            [chName, colIdx] = td.parseIndexedParamChannelName(name);
            cd = td.channelDescriptorsByName.(chName);
            impl = cd.getImpl();

            if isa(cd, 'ParamChannelDescriptor')
                values = {td.data.(cd.dataFields{1})}';
                values = impl.convertDataCellOnAccess(1, values); % convert to access data class

                if ~isnan(colIdx)
                    values = values(:, colIdx);
                end
            elseif isa(cd, 'AnalogChannelDescriptor')
                values = td.getAnalogSample(chName);
            elseif isa(cd, 'EventChannelDescriptor')
                values = td.getEventFirst(chName);
            else
                error('Only valid for parameter or analog channels');
            end
        end

        % Basic access methods, very fast
        function values = getParam(td, name)
            [chName, ~] = td.parseIndexedParamChannelName(name);
            cd = td.channelDescriptorsByName.(chName);
            values = td.getParamRaw(name);
            values = td.replaceInvalidMaskWithValue(values, cd.missingValueByField{1});
        end

        function values = getParamValid(td, name)
            values = td.getParam(name);
            values = values(td.valid);
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
            names = string(names);
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

            names = string(names);
            valueCell = td.getParamMultiAsCell(names);
            units = arrayfun(@(name) td.getChannelUnitsPrimary(name), names);

            if p.Results.includeValidColumn
                valueCell = horzcat( num2cell(td.valid), valueCell );
                names = vertcat("valid", makecol(names));
            end

            numWidth = ceil(log(td.nTrials)/log(10));
            trialNames = arrayfun(@(ind) sprintf('trial%0*d', numWidth, ind), (1:td.nTrials)', 'UniformOutput', false);

            valueTable = cell2table(valueCell, 'VariableNames', names, ...
                'RowNames', trialNames);
            valueTable.Properties.VariableUnits = units;
        end

        function uniqTable = getParamMultiUniqueAsTable(td, names, varargin)
            p = inputParser();
            p.addParameter('includeCounts', true, @islogical);
            p.parse(varargin{:});

            valueTable = td.getParamMultiAsTable(names);
            valueTable = valueTable(td.valid, :);

            [uniqTable, ~, which] = TrialDataUtilities.Data.uniqueSingleNaN(valueTable);
            uniqTable.Properties.RowNames = {};

            counts = hist(which, 1:max(which))'; %#ok<HIST>
            uniqTable.TrialCount = counts;
            uniqTable.Properties.VariableUnits{end} = 'trials';
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
            if ~iscell(vals)
                if isstring(vals)
                    if any(ismissing(vals))
                        values = [unique(vals); string(missing)];
                    else
                        values = unique(vals);
                    end
                elseif islogical(vals)
                    values = unique(vals);
                else
                    values = TrialDataUtilities.Data.uniqueTolSingleNaN(vals);
                end
            end
                
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
            if isnumeric(vals) && cd.collectAsCellByField(1)
                % provided as matrix but must be split along dim 1
                assert(size(vals, 1) == td.nTrials, 'Size of matrix along dim 1 must equal nTrials');
                vals = TensorUtils.splitAlongDimension(vals, 1);
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

            missing = cd.missingValueByField{1};
            if iscell(valuesCurrent)
                filler = repmat({missing}, n, 1);
                values = cat(1, filler, valuesCurrent(1:(td.nTrials-n)));
            else
                sz = size(valuesCurrent);
                sz(1) = n;
                filler = repmat(missing, sz);
                values = cat(1, filler, TensorUtils.selectAlongDimension(valuesCurrent, 1, 1:(td.nTrials-n)));
            end

            % flush values on currently invalid trials
            values = td.replaceInvalidMaskWithValue(values, cd.missingValueByField{1});
        end

        function td = addParamNBack(td, name, n, varargin)
            p = inputParser;
            p.addParameter('as', '', @ischar); % new channel name
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            %cd = td.channelDescriptorsByName.(name);

            if isempty(p.Results.as)
                as = [name sprintf('_%dback', n)];
            else
                as = p.Results.as;
            end

            values = td.getParamNBack(name, n);
            td = td.addParam(as, values, 'like', name);
        end
    end

    methods % Spike channel methods
        function td = addSpikeChannel(td, unitStr, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addOptional('spikes', {}, @isvector);
            p.addParameter('isAligned', true, @isscalar);
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsTime', [], @isvector); % common time vector to be shared for ALL waveforms for this channel
            p.addParameter('waveformsField', sprintf('%s_waveforms', unitStr), @ischar);
            p.addParameter('waveformsScaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsScaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is
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
            spikes = cellfun(@(x, y) makecol(x+y), p.Results.spikes, num2cell(offsets), 'UniformOutput', false);

            channelData = {spikes};

            if ~isempty(p.Results.waveforms)
                waveforms = p.Results.waveforms;

                waveClass = TrialDataUtilities.Data.cellDetermineCommonClass(waveforms);
                cd = cd.addWaveformsField(p.Results.waveformsField, 'time', p.Results.waveformsTime, ...
                    'dataClass', waveClass, 'scaleFromLims', p.Results.waveformsScaleFromLims, ...
                    'scaleToLims', p.Results.waveformsScaleToLims);

                % undo scaling and convert back to memory scale
                if ~p.Results.waveformsInMemoryScale
                    impl = cd.getImpl();
                    waveforms = impl.unscaleWaveforms(waveforms);
                end

                channelData{end+1} = makecol(waveforms);
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
        
        function n = getSpikeArrayChannelCount(td, array)
            if td.hasExplicitSpikeArray(array)
                cd = td.getChannelDescriptor(array);
                n = cd.nChannels;
            else
                % implicit array
                subNames = td.listSpikeChannels('array', array);
                n = numel(subNames);
            end
        end
        
        function [electrodes, units] = listSpikeElectrodesUnitsOnArray(td, array)
            if td.hasExplicitSpikeArray(array)
                cd = td.getChannelDescriptor(array);
                electrodes = cd.electrodes;
                units = cd.units;
            else
                % implicit array
                subNames = td.listSpikeChannels('array', array);
                [~, electrodes, units] = arrayfun(@(subname) SpikeChannelDescriptor.parseArrayElectrodeUnit(subname), string(subNames));
            end
        end

        function td = setSpikeChannel(td, name, times, varargin)
            % td = setSpikeChannel(td, name, times, varargin)
            % options:
            %  isAligned
            %  waveforms
            p = inputParser();
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('blankingRegions', {}, @(x) isempty(x) || iscell(x)); % nTrials x 1 cell of nIntervals x 2 matrices
            p.addParameter('storeBlankingRegions', true, @islogical); % actually mark the blanking regions so that they will
            p.addParameter('sortQualityByTrial', [], @(x) isempty(x) || isvector(x)); % nTrials x 1 vector of per-trial ratings
            p.addParameter('waveformsInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is

            p.KeepUnmatched = true; % we don't want to capture a fieldMask because this would necessitate thinking about the logic of whether to apply the blanking
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            td.assertHasChannel(name);
            cd = td.channelDescriptorsByName.(name);
            assert(isa(cd, 'SpikeChannelDescriptor') || isa(cd, 'SpikeArrayChannelDescriptor'));

            if p.Results.isAligned
                % add the zero offset to the time vector for each trial
                % this is mostly for TDCA, so that alignments info is
                % preserved
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end

            if cd.nChannels == 1
                times = makecol(times);
            end
            if iscell(times)
                times = cellfun(@plus, times, repmat(num2cell(offsets), 1, cd.nChannels), 'UniformOutput', false);
            else
                times = times + offsets;
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

                % undo scaling and convert back to memory scale
                if ~p.Results.waveformsInMemoryScale
                    impl = cd.getImpl();
                    newWaves = impl.unscaleWaveforms(newWaves);
                end

            else
                % no waveforms provided
                if cd.hasWaveforms
                    % no waveforms provided, but channel has waveforms
                    % check that number of spikes isn't changing so that
                    % the correspondence is maintained
                    nSpikesProvided = cellfun(@numel, times);
                    newWaves = td.getRawSpikeWaveforms(name);
                    nSpikesWave = cellfun(@(w) size(w, 1), newWaves);
                    assert(isequaln(nSpikesProvided, nSpikesWave), 'Number of spikes in each trial cannot change unless waveforms are also provided');

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

            % start with empty
            if ~isempty(p.Results.blankingRegions)
                blanking = p.Results.blankingRegions;

                assert(numel(blanking) == td.nTrials);
                %assert(all(cellfun(@(x) isempty(x) || (ismatrix(x) && size(x, 2) == 2), blanking)), 'Blanking regions must be nTrials cellvec with nRegions x 2 matrices');

                % align the blanking regions with the zero offset to
                % unalign the blanking regions specified
                blanking = cellfun(@plus, blanking, num2cell(offsets), 'UniformOutput', false);

                % add field if not found (and we're storing the regions)
                if ~cd.hasBlankingRegions
                    if p.Results.storeBlankingRegions
                        cd = cd.addBlankingRegionsField();
                        td.channelDescriptorsByName.(name) = cd;
                    end

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

                if p.Results.storeBlankingRegions
                    % and add the blanking info as the last channel data
                    channelData = cat(1, channelData, {blanking});
                    channelFieldMask(end+1) = true;
                else
                    if cd.hasBlankingRegions
                        channelData = cat(1, channelData, {{}});
                        channelFieldMask(end+1) = false;
                    end
                end
            else
                if cd.hasBlankingRegions
                    channelData = cat(1, channelData, {{}});
                    channelFieldMask(end+1) = false;
                end
            end

            % make the final update
            td = td.setChannelData(name, channelData, 'fieldMask', channelFieldMask, p.Unmatched);
        end

        function td = addOrUpdateSpikeChannel(td, name, times, varargin)
            % set spike times of channel name if it exists where mask is true.
            % By default mask is non-empty cells in times
            % otherwise create channel
            p = inputParser();
            p.addParameter('mask', ~cellfun(@isempty, times), @isvector);
            p.addParameter('isAligned', true, @islogical); % time vectors reflect the current 0 or should be considered relative to TrialStart?
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsTime', [], @isvector); % common time vector to be shared for ALL waveforms for this channel
            p.addParameter('waveformsScaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsScaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('sortQuality', NaN, @isscalar); % numeric scalar metric of sort quality
            p.addParameter('sortMethod', '', @ischar);
            p.addParameter('sortQualityByTrial', [], @isvector); % nTrials x 1 vector of per-trial ratings
            p.addParameter('blankingRegions', {}, @iscell); % nTrials x 1 cell of nIntervals x 2 matrices

            p.parse(varargin{:});

            mask = TensorUtils.vectorIndicesToMask(makecol(p.Results.mask), td.nTrials) & td.valid;

            td.warnIfNoArgOut(nargout);
            if td.hasSpikeChannel(name)
                td = td.setSpikeChannel(name, times, 'updateMask', mask, ...
                    'isAligned', p.Results.isAligned, 'waveforms', p.Results.waveforms, ...
                    'blankingRegions', p.Results.blankingRegions, ...
                    'sortQualityByTrial', p.Results.sortQualityByTrial);
            else
                % clear masked out cells
                [times{~mask}] = deal([]);
                if ~isempty(p.Results.waveforms)
                    waves = p.Results.waveforms;
                    [waves{~mask}] = deal([]);
                else
                    waves = p.Results.waveforms;
                end
                td = td.addSpikeChannel(name, times, 'updateMask', mask, ...
                    'isAligned', p.Results.isAligned, 'waveforms', waves, ...
                    'blankingRegions', p.Results.blankingRegions, ...
                    'sortQualityByTrial', p.Results.sortQualityByTrial, ...
                    'waveformsTime', p.Results.waveformsTime, ...
                    'waveformsScaleFromLims', p.Results.waveformsScaleFromLims, ...
                    'waveformsScaleToLims', p.Results.waveformsScaleToLims);
            end
        end

        function td = mergeSpikeChannels(td, chList, varargin)
            % TODO update this for spike array descriptors
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

        function td = setSpikeChannelArray(td, spikeCh, array)
            td.warnIfNoArgOut(nargout);
            spikeCh = TrialDataUtilities.Data.wrapCell(spikeCh);

            for iCh = 1:numel(spikeCh)
                cd = td.getChannelDescriptor(spikeCh{iCh});
                newName = cd.getNameWithUpdatedArray(array);
                td = td.renameChannel(spikeCh{iCh}, newName);
            end
        end

        function td = renameSpikeChannelArray(td, arrayCurrent, arrayNew)
            td.warnIfNoArgOut(nargout);
            spikeChList = td.listSpikeChannels('array', arrayCurrent);
            td = td.setSpikeChannelArray(spikeChList, arrayNew);
        end

        function td = blankSpikesInTimeIntervals(td, name, intervalCell, varargin)
            % TODO update this for spike array descriptors

            % adds a blanking region to the spiking data
            % this both removes the spikes from that period of time each
            % trial, and will inform the spike rate filtering that this
            % time interval is not observed, rather than simply has no
            % spikes
            %
            % intervalCell is nTrials x 1 cell with nIntervals(iT) x 2 matrices inside

            p = inputParser();
            p.addParameter('storeBlankingRegions', true, @islogical); % mark as blanked these regions to prevent PSTHs from filtering over them
            p.addParameter('isAligned', true, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            if ischar(name)
                nameList = {name};
            else
                nameList = name;
            end

            for iU = 1:numel(nameList)
                name = nameList{iU};
                td.assertHasChannel(name);
                cd = td.channelDescriptorsByName.(name);
                assert(isa(cd, 'SpikeChannelDescriptor'));

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
                    'blankingRegions', intervalCell, 'storeBlankingRegions', p.Results.storeBlankingRegions, ...
                    'waveforms', rawWaves);
            end
        end

        function intervalCell = getRawSpikeBlankingRegions(td, unitNames, varargin)
            p = inputParser();
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});
            
            [fieldList, fieldIsArray, colIdxEachField, assignIdxEachField, nColumnsPerName] = td.getSpikeChannelMultiAccessInfo(unitNames);
            error('Need to update implementation!')

            if isnumeric(unitNames)
                unitNames = td.lookupSpikeChannelByIndex(unitNames);
            end
            unitNames = string(unitNames);
            if numel(unitNames) == 1 && td.hasExplicitSpikeArray(unitNames)
                intervalCell = td.getRawSpikeArrayBlankingRegions(unitNames(1), p.Results);
                return;
            end
            nU = numel(unitNames);
            emptyInterval = zeros(0, 2);
            
            % check if any units have blanking regions first, we can return {} to save time if none does
            channelDescriptors = td.getChannelDescriptor(unitNames);
            if ~any(arrayfun(@(cd) cd.hasBlankingRegions, channelDescriptors))
                intervalCell = {};
                return
            end
            
            intervalCell = cell(td.nTrials, numel(unitNames));
            
            for iU = 1:nU
                cd = channelDescriptors(iU);
                fld = cd.blankingRegionsField;
                if ~isempty(fld)
                    intervalCell(:, iU) = TrialDataUtilities.SpikeData.removeOverlappingIntervals(makecol({td.data.(fld)}));
                else
                    intervalCell(:, iU) = repmat({emptyInterval}, td.nTrials, 1);
                end
            end

            if p.Results.combine
                % combine the intervals in multiple units
                intervalCellByUnit = cellvec(numel(unitNames));
                for iU = 1:numel(unitNames)
                    intervalCellByUnit{iU} = intervalCell(:, iU);
                end
                intervalCell = TrialDataUtilities.SpikeData.removeOverlappingIntervals(intervalCellByUnit{:});
            end
            end
        
        function intervalCell = getRawSpikeArrayBlankingRegions(td, arrayName, varargin)
            p = inputParser();
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

%             emptyInterval = zeros(0, 2);
            cd = td.getChannelDescriptor(arrayName);
%             nU = cd.nChannels;
            
            if ~cd.hasBlankingRegions
                intervalCell = {};
%                 if p.Results.combine
%                     intervalCell = repmat({emptyInterval}, td.nTrials);
%                 else
%                     intervalCell = repmat({emptyInterval}, td.nTrials, nU);
%                 end
                return;
            end
            
            fld = cd.blankingRegionsField;
            intervalCell = TrialDataUtilities.SpikeData.removeOverlappingIntervals(makecol({td.data.(fld)}));

            if p.Results.combine
                % combine the intervals in multiple units
                intervalCellByUnit = cellvec(numel(unitNames));
                for iU = 1:numel(unitNames)
                    intervalCellByUnit{iU} = intervalCell(:, iU);
                end
                intervalCell = TrialDataUtilities.SpikeData.removeOverlappingIntervals(intervalCellByUnit{:});
            end
        end  

        function intervalCell = getSpikeBlankingRegions(td, unitName, varargin)
            intervalCell = td.getRawSpikeBlankingRegions(unitName, varargin{:});
            if ~isempty(intervalCell)
                intervalCell = td.replaceInvalidMaskWithValue(intervalCell, []);
            end
        end
        
        function nUnitsTotal = computeSpikeUnitCountFromName(td, unitNames)
            unitNames = string(unitNames);
            if isempty(unitNames)
                nUnitsTotal = 0;
            elseif strlength(unitNames) == 0
                nUnitsTotal = 0;
            else
                [~, ~, ~, ~, nColumnsPerName] = td.getSpikeChannelMultiAccessInfo(unitNames);
                nUnitsTotal = sum(nColumnsPerName);
            end
        end

        function timesCell = getRawSpikeTimes(td, unitNames, varargin)
            % timesCell is nTrials x nUnits (if unitNames is cellstr)
            % if combine is true, all spikes will be interleaved
            p = inputParser();
            p.addParameter('combine', false, @islogical);
            p.addParameter('slice', [], @(x) true); % for subselecting units from array
            p.parse(varargin{:});

            [fieldList, fieldIsArray, colIdxEachField, assignIdxEachField, nColumnsPerName] = td.getSpikeChannelMultiAccessInfo(unitNames, 'times', p.Results.slice);
            nUnits = sum(nColumnsPerName);
            nFields = numel(fieldList);
            timesCellByUnit = cell(td.nTrials, nUnits);
            for iF = 1:nFields
                fld = fieldList{iF};
                col = colIdxEachField{iF};
                idx = assignIdxEachField{iF};

                if td.nTrialsValid > 0
                    if fieldIsArray(iF)
                        catData = cat(1, td.data.(fld));
                        assert(size(catData, 1) == td.nTrials);
                    else
                        catData = {td.data.(fld)}';
                    end
                    timesCellByUnit(:, idx) = catData(:, col);
                else
                    timesCellByUnit(:, idx) = {zeros(0, numel(col))};
                end
            end

            if p.Results.combine
                timesCell = cellvec(td.nTrials);
                for iT = 1:td.nTrials
                    timesCell{iT} = cat(1, timesCellByUnit{iT, :});
                end
            else
                timesCell = timesCellByUnit;
            end
        end

        function timesCell = getSpikeTimes(td, unitNames, varargin)
            timesCell = td.getSpikeTimesUnaligned(unitNames, varargin{:});
        end

        function timesCell = getSpikeTimesUnaligned(td, unitNames, varargin)
            timesCell = td.getRawSpikeTimes(unitNames, varargin{:});
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
        end

        function [counts, hasSpikes] = getSpikeCounts(td, unitName, varargin)
            counts = cellfun(@numel, td.getSpikeTimes(unitName, varargin{:}));
            counts = td.replaceInvalidMaskWithValue(counts, NaN);
            hasSpikes = counts > 0;
        end

        function [rates, durations, containsBlanked, poissonCountMultipliers] = getSpikeMeanRate(td, unitName, varargin)
            p = inputParser();
            p.addParameter('invalidIfBlanked', false, @islogical); % if true, any trial that is partially blanked will be NaN, if false, the blanked region will be ignored and will not contribute to the time window used as the denominator for the rate calculation
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            counts = td.getSpikeCounts(unitName, 'combine', p.Results.combine);
            [durations, containsBlanked] = td.getValidDurationsForSpikeChannel(unitName, 'combine', p.Results.combine);

            if p.Results.invalidIfBlanked
                durations(containsBlanked) = NaN;
            end
            rates = counts ./ durations * td.timeUnitsPerSecond;
            poissonCountMultipliers = 1 ./ durations * td.timeUnitsPerSecond;
        end

        function [tf, blankingRegionsField] = hasBlankingRegions(td, unitName)
            if ischar(unitName)
                [tf, blankingRegionsField] = inner(unitName);
            else
                [tf, blankingRegionsField] = cellfun(@inner, unitName, 'UniformOutput', false);
                tf = cell2mat(tf);
            end

            function [tf, field] = inner(name)
                if ~td.hasSpikeChannel(name)
                    tf = false;
                    return;
                end
                tf = td.channelDescriptorsByName.(name).hasBlankingRegions;
                field = td.channelDescriptorsByName.(name).blankingRegionsField;
            end
        end

        function tf = hasSpikeWaveforms(td, unitNames)
            unitNames = string(unitNames);
            cds = td.getChannelDescriptor(unitNames);
            tf = arrayfun(@(cd) ~isempty(cd.waveformsField), cds);
        end

        function td = dropSpikeWaveforms(td, unitNames)
            td.warnIfNoArgOut(nargout);

            if nargin < 2
                unitNames = td.listSpikeChannels();
            end

            if ischar(unitNames)
                unitNames = {unitNames};
            end

            wavefields = cellvec(numel(unitNames));
            mask = falsevec(numel(unitNames));
            for iU = 1:numel(unitNames)
                unitName = unitNames{iU};
                if td.hasSpikeChannel(unitName) || td.hasSpikeArrayChannel(unitName)
                    wavefields{iU} = td.channelDescriptorsByName.(unitName).waveformsField;
                    if isempty(wavefields{iU}), continue; end
                    mask(iU) = true;
                    td.channelDescriptorsByName.(unitName) = td.channelDescriptorsByName.(unitName).removeWaveformsField();
                end
            end

            td = td.dropChannelFields(wavefields(mask));
        end

        function td = dropSpikeChannels(td)
            td.warnIfNoArgOut(nargout);
            td = td.dropChannels(td.listSpikeChannels('includeArraySubChannels', false));
            td = td.dropChannels(td.listExplicitSpikeArrays());
        end

        function td = maskSpikeChannelSpikesRaw(td, unitName, mask, varargin)
            td.warnIfNoArgOut(nargout);
            assert(ismatrix(mask) && size(mask, 1) == td.nTrials);

            %td.assertHasSpikeChannel(unitName);
            cd = td.channelDescriptorsByName.(unitName);
            assert(size(mask, 2) == cd.nChannels);

            times = cat(1, td.data.(unitName));
            notEmpty = ~cellfun(@isempty, times);
            times(notEmpty) = cellfun(@(t, m) t(m), times(notEmpty), mask(notEmpty), 'UniformOutput', false);

            hasWaves = cd.hasWaveforms;
            if hasWaves
                waves = cat(1, td.data.(cd.waveformsField));
                waves(notEmpty) = cellfun(@(w, m) w(m, :), waves(notEmpty), mask(notEmpty), 'UniformOutput', false);
            else
                waves = {};
            end

            td = td.setSpikeChannel(unitName, times, 'isAligned', false, 'waveforms', waves, varargin{:});
        end

        function td = trimSpikeChannelRaw(td, unitNames, varargin)
            % delete time points outside of a certain start stop interval
            % defaults to TrialStart:TrialStop
            p = inputParser();
            p.addOptional('startTimes', {}, @isvector);
            p.addOptional('stopTimes', {}, @isvector);
            p.addParameter('clearInvalidTrials', true, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            if ischar(unitNames)
                unitNames = {unitNames};
            end

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
            if p.Results.clearInvalidTrials
                startTimes(~td.valid) = NaN;
                stopTimes(~td.valid) = NaN;
            end

            tf = td.hasSpikeChannel(unitNames, 'includeArraySubChannels', false) | td.hasExplicitSpikeArray(unitNames);
            if any(~tf)
                error('Spike channels must be explicit spike arrays or separate explicit channels: %s', strjoin(unitNames(~tf)));
            end

            prog = ProgressBar(numel(unitNames), 'Trimming spike channels');
            for iU = 1:numel(unitNames)
                unitName = unitNames{iU};
                prog.update(iU, 'Trimming spike channel %s', unitName);

                cd = td.channelDescriptorsByName.(unitName);
                spikeField = cd.dataFieldPrimary;
                if isa(cd, 'SpikeArrayChannelDescriptor')
                    nUnits = cd.nChannels;
                else
                    nUnits = 1;
                end

                % trying this for spike arrays, will fail if empty on some
                % trials, but this shouldn't be the cease
                times = cat(1, td.data.(spikeField));
                timesMask = cell(td.nTrials, nUnits);

                for iT = 1:td.nTrials
                    for iV = 1:nUnits
                        timesMask{iT, iV} = times{iT, iV} >= startTimes(iT) & times{iT, iV} <= stopTimes(iT);
                    end
                end

                maskNeedsModification = ~cellfun(@all, timesMask);
                if ~any(maskNeedsModification)
                    continue;
                end

                td = td.maskSpikeChannelSpikesRaw(unitName, timesMask);
            end
            prog.finish();
        end

        function td = trimSpikeChannelToTrialStartEnd(td, unitNames)
            % Timepoints that lie outside of TrialStart and TrialStop will
            % never be accessible via getTimes since they will be filtered
            % out by the AlignInfo

            td.warnIfNoArgOut(nargout);
            % default is TrialStart and TrialEnd, so just pass it along
            td = td.trimSpikeChannelRaw(unitNames);
        end
        
        function tf = getSpikeChannelHasWaveforms(td, unitNames)
            unitNames = string(unitNames);

            tf = false(numel(unitNames), 1);
            for iU = 1:numel(unitNames)
                tf(iU) = strlength(td.getSpikeChannelMultiAccessInfo(unitNames(iU), 'waveforms')) > 0;
            end
        end

        function [wavesCell, waveTvec, timesCell, whichUnitCell] = getRawSpikeWaveforms(td, unitNames, varargin)
            % wavesCell is nTrials x nUnits of nSpikes x nWaveTime x nWaveChannels
            % timesCell is nTrials x nUnits
            % if combine is true, all spikes will be interleaved

            p = inputParser();
            p.addParameter('combine', false, @islogical);
            p.addParameter('applyScaling', true, @islogical);
            p.addParameter('fillMissingWithNaN', false, @islogical); % if any of the channels doesn't have waveforms, fill with nan
            p.parse(varargin{:});

            unitNames = string(unitNames);

            [fieldList, fieldIsArray, colIdxEachField, assignIdxEachField, nColumnsPerName, cdByField] = td.getSpikeChannelMultiAccessInfo(unitNames, 'waveforms');
            nUnits = sum(nColumnsPerName);
            nFields = numel(fieldList);
            wavesCellByUnit = cell(td.nTrials, nUnits);

            if nargout > 2 || p.Results.fillMissingWithNaN
                timesCell = td.getRawSpikeTimes(unitNames, 'combine', p.Results.combine);
            end
            for iF = 1:nFields
                fld = fieldList(iF);
                col = colIdxEachField{iF};
                idx = assignIdxEachField{iF};

                if string(fld) == ""
                    % no waveforms!
                    if p.Results.fillMissingWithNaN
                        nTime = numel(cdByField(iF).waveformsTime);
                        nWaveChan = cdByField(iF).waveformsNumChannels;
                        waveClass = cdByField(iF).waveformsOriginalDataClass;

                        % create a filler waveforms cell with nans of the
                        % right size
                        catData = cellfun(@(times) nan(numel(times), nTime, nWaveChan, waveClass), ...
                            timesCell(:, col), 'UniformOutput', false);
                    else
                        error('Unit %s does not have waveforms', cdByField(iF).name);
                    end
                    wavesCellByUnit(:, idx) = catData;

                else
                    % use actual waveforms
                    if fieldIsArray(iF)
                        catData = cat(1, td.data.(fld));
                    else
                        catData = {td.data.(fld)}';
                    end

                    wavesCellByUnit(:, idx) = catData(:, col);
                end
                if p.Results.applyScaling
                   wavesCellByUnit(:, idx) = cdByField(iF).getImpl().scaleWaveforms(wavesCellByUnit(:, idx));
                end
            end

            waveTvec = cdByField(1).waveformsTime;

            if p.Results.combine
                [wavesCell, whichUnitCell] = cellvec(td.nTrials);
                for iT = 1:td.nTrials
                    [wavesCell{iT}, whichUnitCell{iT}] = TensorUtils.catWhich(1, wavesCellByUnit{iT, :});
                end
            else
                wavesCell = wavesCellByUnit;
                if nargout > 3
                    whichUnitCell = cellfun(@(x) ones(size(x, 1), 1), wavesCell, 'UniformOutput', false);
                end
            end
        end

        function [wavesCell, waveTvec, timesCell] = getSpikeWaveforms(td, unitName, varargin)
            [wavesCell, waveTvec, timesCell] = td.getRawSpikeWaveforms(unitName, varargin{:});
            wavesCell = td.replaceInvalidMaskWithValue(wavesCell, []);
            timesCell = td.replaceInvalidMaskWithValue(timesCell, []);
        end

        function td = equalizeSpikeWaveformTimeVectors(td, unitNames)
            td.warnIfNoArgOut(nargout);
            if nargin < 2
                unitNames = td.listSpikeChannels();
            end
            waveTvecCell = cellfun(@(ch) makecol(td.channelDescriptorsByName.(ch).waveformsTime), ...
                unitNames, 'UniformOutput', false);

            delta = cellfun(@(x) nanmedian(diff(x)), waveTvecCell);
            if max(delta) - min(delta) > mean(delta) / 1000
                warning('Waveforms have different sampling rates, ignoring');
            end
            delta = median(delta);

            indZero = cellfun(@(x) TensorUtils.argMin(abs(x)), waveTvecCell);
            nSamples = cellfun(@numel, waveTvecCell);

            nSamplesPre = indZero-1;
            nSamplesPost = nSamples - indZero;

            N = max(nSamplesPre) + 1 + max(nSamplesPost);
            idxStart = max(nSamplesPre) - nSamplesPre + 1;
            idxStop = max(nSamplesPre) + 1 + nSamplesPost;
            nPadPre = max(nSamplesPre) - nSamplesPre;
            nPadPost = max(nSamplesPost) - nSamplesPost;

            waveTvec = (-max(nSamplesPre):max(nSamplesPost))' * delta;

            for iU = 1:numel(unitNames)
                if nPadPre(iU) == 0 && nPadPost(iU) == 0
                    continue;
                end
                cd = td.channelDescriptorsByName.(unitNames{iU});
                waveField = cd.waveformsField;
                waveData = {td.data.(waveField)}';

                newWaveData = cell(size(waveData));
                for iT = 1:size(newWaveData, 1)
                    if isfloat(waveData{iT})
                        newWaveData{iT} = nan(size(waveData{iT}, 1), N, 'like', waveData{iT});
                        newWaveData{iT}(:, idxStart(iU):idxStop(iU)) = waveData{iT};
                    else
                        newWaveData{iT} = padarray(...
                            padarray(waveData{iT}, [0 nPadPre(iU)], 'replicate', 'pre'), ...
                            [0 nPadPost(iU)], 'replicate', 'post');
                    end
                end

                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, waveField, newWaveData);
                cd.waveformsTime = waveTvec;
                td = td.setChannelDescriptor(unitNames{iU}, cd);
            end
        end
    end

    methods % Spike array channel methods
        function td = addSpikeArrayChannel(td, arrayName, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addOptional('electrodes', @(x) isnumeric(x) && isvector(x));
            p.addOptional('units', @(x) isnumeric(x) && isvector(x));
            p.addOptional('spikes', {}, @(x) ismatrix(x) && iscell(x));
            p.addParameter('isAligned', true, @isscalar);
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsTime', [], @(x) isempty(x) || isvector(x)); % common time vector to be shared for ALL waveforms for this channel
            p.addParameter('waveformsField', sprintf('%s_waveforms', arrayName), @ischar);
            p.addParameter('waveformsScaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsScaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is
            p.addParameter('sortQuality', NaN, @isscalar); % numeric scalar metric of sort quality
            p.addParameter('sortMethod', '', @ischar);
            p.addParameter('sortQualityEachTrial', [], @(x) isempty(x) || isvector(x));
            p.addParameter('blankingRegions', {}, @iscell); % nTrials x 1 cell of nIntervals x 2 matrices

            p.addParameter('raw', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            electrodes = p.Results.electrodes;
            units = p.Results.units;
            spikes = p.Results.spikes;

            % figure out nUnits and fill in missing units, electrodes, or
            % spikes
            specStr = 'At least one of units, electrodes, or spikes must be provided';
            if isempty(spikes)
                if isempty(electrodes)
                    assert(~isempty(units), specStr);
                    nUnits = numel(units);
                else
                    assert(~isempty(electrodes), specStr);
                    nUnits = numel(electrodes);
                end
            else
                assert(size(spikes, 1) == td.nTrials, 'Spikes must be nTrials x nUnits');
                nUnits = size(spikes, 2);
            end

            if isempty(electrodes)
                electrodes = nan(nUnits, 1);
            end
            if isempty(units)
                units = nan(nUnits, 1);
            end
            if isempty(spikes)
                spikes = cell(td.nTrials, nUnits);
            end

            cd = SpikeArrayChannelDescriptor.build(arrayName, electrodes, units);
            cd.sortQuality = p.Results.sortQuality;
            cd.sortMethod = p.Results.sortMethod;

            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial('raw', p.Results.raw);
            else
                % consider it aligned to trial start
                offsets = zerosvec(td.nTrials);
            end
            offsets = repmat(offsets, 1, nUnits);
            spikes = cellfun(@(x, y) makecol(x+y), spikes, num2cell(offsets), 'UniformOutput', false);

            % convert nTrials x nUnits cell array to nTrials x 1 of 1 x
            % nUnits fields
            splitCellIntoTrials = @(x) mat2cell(x, ones(td.nTrials, 1), nUnits);
            spikes = splitCellIntoTrials(spikes);
            channelData = {spikes};

            if ~isempty(p.Results.waveforms)
                waveforms = p.Results.waveforms;
                assert(isequal(size(waveforms), [td.nTrials, nUnits]), 'Size of Waveforms must be nTrials x nUnits');
                waveClass = TrialDataUtilities.Data.cellDetermineCommonClass(waveforms);
                cd = cd.addWaveformsField(p.Results.waveformsField, 'time', p.Results.waveformsTime, ...
                    'dataClass', waveClass, 'scaleFromLims', p.Results.waveformsScaleFromLims, ...
                    'scaleToLims', p.Results.waveformsScaleToLims);

                % undo scaling and convert back to memory scale
                if ~p.Results.waveformsInMemoryScale
                    impl = cd.getImpl();
                    waveforms = impl.unscaleWaveforms(waveforms);
                end

                channelData{end+1} = splitCellIntoTrials(waveforms);
            end

            if ~isempty(p.Results.sortQualityEachTrial)
                assert(numel(p.Results.sortQualityEachTrial) == td.nTrials);
                cd = cd.addSortQualityEachTrialField();
                channelData{end+1} = splitCellIntoTrials(p.Results.sortQualityEachTrial);
            end

            if ~isempty(p.Results.blankingRegions)
                blanking = p.Results.blankingRegions;
                assert(isequal(size(blanking), [td.nTrials nUnits]), 'Size of blankingRegions must be nTrials x nUnits');
                assert(all(cellfun(@(x) isempty(x) || (ismatrix(x) && size(x, 2) == 2), blanking)), ...
                    'Blanking regions must be nTrials cellvec with nRegions x 2 matrices');

                cd = cd.addBlankingRegionsField();

                % align the blanking regions
                blanking = cellfun(@(x,y) makecol(x+y), blanking, num2cell(offsets), 'UniformOutput', false);

                blanking = splitCellIntoTrials(blanking);

                % manually blank the spike times in the appropriate trials
                channelData{1} = TrialDataUtilities.SpikeData.blankRegionsEachTrial(channelData{1}, channelData{1}, blanking);

                channelData{end+1} = blanking;
            end

            td = td.addChannel(cd, channelData, 'updateValidOnly', ~p.Results.raw);
        end

        function td = incorporateSpikeChannelsIntoArrays(td, arrayNames)
            td.warnIfNoArgOut(nargout);

            if nargin < 2
                arrayNames = td.listImplicitSpikeArrays();
            end
            arrayNames = string(arrayNames);

            bad = td.hasExplicitSpikeArray(arrayNames);
            assert(~any(bad), 'Failed. Explicit spike array(s) %s already exists', strjoin(arrayNames(bad), ','));

            prog = ProgressBar(numel(arrayNames), 'Incorporating spike array channels');
            for iA = 1:numel(arrayNames)
                prog.update(iA, 'Incorporating spike array %s channels', arrayNames{iA});
                % find elec/unit combos in this array
                subNames = td.listSpikeChannels('array', arrayNames{iA});

                hasWaves = td.hasSpikeWaveforms(subNames);
                if any(hasWaves)
                    if ~all(hasWaves)
                        warning('Some but not all spike channels on array %s had waveforms; the rest will be nan-filled.', arrayNames{iA});
                    end
                    cd = td.getChannelDescriptor(subNames(1));

                    [waveforms, waveTvec] = td.getRawSpikeWaveforms(subNames, 'combine', false, 'fillMissingWithNaN', true, 'applyScaling', false);
                    waveArgs = {'waveforms', waveforms, 'waveformsTime', waveTvec, ...
                        'waveformsScaleFromLims', cd.waveformsScaleFromLims, 'waveformsScaleToLims', cd.waveformsScaleToLims, ...
                        'waveformsOriginalDataClass', cd.waveformsOriginalDataClass};
                else
                    waveArgs = {};
                end

                % get the times and add the array
                [~, electrodes, units] = arrayfun(@(subname) SpikeChannelDescriptor.parseArrayElectrodeUnit(subname), string(subNames));
                rawTimes = td.getRawSpikeTimes(subNames, 'combine', false);
                td = td.addSpikeArrayChannel(arrayNames{iA}, electrodes, units, rawTimes, waveArgs{:});

                % and drop the individual channels
                td = td.dropChannels(subNames);
            end
            prog.finish();
        end
    end

    methods % spike channel listing by arrays
        % there are two types of spike channels:
        %   direct: those with their own SpikeChannelDescriptor directly in channelDescriptors
        %   indirect: those that are part of a SpikeArrayChannelDescriptor
        % these functions below try to seamlessly blend the two so they are indistiguishable to the end user

        function arrays = listSpikeArrays(td)
            arrays = cat(1, td.listExplicitSpikeArrays(), td.listImplicitSpikeArrays());
        end

        function arrays = listImplicitSpikeArrays(td)
            channelDescriptors = td.getChannelDescriptorArray();
            if isempty(channelDescriptors)
                arrays = string([]);
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'SpikeChannelDescriptor'), channelDescriptors);
            arrays = arrayfun(@(cd) string(cd.array), channelDescriptors(mask));

            arrays = unique(arrays);
        end

        function names = listExplicitSpikeArrays(td)
            % internal use only, return SpikeArrayChannelDescriptors in
            % table
            channelDescriptors = td.getChannelDescriptorArray();
            if isempty(channelDescriptors)
                names = string([]);
                return;
            end
            mask = arrayfun(@(cd) isa(cd, 'SpikeArrayChannelDescriptor'), channelDescriptors);
            names = string({channelDescriptors(mask).name})';
        end

        function tf = hasExplicitSpikeArray(td, arrayName)
            tf = ismember(arrayName, td.listExplicitSpikeArrays());
        end

        function [names, channelDescriptors] = listSpikeChannels(td, varargin)
            p = inputParser();
            p.addParameter('includeDerivedChannelTypes', true, @islogical);

            p.addParameter('includeNonArraySubChannels', true, @islogical); % include those not in explicit arrays
            p.addParameter('includeArraySubChannels', true, @islogical); % include those in explicit arrays
            p.addParameter('array', '', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.addParameter('electrode', [], @(x) isvector(x));
            p.addParameter('unit', [], @isvector);
            p.addParameter('excludeUnit', [], @isvector); % for convenience, match unit and exclude these
            p.parse(varargin{:});

            % this should grab array sub channels if includeArraySubChannels
            [names, channelDescriptors] = td.listChannels('includeNamedSubChannels', p.Results.includeArraySubChannels, ...
                'channelDescriptorClass', ["SpikeChannelArrayDescriptor", "SpikeChannelDescriptor"], 'includeDerivedChannelTypes', true);
            mask = arrayfun(@(cd) isa(cd, 'SpikeChannelDescriptor'), channelDescriptors);
            names = names(mask);
            channelDescriptors = channelDescriptors(mask);

            mask = true(numel(names), 1);
            if ~p.Results.includeNonArraySubChannels
                mask = mask & arrayfun(@(cd) cd.isisColumnOfArray, channelDescriptors);
            end

            if nnz(mask) > 0
                if ~isempty(p.Results.array)
                    array = arrayfun(@(cd) string(cd.array), channelDescriptors);
                    mask = mask & ismember(array, string(p.Results.array));
                end
                if ~isempty(p.Results.electrode)
                    electrode = arrayfun(@(cd) cd.electrode, channelDescriptors);
                    mask = mask & ismember(electrode, p.Results.electrode);
                end
                if ~isempty(p.Results.unit)
                    unit = arrayfun(@(cd) cd.unit, channelDescriptors);
                    mask = mask & ismember(unit, p.Results.unit);
                end
                if ~isempty(p.Results.excludeUnit)
                    unit = arrayfun(@(cd) cd.unit, channelDescriptors);
                    mask = mask & ~ismember(unit, p.Results.excludeUnit);
                end
            end
            names = names(mask);
            channelDescriptors = channelDescriptors(mask);
        end

        function [subNames, parentNames, subChIdx, channelDescriptors]  = listSpikeNamedSubChannels(td)
            groups = td.listExplicitSpikeArrays();
            [subNames, parentNames, subChIdx, channelDescriptors] = td.listNamedSubChannels(groups);
        end

        function tf = hasSpikeChannel(td, name, varargin)
            tf = ismember(name, td.listSpikeChannels(varargin{:}));
        end

        function tf = hasSpikeArrayChannel(td, name, varargin)
            tf = ismember(name, td.listSpikeArrays(varargin{:}));
        end

        function td = remapSpikeArrayElectrodes(td, arrayChName, new_electrodes, varargin)
            p = inputParser();
            p.addParameter('resort', false, @islogical);
            p.parse(varargin{:});
            % works for spike arrays and continuous neural channel groups
            td.warnIfNoArgOut(nargout);

            cd = td.getChannelDescriptor(arrayChName);
            assert(isa(cd, 'SpikeArrayChannelDescriptor'));

            nElectrodes = numel(cd.electrodes);
            assert(nElectrodes == numel(new_electrodes));
            cd.electrodes = new_electrodes;
            td = td.setChannelDescriptor(arrayChName, cd);

            if p.Results.resort
                td = td.sortSpikeChannelArrayColumns(arrayChName);
            end
        end
        
        function td = mergeSpikeArraysEachUnit(td, arrayList, varargin)
            % TODO update this for spike array descriptors
            arrayList = string(arrayList);
            nA = numel(arrayList);
            assert(nA > 1);
            defaultAs = append("fuse_", strjoin(arrayList, '_'));
            td.warnIfNoArgOut(nargout);

            p = inputParser();
            p.addParameter('as', defaultAs, @TrialDataUtilities.String.isstringlike);
            p.addParameter('dropSourceArrays', true, @islogical);
            p.parse(varargin{:});

            % check that arrays have the same units and electrodes
            cds = td.getChannelDescriptorArray(arrayList);
            for a = 1:nA
                assert(isa(cds(a), 'SpikeArrayChannelDescriptor'), 'All arrays must be explicit SpikeArrayChannelDescriptors');
            end
            nUnitsEach = arrayfun(@(cd) cd.nChannels, cds);
            assert(all(nUnitsEach == nUnitsEach(1)));
            nU = nUnitsEach(1);
            
            electrodes = cds(1).electrodes;
            units = cds(1).units;
            for a = 2:nA
                assert(isequaln(cds(a).electrodes, electrodes) && isequaln(cds(a).units, units), ...
                    'Array %d electrodes and units do not match', a);
            end
            
            % gather all spike times and waveforms and other fields
            hasWaveforms = cds(1).hasWaveforms;
            hasBlankingRegions = cds(1).hasBlankingRegions;
            timesCell = cell(nA, td.nTrials, nU);
            
            if hasWaveforms
                wavesCell = cell(nA, td.nTrials);
                waveTvecCell = cell(nA, 1);
                [idxZero, nTimeWave, deltaTimeWave] = nanvec(nA);
                
                for a = 1:nA
                    [wavesCell(a, :, :), waveTvecCell{a}, timesCell(a, :, :)] = td.getRawSpikeWaveforms(arrayList{a});
                    [~, idxZero(a)] = min(abs(waveTvecCell{a}));
                    nTimeWave(a) = numel(waveTvecCell{a});
                    deltaTimeWave(a) = mean(diff(waveTvecCell{a}));
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
                for a = 1:nA
                    padLeft = zeroInd - idxZero(a);
                    padRight = nRightZeroGlobal - (nTimeWave(a) - idxZero(a));
                    if padLeft > 0 || padRight > 0
                        debug('Padding waveforms for %s to match new combined waveforms time vector\n', arrayList{a});
                        for t = 1:size(wavesCell, 2)
                            nWaves = size(wavesCell{a, t}, 1);
                            wavesCell{a, t} = cat(2, zeros(nWaves, padLeft), wavesCell{c, t}, zeros(nWaves, padRight));
                        end
                    end
                end
                
            else
                % no waveforms
                waveTvecGlobal = [];
                for a = 1:nA
                    timesCell(a, :, :) = td.getRawSpikeTimes(arrayList{a});
                end
            end
            
            if hasBlankingRegions
                blankingRegionsByArray = cell(nA, 1);
                for a = 1:nA
                    blankingRegionsByArray{a} = td.getRawSpikeBlankingRegions(arrayList{a});
                end
                
                blankingRegions = TrialDataUtilities.SpikeData.removeOverlappingIntervals(blankingRegionsByArray{:});
            else
                blankingRegions = {};
            end
            
            % combine the spiking data across channels
            times = cellvec(td.nTrials);
            if hasWaveforms
                waves = cellvec(td.nTrials, 1);
            else
                waves = {};
            end
            prog = ProgressBar(td.nTrials, 'Combining spike data');
            for iT = 1:td.nTrials
                prog.update(iT);
                for iU = 1:nU
                    [times{iT, iU}, sort_idx] = sort(cat(1, timesCell{:, iT, iU}));
                    if hasWaveforms
                        temp = cat(1, wavesCell{:, iT, iU});
                        waves{iT, iU} = temp(sort_idx, :);
                    end
                end
                
            end
            prog.finish();

            % add the new channel with combined data
            if p.Results.dropSourceArrays
                td = td.dropChannels(arrayList);
            end
            td = td.addSpikeArrayChannel(p.Results.as, electrodes, units, times, 'waveforms', waves, ...
                'waveformsTime', waveTvecGlobal, ...
                'isAligned', false, ...
                'blankingRegions', blankingRegions);
        end

        function td = sortSpikeChannelArrayColumns(td, array)
            td.warnIfNoArgOut(nargout);
            cd = td.getChannelDescriptor(array);
            assert(isa(cd, 'SpikeArrayChannelDescriptor'));

            [sorted, sortidx] = sortrows([cd.electrodes, cd.units]);
            electrodes = sorted(:, 1);
            units = sorted(:, 2);
            nUnits = numel(units);

            prog = ProgressBar(td.nTrials, 'Sorting spike array columns for %s', array);
            for iT = 1:td.nTrials
                prog.update(iT);
                if ~isempty(td.data(iT).(array))
                    td.data(iT).(array) = td.data(iT).(array)(:, sortidx);
                else
                    td.data(iT).(array) = cell(0, nUnits);
                end
            end
            prog.finish();

            cd.electrodes = electrodes;
            cd.units = units;
            td = td.setChannelDescriptor(array, cd);
            td = td.postDataChange(array);
        end
        
        function td = selectSpikeArrayChannels(td, array, mask, varargin)
            p = inputParser();
            p.addParameter('as', array, @TrialDataUtilities.String.isstringlike)
            p.parse(varargin{:});
            as = string(p.Results.as);
            
            td = td.reset();
            td.warnIfNoArgOut(nargout);
            cd = td.getChannelDescriptor(array);
            assert(isa(cd, 'SpikeArrayChannelDescriptor'));

            electrodes = cd.electrodes(mask);
            units = cd.units(mask);
            
            if cd.hasWaveforms
                [waves, waveTvecGlobal, times] = td.getRawSpikeWaveforms(array);
                waves = waves(:, mask);
            else
                waves = {};
                waveTvecGlobal = [];
                times = td.getRawSpikeTimes(array);
            end
            times = times(:, mask);
            
            if cd.hasBlankingRegions
                blankingRegions = td.getRawSpikeBlankingRegions(array);
                blankingRegions = blankingRegions(:, mask);
            else
                blankingRegions = {};
            end
            if cd.hasSortQualityEachTrial
                error('Not implemented');
            end
            
            if strcmp(array, as)
                td = td.dropChannels(array);
            end
                
            td = td.addSpikeArrayChannel(as, electrodes, units, times, 'waveforms', waves, ...
                'waveformsTime', waveTvecGlobal, ...
                'isAligned', false, ...
                'blankingRegions', blankingRegions);
        end

        function [fieldList, fieldIsArray, colIdxEachField, assignIdxEachField, nColumnsPerName, cdsByField] = ...
                getSpikeChannelMultiAccessInfo(td, names, type, slice)
            % internal method for arranging efficient access to multiple
            % array or non array spike fields at once. also supports
            % event channels and analog channels (marked as the time samples)
            %
            % type is either times (spike times) or waveforms (waveforms
            % field)
            %
            % K is the number of separate fields that need to be accessed
            % assignIdx is K x 1 cell of vectors of indices into names
            %   indicating which requested channel is located according to
            %   the corresponding location in fieldList and colIdxEachField
            %
            % fieldList is K x 1 strings of these fields
            %
            % colIdxEachField is K x 1 cell of vectors of column indices

            if nargin < 3
                type = 'times';
            end
            if nargin < 4
                slice = [];
            end

            if isnumeric(names)
                names = td.lookupSpikeChannelByIndex(names);
            end
            if ischar(names) || iscellstr(names) || isstring(names)
                names = string(names);
            end
            [fieldsByName, fieldIsArrayByName, colByName] = deal(cell(numel(names), 1));
            nColumnsPerName = nan(numel(names), 1);

            switch type
                case 'times'
                    cdField = 'dataFieldPrimary';
                case 'waveforms'
                    cdField = 'waveformsField';
                otherwise
                    error('Unknown type %s', type);
            end
            
            if ~isempty(slice)
                if iscell(slice)
                    assert(numel(slice) == numel(names));
                elseif numel(names) == 1
                    slice = {slice};
                else
                    error('Must provide one slice cell element for each unitNames element');
                end
            end
            sliceErrMsg = 'Slice argument not supported unless accessing SpikeArrayChannelDescriptor';

            % first lookup by names
            for iN = 1:numel(names)
                name = names{iN};
                if isnumeric(name)
                    % global index
                    name = td.lookupSpikeChannelByIndex(name);
                end
                if isfield(td.channelDescriptorsByName, name)
                    cd = td.channelDescriptorsByName.(name);
                    if isa(cd, 'SpikeChannelDescriptor')
                        fieldsByName{iN} = {cd.(cdField)};
                        colByName{iN} = cd.primaryDataFieldColumnIndex;
                        fieldIsArrayByName{iN} = false;
                        nColumnsPerName(iN) = 1;
                        assert(isempty(slice), sliceErrMsg);
                        
                    elseif isa(cd, 'EventChannelDescriptor')
                        fieldsByName{iN} = {cd.(cdField)};
                        colByName{iN} = 1;
                        fieldIsArrayByName{iN} = false;
                        nColumnsPerName(iN) = 1;
                        assert(isempty(slice), sliceErrMsg);

                    elseif isa(cd, 'AnalogChannelDescriptor') || isa(cd, 'AnalogChannelGroupDescriptor')
                        % refers to analog channel, use timestamps as
                        % "spikes" to make a sample time raster
                        fieldsByName{iN} = cd.timeField;
                        colByName{iN} = 1;
                        fieldIsArrayByName{iN} = false;
                        nColumnsPerName(iN) = 1;
                        assert(isempty(slice), sliceErrMsg);

                    elseif isa(cd, 'SpikeArrayChannelDescriptor')
                        % refers to array directly, include all channels
                        fieldsByName{iN} = repmat({cd.(cdField)}, cd.nChannels, 1);
                        colByName{iN} = (1:cd.nChannels)';
                        fieldIsArrayByName{iN} = true(cd.nChannels, 1);
                        nColumnsPerName(iN) = cd.nChannels;
                        
                        if ~isempty(slice) && ~isempty(slice{iN})
                            colByName{iN} = colByName{iN}(slice{iN});
                            fieldsByName{iN} = fieldsByName{iN}(slice{iN});
                            nColumnsPerName(iN) = numel(colByName{iN});
                        end 
                    else
                        error('Channel type %s not supported', class(cd));
                    end
                    cdByName(iN, 1) = cd; %#ok<AGROW>
                    
                else
                    assert(isempty(slice), sliceErrMsg);
                    
                    if contains(name, '(')
                        % parse array(idx) channel name
                        [tf, array, chidx] = td.parseIndexedSpikeArrayChannelName(name);
                        if ~tf
                            error('Channel %s was not not found or parsed successfully', name);
                        end
                        if ~isfield(td.channelDescriptorsByName, array)
                            error('Channel %s not found', array);
                        end
                        cd = td.channelDescriptorsByName.(array);
                        if ~isa(cd, 'SpikeArrayChannelDescriptor')
                            error('Channel %s is not a spike array channel', array);
                        end
                    else
                        % parse as arrayElec_unit name
                        [array, electrode, unit] = SpikeChannelDescriptor.parseArrayElectrodeUnit(name);
                        if isempty(array)
                             error('Channel %s was not not found or parsed successfully', name);
                        end
                        if ~isfield(td.channelDescriptorsByName, array)
                            error('Channel %s not found', array);
                        end
                        cd = td.channelDescriptorsByName.(array);
                        if ~isa(cd, 'SpikeArrayChannelDescriptor')
                            error('Channel %s is not a spike array channel', array);
                        end
                        chidx = cd.findSubChannelByElectrodeUnit(electrode, unit);
                        if isnan(chidx)
                            error('Channel %s not found within array %s', name, array);
                        end
                    end
                    if ~isfield(td.channelDescriptorsByName, array)
                        error('Channel %s not found', array);
                    end

                    fieldsByName{iN} = repmat(string(cd.(cdField)), numel(chidx), 1);
                    colByName{iN} = chidx;
                    fieldIsArrayByName{iN} = true(numel(chidx), 1);
                    cdByName(iN, 1) = cd; %#ok<AGROW>
                    nColumnsPerName(iN) = numel(chidx);
                end
            end
            
            [catFields, catWhichName] = TensorUtils.catWhich(1, fieldsByName{:});
            catFields = string(catFields);
            catCols = cat(1, colByName{:});
            catFieldIsArray = cat(1, fieldIsArrayByName{:});
            assignIntoList = (1:numel(catFields))';

            [fieldList, exemplarEntry, whichUniqueField] = unique(catFields);
            [assignIdxEachField, colIdxEachField] = deal(cell(numel(fieldList), 1));
            fieldIsArray = false(numel(fieldList), 1);
            for iF = 1:numel(fieldList)
                assignIdxEachField{iF} = assignIntoList(whichUniqueField == iF);
                colIdxEachField{iF} = catCols(whichUniqueField == iF);
                fieldIsArray(iF) = catFieldIsArray(exemplarEntry(iF));
                cdsByField(iF, 1) = cdByName(catWhichName(exemplarEntry(iF))); %#ok<AGROW>
            end
        end

        % this was replaced by getSpikeChannelMultiAccessInfo
%         function [cd, array, timesField, timesFieldColIdx, cdArray] = getSpikeChannelInfo(td, name)
%             % internal method for helping access spike channels by name
%             if isfield(td.channelDescriptorsByName, name)
%                 cd = td.channelDescriptorsByName.(name);
%                 if ~isa(cd, 'Sp ikeChannelDescriptor')
%                     error('Channel %s is not a spike channel', name);
%                 end
%                 array = cd.array;
%                 cdArray = [];
%
%             else
%                 % check for 'array(idx)' name
%                 [tf, array, ch_idx] = td.parseIndexedSpikeArrayChannelName(name);
%                 if tf
%                     if ~isfield(td.channelDescriptorsByName, array)
%                         error('Spike channel %s and spike channel array %s not found', name, array);
%                     end
%                     cdArray = td.channelDescriptorsByName.(array);
%                     if ~isa(cdArray, 'SpikeArrayChannelDescriptor')
%                         error('Channel %s is not a spike array channel', name);
%                     end
%
%                     cd = cdArray.buildIndividualSubChannel(ch_idx);
%                 else
%                     [array, elec, unit] = SpikeChannelDescriptor.parseArrayElectrodeUnit(name);
%                     if ~isfield(td.channelDescriptorsByName, array)
%                         error('Spike channel %s and spike channel array %s not found', name, array);
%                     end
%                     cdArray = td.channelDescriptorsByName.(array);
%                     if ~isa(cdArray, 'SpikeArrayChannelDescriptor')
%                         error('Channel %s is not a spike array channel', name);
%                     end
%
%                     cd = cdArray.buildIndividualSubChannelByElectrodeUnit(elec, unit);
%                 end
%             end
%
%             timesField = cd.dataFieldPrimary;
%             timesFieldColIdx = cd.primaryDataFieldColumnIndex;
%         end
%
        function [tf, array, idx] = parseIndexedSpikeArrayChannelName(td, name) %#ok<INUSL>
            % if name is a normal spike channel, do nothing
            % if it is a name like spikeArray(5), make a sub channel
            % channelDescriptor for analogGroup referring to column 5

            tokens = regexp(name, '(?<ch>[\w_]+)\((?<idx>\d+)\)', 'names');
            if isempty(tokens)
                tf = false;
                array = '';
                idx = [];
            else
                tf = true;
                array = tokens.ch;
                idx = str2double(tokens.idx);
            end
        end

        function assertHasSpikeChannel(td, name)
            assert(td.hasSpikeChannel(name), 'No spike channel %s found', name);
        end

        function names = lookupSpikeChannelByIndex(td, idx)
            allUnits = td.listSpikeChannels();
            names = allUnits(idx);
        end

        function info = listArrayElectrodesWithSpikeChannelsAsTable(td)
            % info is table with fields array (char) and electrode
            % (numeric)
            [channels, cdList] = td.listSpikeChannels();

            array = cellfun(@(cd) cd.array, cdList, 'UniformOutput', false);
            electrode = cellfun(@(cd) cd.electrode, cdList);

            t = table(array, electrode);
            [info, ~, idx] = unique(t);

            chMatch = cellvec(height(info));
            for r = 1:height(info)
                chMatch{r} = channels(idx == r);
            end
            info.channelList = chMatch;
        end

        function info = listSpikeChannelsAsTable(td, varargin)
            % info is table
            [name, cdList] = td.listSpikeChannels(varargin{:});
            array = arrayfun(@(cd) string(cd.array), cdList);
            electrode = arrayfun(@(cd) cd.electrode, cdList);
            unit = arrayfun(@(cd) cd.unit, cdList);
            info = table(name, array, electrode, unit);
        end

        function [nameCell, array, electrode, unitCell] = listSpikeChannelsGroupedByArrayElectrode(td)
            tbl = td.listSpikeChannelsAsTable;
            tae = tbl(:, {'array', 'electrode'});
            [uniq, ~, which] = unique(tae);

            nUniq = height(uniq);
            array = uniq.array;
            electrode = uniq.electrode;

            [nameCell, unitCell] = cellvec(nUniq);
            for iU = 1:nUniq
                nameCell{iU} = tbl.name(which == iU);
                unitCell{iU} = tbl.unit(which == iU);
            end
        end

        function [names, cds] = listSpikeChannelsOnArray(td, arrayName, varargin)
            [names, cds] = td.listSpikeChannels('array', arrayName, varargin{:});
        end

        function [names, cds] = listSpikeChannelsOnArrayElectrode(td, arrayName, electrodeNum, varargin)
            [names, cds] = td.listSpikeChannels('array', arrayName, 'electrode', electrodeNum, varargin{:});
        end

        function [names, cds] = listSpikeChannelsOnSameArrayElectrodeAs(td, chName, varargin)
            % [names, units] = listSpikeChannelsOnSameArrayElectrodeAs(td, chName)
            % chName is spike channel or continuous neural channel name
            cd = td.getChannelDescriptor(chName);
            [names, cds] = td.listSpikeChannelsOnArrayElectrode(cd.array, cd.electrode, varargin{:});
        end

        function [names, cds] = listSpikeChannelsWithUnitNumber(td, unit)
            [names, cds] = td.listSpikeChannels('unit', unit);
        end
    end

    methods % Generic add data methods
        function td = postDataChange(td, fieldsAffected) %#ok<INUSD>
            % call this after any changes to td.data.(fieldsAffected)
            td.warnIfNoArgOut(nargout);
        end

        function td = updatePostChannelDataChange(td, chName)
            td.warnIfNoArgOut(nargout);
            fields = td.channelDescriptorsByName.(chName).dataFields;
            td = td.postDataChange(fields);
        end

        function offsets = getTimeOffsetsFromZeroEachTrial(td, varargin)
            % when adding new data to the trial, all times are stored relative
            % to the current time zero. This will be overridden in
            % TrialDataConditionAlign. This will be added automatically
            % to all new channel time data to match the offsets produced when
            % getting data
            p = inputParser();
            p.addParameter('raw', false, @islogical);
            p.parse(varargin{:})

            offsets = zerosvec(td.nTrials);
        end

        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrialRaw(td)
            tMinByTrial = td.getEventFirst('TrialStart');
            tMaxByTrial = td.getEventFirst('TrialEnd');
        end

        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrial(td)
            % overriden by TDCA to be zero-relative

            [tMinByTrial, tMaxByTrial] = td.getTimeStartStopEachTrialRaw();
            tMinByTrial(~td.valid) = NaN;
            tMaxByTrial(~td.valid) = NaN;
        end

        function tf = alignIncludesFullTrial(td) %#ok<MANU>
            % just a standin b/c this is redefined in TDCA
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
                %has = structfun(@(cd) td.hasChannel(cd.name), cds);
                has = structfun(@(cd) isfield(td.channelDescriptorsByName, cd.name), cds);
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
            p.addOptional('valueCell', {}, @(x) ~ischar(x));
            p.addParameter('clearIfPresent', false, @islogical);
            % p.addParameter('ignoreOverwriteChannel', false, @islogical);
            p.addParameter('updateValidOnly', true, @islogical);

            p.addParameter('ignoreDataFields', false, @islogical); % if true, neither clear nor set channel data fields, used by addAnalogChannelGroup initially
            p.addParameter('ignoreExisting', false, @islogical)
            p.parse(varargin{:});
            valueCell = makecol(p.Results.valueCell);

            td.warnIfNoArgOut(nargout);

            % check for overwrite if requested
            assert(isa(cd, 'ChannelDescriptor'), 'Argument cd must be ChannelDescriptor');

            % give the channelDescriptor a chance to initialize itself
            cd = cd.initialize();

            alreadyHasChannel = td.hasChannel(cd.name);
            if alreadyHasChannel && ~p.Results.ignoreExisting
                warning('Overwriting existing channel with name %s', cd.name);
                td = td.dropChannels(cd.name);

                %                 % check that existing channel matches
                %                 assert(isequaln(cd, td.channelDescriptorsByName.(cd.name)), ...
                %                     'ChannelDescriptor for channel %s does not match existing channel', cd.name);
            end

            td.channelDescriptorsByName.(cd.name) = cd;

            if ~p.Results.ignoreDataFields
                % touch each of the value fields to make sure they exist
                for iF = 1:cd.nFields
                    missing = cd.missingValueByField{iF};
                    if ~isfield(td.data, cd.dataFields{iF}) || ~cd.isShareableByField(iF)
                        % clear if it's missing, or if its there but not
                        % shared, since we're overwriting it
                        td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, cd.dataFields{iF}, missing, [], true); % true ensures missing is treated as scalar if it is cell
                    end
                end

                if isempty(valueCell) && (~alreadyHasChannel || p.Results.clearIfPresent)
                    td = td.clearChannelData(cd.name, 'fieldMask', ~cd.isShareableByField, ...
                        'updateOnlyValidTrials', false);
                elseif ~isempty(valueCell)
                    % clear on fields where no values provided and it's not shared,
                    % set on fields where values are provided
                    nonEmptyMask = ~cellfun(@isempty, valueCell);
                    td = td.clearChannelData(cd.name, 'fieldMask', ~nonEmptyMask & ~cd.isShareableByField, ...
                        'updateOnlyValidTrials', false);
                    td = td.setChannelData(cd.name, valueCell, 'fieldMask', nonEmptyMask, ...
                        'updateValidOnly', p.Results.updateValidOnly);
                end
            end
        end

        function td = copyChannel(td, oldName, newName, varargin)
            p = inputParser();
            p.addParameter('clearForInvalid', false, @islogical);
            p.parse(varargin{:});

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

            if p.Results.clearForInvalid
                td = td.clearChannelData(newName, 'updateMask', ~td.valid, 'updateOnlyValidTrials', false);
            end

            td = td.updatePostChannelDataChange(newName);
        end

        function tdDest = copyChannelToSecondInstance(td, tdDest, oldName, varargin)
            p = inputParser();
            p.addParameter('as', oldName, @TrialDataUtilities.String.isstringlike);
            p.addParameter('clearForInvalid', false, @islogical);
            p.parse(varargin{:});

            % copy channel oldName to newName
            td.warnIfNoArgOut(nargout);
            newName = p.Results.as;

            assert(td.nTrials == tdDest.nTrials, 'Trial counts of td and tdDest differ');

            alreadyHasChannel = tdDest.hasChannel(newName);
            if alreadyHasChannel
                warning('Overwriting existing channel with name %s', newName);
                tdDest = tdDest.dropChannels(newName);
            end

            % ask the channel descriptor to rename itself and store the copy
            [cd, dataFieldRenameMap] = td.channelDescriptorsByName.(oldName).rename(newName);
            tdDest.channelDescriptorsByName.(newName) = cd;
            flds = fieldnames(dataFieldRenameMap);

            % check that no other channels in tdDest (besides newName) were using that field
            for iF = 1:numel(flds)
                oldField = flds{iF};
                newField = dataFieldRenameMap.(oldField);
                otherCh = setdiff(tdDest.getChannelsReferencingFields(newField), newName);
                if ~isempty(otherCh)
                    error('Channel %s in destination td also references field %s', ...
                        TrialDataUtilities.String.strjoin(otherCh), newField);
                end
            end

            % then copy new channel fields
            for iF = 1:numel(flds)
                oldField = flds{iF};
                newField = dataFieldRenameMap.(oldField);
                tdDest.data = copyStructField(td.data, tdDest.data, oldField, newField);
            end

            % if analog channel in group, separate it from the group
            % this will copy the data, but it won't get rid of the old data
            % since the original channel still references that column
            if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix
                tdDest = tdDest.separateAnalogChannelFromGroup(newName);
            end

            if p.Results.clearForInvalid
                % we use td.valid since we want td source to control which data points are considered valid, not tdDest
                tdDest = tdDest.clearChannelData(newName, 'updateMask', ~td.valid, 'updateOnlyValidTrials', false);
            end

            tdDest = tdDest.updatePostChannelDataChange(newName);
        end

        function td = renameChannel(td, oldName, newName)
            % rename channel name to newName
            % if channel name's primary data field is also name, rename
            % that field too to newName
            td.warnIfNoArgOut(nargout);

            cd = td.channelDescriptorsByName.(oldName);

            % update channel descriptor directly
            [cd, dataFieldRenameMap] = cd.rename(newName);

            % then rename the actual channel
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

        function td = renameAnalogTimeField(td, name, timeField)
            td.warnIfNoArgOut(nargout);
            cd = td.channelDescriptorsByName.(name);

            if nargin < 3
                timeField = cd.suggestFieldName('time');
            end
            oldTimeField = cd.timeField;

            if strcmp(oldTimeField, timeField)
                return;
            end

            td.channelDescriptorsByName.(name) = cd.renameDataField('time', timeField);
            td = td.renameDataField(oldTimeField, timeField, name);
        end

        function td = setChannelDescriptor(td, name, cd)
            td.warnIfNoArgOut(nargout);
            assert(isa(cd, 'ChannelDescriptor'));
            assert(isfield(td.channelDescriptorsByName, name));
            td.channelDescriptorsByName.(name) = cd;
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

            conflictChannels = setdiff(td.getChannelsReferencingFields(newFieldName), ignoreChannelList);
            if ~isempty(conflictChannels)
                error('Channels %s are referencing field name %s', TrialDataUtilities.String.strjoin(conflictChannels));
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

        function td = dropUnusedDataFields(td)
            td.warnIfNoArgOut(nargout);

            flds = fieldnames(td.data);
            maskRemove = falsevec(numel(flds));
            for iF = 1:numel(flds)
                maskRemove(iF) = isempty(td.getChannelsReferencingFields(flds{iF}));
            end

            if any(maskRemove)
                td.data = rmfield(td.data, flds(maskRemove));
            end

            td.data = orderfields(td.data);
        end

        function td = copyRenameSharedChannelFields(td, names, fieldMask)
            % copy all fields which belong to this channel which are shared
            % with other channes, rename those fields, and mark the
            % updated fields inside that fields channelDescriptor. if name
            % is a cellstr of channel names, only consider fields shared
            % with channels not in that set of channels.
            td.warnIfNoArgOut(nargout);

            if nargin < 3
                fieldMask = [];
            end
            names = string(names);
            fieldsRenamed = struct();

            for iC = 1:numel(names)
                fixConflictsWithChannel(names(iC), names);
            end

%             [groupNames, singleNames] = sortOutGroups(names);
%             for iC = 1:numel(groupNames)
%                 fixConflictsWithChannel(groupNames{iC}, union(groupNames, names));
%             end
%             for iC = 1:numel(singleNames)
%                 fixConflictsWithChannel(singleNames{iC}, union(groupNames, names));
%             end
%
%             function [groupNames, singleNames] = sortOutGroups(names)
%                 % check for analog channel groups: if one channels is included, all
%                 % must be included
%                 inGroup = falsevec(numel(names));
%                 groupNames = cellvec(numel(names));
%                 for iN = 1:numel(names)
%                     if td.hasAnalogChannelGroup(names{iN})
%                         groupNames{iN} = names{iN};
%                         inGroup(iN) = true;
%                     elseif td.hasAnalogChannel(names{iN})
%                         [inGroup(iN), groupNames{iN}] = td.isAnalogChannelInGroup(names{iN});
%                     end
%                 end
%
%                 if any(inGroup)
%                     groupNames = unique(groupNames(inGroup));
%                     for iG = 1:numel(groupNames)
%                         other = td.listAnalogChannelsInGroup(groupNames{iG});
%                         if ~all(ismember(other, names))
%                             error('If renaming channel fields referenced by channels in analog channel group (%s) all channels in that group must be referenced', groupNames{iG});
%                         end
%                     end
%                 else
%                     groupNames = {};
%                 end
%
%                 singleNames = names(~inGroup);
%             end

            function fixConflictsWithChannel(chName, exclude)
                % resolve conflicts shared between channel chName and any
                % channels not in exclude. Either by renaming fields in
                % chName, or by renaming other channels to avoid it.

                cd = td.getChannelDescriptor(chName);

                if isempty(fieldMask)
                    fieldMaskThis = true(cd.nFields, 1);
                else
                    fieldMaskThis = TensorUtils.vectorIndicesToMask(fieldMask, cd.nFields);
                end

                for iF = 1:cd.nFields
                    if ~cd.isShareableByField(iF) || ~fieldMaskThis(iF), continue; end
                    if ismember(cd.dataFields{iF}, fieldnames(fieldsRenamed))
                        % already renamed earlier, change cd to reflect that
                        doFieldRename(chName, iF);
                    end

                    [fullList, whichField] = td.getChannelsReferencingFields(cd.dataFields{iF});
                    [otherChannels, mask] = setdiff(fullList, exclude);
                    if isempty(otherChannels), continue; end
                    whichField = whichField(mask);

                    % check whether the field matches my suggested name for
                    % it, which suggests it was mine and the other channels
                    % are sharing it.
                    template = cd.suggestFieldName(iF);
                    myFieldName = strcmp(template, cd.dataFields{iF});

                    if myFieldName
                        % have all the other channels rename their fields
                        % to something else. This will handle group sub
                        % channels intelligently
                        for iO = 1:numel(otherChannels)
                            doFieldRename(otherChannels{iO}, whichField(iO));
                        end
                    else
                        newName = doFieldRename(chName, iF);
                        cd = td.getChannelDescriptor(chName);

                        % rename all sub channels if any
                        if isa(cd,'AnalogChannelGroupDescriptor')
                            groupSubChList = td.listAnalogChannelsInGroup(chName);
                            for iSub = 1:numel(groupSubChList)
                                td.channelDescriptorsByName.(groupSubChList{iSub}) = ...
                                    td.channelDescriptorsByName.(groupSubChList{iSub}).renameDataField(iF, newName);
                            end
                        end
                    end
                end
            end

            function newName = doFieldRename(chName, iF)
                % rename field iF of channel chName to a new field name

                cd = td.getChannelDescriptor(chName);
                if ~cd.isShareableByField(iF)
                    return
                end
                if ismember(cd.dataFields{iF}, fieldnames(fieldsRenamed))
                    % already renamed earlier, change cd to reflect that
                    newName = fieldsRenamed.(cd.dataFields{iF});
                else
                    % let the group suggest the name if in group
                    cdSuggest = cd;
                    if isa(cd, 'AnalogChannelGroup')
                        % in group
                        [tf, group] = td.isAnalogChannelInGroup(cd);
                        if tf
                            cdSuggest = td.getChannelDescriptor(group);
                        end
                    end

                    % not yet renamed, come up with new field name and copy
                    % data over
                    template = cdSuggest.suggestFieldName(iF);
                    template = template(1:min(numel(template), namelengthmax()-10));
                    newName = matlab.lang.makeUniqueStrings(template, fieldnames(td.data));
                    td.data = copyStructField(td.data, td.data, cd.dataFields{iF}, newName);
                    fieldsRenamed.(cd.dataFields{iF}) = newName;
                end

                cd = cd.renameDataField(iF, newName);
                td.channelDescriptorsByName.(cd.name) = cd;
            end
        end

        function td = clearChannelData(td, name, varargin)
            p = inputParser();
            p.addParameter('fieldMask', [], @islogical);
            p.addParameter('updateMask', truevec(td.nTrials), @islogical);
            p.addParameter('updateOnlyValidTrials', false, @islogical);
            p.parse(varargin{:});

            cd = td.channelDescriptorsByName.(name);
            fieldMask = p.Results.fieldMask;
            if isempty(fieldMask)
                fieldMask = true(cd.nFields, 1);
            end
            if ~any(fieldMask)
                return;
            end

            updateMask = p.Results.updateMask;
            updateMask = TensorUtils.vectorIndicesToMask(updateMask, td.nTrials);
            if p.Results.updateOnlyValidTrials
                updateMask = updateMask & td.valid;
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
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, dataFields{iF}, val, updateMask);
                %                 for iT = 1:td.nTrials
                %                     td.data(iT).() = val;
                %                 end
            end

        end

        function td = setChannelDisplayGroup(td, names, displayGroup)
            td.warnIfNoArgOut(nargout);

            names = string(names);

            for iN = 1:numel(names)
                name = names(iN);
                cd = td.channelDescriptorsByName.(name);
                cd.displayGroup = displayGroup;
                td.channelDescriptorsByName.(name) = cd;
            end
        end

        function [td, trialsUpdated, fieldsUpdated] = setChannelData(td, name, valueCell, varargin)
            % this will replace each trial's entire data, rather than just the aligned window
            % note that by default, updateValidOnly is true, meaning that
            % the values on invalid trials will not be updated
            % trialsUpdated is a mask indicating which trials were changed
            % (or cleared)
            % fieldsUpdated is a mask indicating which of the data fields
            % were updated as well
            %
            % data passed into setChannelData should already be in memory data type and scale,
            % so the caller should handle any access to memory conversion
            %
            p = inputParser();
            p.addParameter('fieldMask', [], @islogical);
            p.addParameter('clearForInvalid', false, @islogical);
            p.addParameter('updateValidOnly', true, @islogical);
            p.addParameter('updateMask', [], @isvector);
            p.addParameter('deferPostDataChange', false, @islogical);
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
            td = td.clearChannelData(cd.name, 'fieldMask', fieldMissing & fieldMask, ...
                'updateOnlyValidTrials', false);

            % avoid overwriting data shared by other channels
            td = td.copyRenameSharedChannelFields(cd.name, fieldMask);
            cd = td.channelDescriptorsByName.(name);
            dataFields = cd.dataFields;
            impl = cd.getImpl();

            updateMask = updateMaskManual;
            if p.Results.updateValidOnly
                updateMask = updateMask & td.valid;
            end

            for iF = 1:nFields
                % only touch specified fields
                if ~fieldMask(iF), continue; end

                if size(valueCell{iF}, 2) > 1 && size(valueCell{iF}, 1) == td.nTrials
                    % split along trial axis into nTrials separate cells
                    valueCell{iF} = mat2cell(valueCell{iF}, ones(td.nTrials, 1), size(valueCell{iF}, 2));
                end


                % fields that have CELL values could be passed in as
                % nTrials x nFields, here we're choosing to have this
                % conversion be handled by the caller
                assert(numel(valueCell{iF}) == td.nTrials, 'Each field''s value must be nTrials x 1 cell');

                % validate the overall structure of the data (scalar vs.
                % numeric vs. vector, etc.) and update the memory class to
                % match the data passed in (e.g. uint16 --> double)
                [cd, valueCell{iF}] = impl.checkConvertDataAndUpdateMemoryClassToMakeCompatible(iF, valueCell{iF});
                impl = cd.getImpl();

                if p.Results.clearForInvalid
                    % here we want the update mask to stay the same as
                    % updateMaskManual so that everything gets updated with
                    % the cleared value
                    valueCell{iF} = td.replaceInvalidMaskWithValue(valueCell{iF}, cd.missingValueByField{iF});
                end

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
                    vals = valueCell{iF};
                    if ~iscell(vals)
                        vals = TensorUtils.splitAlongDimension(vals, 1);
                    end
                    td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, dataFields{iF}, TensorUtils.selectAlongDimension(vals, 1, updateMask), updateMask);

                    %                   for iT = 1:numel(td.data)
                    %                       if updateMask(iT)
                    %                            td.data(iT).(dataFields{iF}) = valueCell{iF}{iT};
                    %                       end
                    %                  end
                end
            end

            td.channelDescriptorsByName.(cd.name) = cd;

            if ~p.Results.deferPostDataChange
                td = td.postDataChange(dataFields(fieldMask));
            end

            fieldsUpdated = dataFields(fieldMask);
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
                if isempty(dataCell{i}), continue, end

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
            hold(axh, 'off');
        end
    end

    % save / load wrappers for CacheCustomSaveLoad
    properties
        cacheWithSaveFast = false;
        saveFastPartitionInfo = struct();
        saveFastPartitionWaveforms = false;
        saveFastTrialsPerChunk = 1;
    end

    % CacheCustomSaveLoad impl
    methods
        function tf = getUseCustomSaveLoad(td, info) %#ok<INUSD>
            tf = td.cacheWithSaveFast;
        end

        function token = saveCustomToLocation(td, location)
            td.saveFast(location, 'partitions', td.saveFastPartitionInfo, ...
                'partitionWaveforms', td.saveFastPartitionWaveforms);
            token = [];
        end
    end

    % CacheCustomSaveLoad load impl
    methods(Static)
        function data = loadCustomFromLocation(location, token, varargin)  %#ok<INUSL>
            data = TrialData.loadFast(location, varargin{:});
        end
    end
end
