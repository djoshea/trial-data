classdef TrialDataConditionAlign < TrialData

    % ConditionInfo and AlignInfo stores
    properties(SetAccess=protected)
        % ConditionInfo instance
        conditionInfo

        % AlignInfo instance
        alignInfoSet

        % for functions that operate on one align info, this index is the one which is active
        % use .useAlign(idx) to set
        alignInfoActiveIdx

        % time units elapsing between successive alignments, used mainly
        % for the purposes of plotting
        interAlignGaps
        
        % alternatively specify the offset where each align occurs (each align's zero will occur at this offset)
        manualAlignTimeOffsets
    end

    % properties which are stored in odc, see get/set below
    properties(Transient, Dependent, SetAccess=protected)
        % scalar struct with fields .event
        % each field containing an nTrials x nOccurrence table of event times
        eventData

        % scalar struct with fields .event
        % each field containing the number of occurrences of each event
        eventCounts

        % nAlign x 1 cell of AlignSummary instances
        alignSummarySet

        % size(conditions) by nRandomized
        randomizedListsByCondition

        conditionInfoRandomized
    end

    properties(SetAccess=protected)
        nRandomized = 10; % number of cached copies of the randomized condition lists to keep around
        % set with setNumRandomized
    end

    % Properties which read through to ConditionInfo
    properties(Dependent, SetAccess=protected)
        attributeParams % feeds thru to .conditionInfo.attributeNames

        isGrouped
        nConditions
        listByCondition
        listByConditionWeights
        nTrialsByCondition
        conditionIdx
        conditionSubs
        conditionMembership
        conditionNamesEachTrial
        conditions
        conditionsAxisAttributesOnly
        conditionsAsStrings
        conditionNames
        conditionNamesShort
        conditionNamesMultiline
        conditionNamesMultilineShort
        conditionAppearances
        conditionColors
        conditionsSize
        conditionsSizeNoExpand
        conditionIncludeMask

        nAxes
        axisValueLists
        axisValueListsAsStrings
        axisValueListsAsStringsShort
        axisNames
        axisDescriptions

        conditionAppearanceFn
        conditionNameFn

        hasRandomizationActive
        hasRandomizationSpecified
        randomizationDescription
        randomSeed % use setRandomSeed
    end

    % Simple dependent properties
    properties(Dependent, SetAccess=protected)
        nAlign

        alignInfoActive % alignInfo currently active

        alignSummaryActive % alignSummary currently active

        conditionDescriptor % conditionInfo reduced to conditionDescriptor
    end

    % Properties which read through to AlignInfo
    properties(Dependent, SetAccess=protected)
        minTimeDelta
    end

    % Initializing and building
    methods
        function td = TrialDataConditionAlign(varargin)
            td = td@TrialData(varargin{:});
            td = td.rebuildOnDemandCache();
            td = td.initializeConditionInfo();
            td = td.initializeAlignInfo();
            td = td.invalidateValid();
        end

        function td = rebuildOnDemandCache(td)
            td.warnIfNoArgOut(nargout);
            td.odc = TrialDataConditionAlignOnDemandCache();
        end
    end
    
    methods(Static)
        function td = loadFast(location, varargin)
            td = loadFast@TrialData(location, varargin{:});
            td = TrialDataConditionAlign(td);
        end
    end

    % Simple dependent getters
    methods
        function v = get.nAlign(td)
            v = numel(td.alignInfoSet);
        end

        function cd = get.conditionDescriptor(td)
            cd = ConditionDescriptor.fromConditionDescriptor(td.conditionInfo.fixAllValueLists);
        end
    end

    % get / set accessors that read / write through to ODC
    methods
        function v = get.eventData(td)
            v = td.odc.eventData;
            if isempty(v)
                td.buildEventData();
                v = td.odc.eventData;
            end
        end

        function td = set.eventData(td, v)
            td.odc = td.odc.copy();
            td.odc.eventData = v;
        end

        function v = get.eventCounts(td)
            v = td.odc.eventCounts;
            if isempty(v)
                td.buildEventData();
                v = td.odc.eventCounts;
            end
        end

        function td = set.eventCounts(td, v)
            td.odc = td.odc.copy();
            td.odc.eventCounts = v;
        end

        function v = get.alignSummarySet(td)
            v = td.odc.alignSummarySet;
            if isempty(v)
                td.buildAlignSummarySet();
                v = td.odc.alignSummarySet;
            end
        end

        function td = set.alignSummarySet(td, v)
            td.odc = td.odc.copy();
            td.odc.alignSummarySet = v;
        end

        function v = get.randomizedListsByCondition(td)
            v = td.odc.randomizedListsByCondition;
            if isempty(v)
                td.odc.randomizedListsByCondition = td.buildRandomizedListsByCondition();
                v = td.odc.randomizedListsByCondition;
            end
        end

        function td = set.randomizedListsByCondition(td, v)
            td.odc = td.odc.copy();
            td.odc.randomizedListsByCondition = v;
        end

        function v = get.conditionInfoRandomized(td)
            v = td.odc.conditionInfoRandomized;
            if isempty(v)
                % we just copy the condition Info on first access, will be
                % subsequently modified
                td.odc.conditionInfoRandomized = td.conditionInfo;
                v = td.odc.conditionInfoRandomized;
            end
        end

        function td = set.conditionInfoRandomized(td, v)
            td.odc = td.odc.copy();
            td.odc.conditionInfoRandomized = v;
        end

        function v = get.randomSeed(td)
            v = td.conditionInfoRandomized.randomSeed;
        end
    end

    % build* methods for properties stored in odc
    methods(Access=protected)
        function buildEventData(td)
            % build a cached version of event times into a matrix for easy alignment
            % sets .eventData and .eventCounts

            evStruct = td.getRawEventFlatStruct();
            evList = fieldnames(evStruct);
            nEvents = numel(evList);
            eventCounts = struct(); %#ok<*PROP>
            eventData = struct();

            for iE = 1:nEvents
                ev = evList{iE};
                times = evStruct.(ev);
                if iscell(times)
                    % event may happen zero, one, or multiple times
                    % convert to nTrials x nOccur matrix
                    counts = cellfun(@(x) nnz(~isnan(x)), times);
                    maxCount = max(counts);
                    timeMat = nan(td.nTrials, maxCount);
                    for iT = 1:td.nTrials
                        timeMat(iT, 1:counts(iT)) = times{iT};
                    end
                    eventCounts.(ev) = counts;
                    eventData.(ev) = timeMat;
                else
                    eventCounts.(ev) = makecol(double(~isnan(times)));
                    eventData.(ev) = makecol(times);
                end
            end

            c = td.odc;
            c.eventCounts = eventCounts;
            c.eventData = eventData;
        end

        function td = updateEventData(td, eventsUpdate)
            % similar to buildEventData, but updates a specific subset of
            % fields to save time when event data is changed

            % since this isn't computed on demand, a copy is needed
            c = td.odc.copy();

            if isstringlike(eventsUpdate)
                eventsUpdate = {eventsUpdate};
            end

            if isempty(td.odc.eventCounts) || isempty(td.odc.eventData)
                % hasn't been computed yet, no need to update
                return;
            else
                evPartialStruct = td.getRawEventFlatStruct(eventsUpdate);
                nEvents = numel(eventsUpdate);
                for iE = 1:nEvents
                    ev = eventsUpdate{iE};
                    times = evPartialStruct.(ev);
                    if iscell(times)
                        % event may happen zero, one, or multiple times
                        % convert to nTrials x nOccur matrix
                        counts = cellfun(@(x) nnz(~isnan(x)), times);
                        maxCount = max(counts);
                        timeMat = nan(td.nTrials, maxCount);
                        for iT = 1:td.nTrials
                            timeMat(iT, 1:counts(iT)) = times{iT};
                        end
                        c.eventCounts.(ev) = counts;
                        c.eventData.(ev) = timeMat;
                    else
                        c.eventCounts.(ev) = makecol(double(~isnan(times)));
                        c.eventData.(ev) = makecol(times);
                    end
                end
            end
            td.odc = c;

            % update the cached events inside the align info as well
            td = td.applyAlignInfoSet();
        end

        function buildAlignSummarySet(td)
            alignSummarySet = cell(td.nAlign, 1);
            for i = 1:td.nAlign
                alignSummarySet{i} = AlignSummary.buildFromConditionAlignInfo(td.conditionInfo, td.alignInfoSet{i});
                alignSummarySet{i}.timeUnitName = td.timeUnitName;
            end

            c = td.odc;
            c.alignSummarySet = alignSummarySet;
        end

        function alignSummarySet = buildAlignSummarySetWithTrials(td, trialIdx)
            td = td.withTrials(trialIdx);
            alignSummarySet = cell(td.nAlign, 1);
            for i = 1:td.nAlign
                alignSummarySet{i} = AlignSummary.buildFromConditionAlignInfo(td.conditionInfo, td.alignInfoSet{i});
                alignSummarySet{i}.timeUnitName = td.timeUnitName;
            end
        end

        function buildValid(td)
            % builds .valid and .invalidCause, no odc copy
            % be very careful not to request .valid or .invalidCause here
            % as this will lead to infinite recursion.
            %
            % valid is the intersection of manualValid, conditionInfo valid,
            % and all AlignInfo valid
            [tempValid, tempCause] = td.getTemporaryValid();
            [manualValid, manualCause] = td.getManualValid();
            cvalid = td.conditionInfo.computedValid;
            avalid = truevec(td.nTrials);
            for iA = 1:td.nAlign
                avalid = avalid & td.alignInfoSet{iA}.computedValid;
            end

            valid = manualValid & tempValid & cvalid & avalid;

            cause = manualCause;
            explained = ~manualValid;

            % invalidated temporarily
            tmask = ~tempValid & ~explained;
            cause(tmask) = cellfun(@(s) ['(temporary) ' s], ...
                tempCause(tmask), 'UniformOutput', false);

            % invalid by condition info
            cmask = ~td.conditionInfo.computedValid & ~explained;
            cause(cmask) = cellfun(@(s) ['(temporary) ConditionInfo: ' s], ...
                td.conditionInfo.invalidCause(cmask), 'UniformOutput', false);

            % invalid by each align info
            for iA = 1:td.nAlign
                amask = ~td.alignInfoSet{iA}.computedValid & ~explained;
                cause(amask) = cellfun(@(s) sprintf('(temporary) AlignInfo %d: %s', ...
                    iA, s), td.alignInfoSet{iA}.invalidCause(amask), 'UniformOutput', false);
            end

            cause(valid) = {''};

            % don't copy in build methods
            c = td.odc;
            c.valid = valid;
            c.invalidCause = string(cause);
            td.odc = c;
        end
    end

    % General utilites
    methods
        function td = postAddNewTrials(td)
            td.warnIfNoArgOut(nargout);

            % copy and flush valid
            td.odc = td.odc.copy();
            td.odc.flush();

            td = td.reinitializeConditionAlign();

            td = postAddNewTrials@TrialData(td); % will call update valid
        end

        function td = postDataChange(td, fieldsAffected)
            td.warnIfNoArgOut(nargout);

            if isempty(fieldsAffected)
                return;
            end
            needUpdate = false;

            % check if any event fields were affected. if so, all align
            % info's need to be reapplied since they cache all events, not
            % just the ones they need
            eventFields = td.listEventChannels();
            if any(ismember(fieldsAffected, eventFields))
                needUpdate = true; % force update of .valid below
                td = td.updateEventData(intersect(fieldsAffected, eventFields));
                % this will also flush alignSummarySet and reapply the
                % alignDescriptors to trial data

            end

            % check whether the affected fields matter for the
            if any(ismember(fieldsAffected, td.conditionInfo.attributeNames))
                % condition descriptor affected
                td = td.setConditionDescriptor(td.conditionInfo);
                % this already calls update valid so we don't need to do it again

            elseif needUpdate
                td = td.invalidateValid();
            end
        end

        function td = invalidateEventCache(td)
            td.warnIfNoArgOut(nargout);
            % flush event data, counts, align summary
            td.odc = td.odc.copy();
            td.odc.flushEventData();
        end

        function td = reinitializeConditionAlign(td)
            td.warnIfNoArgOut(nargout);

            td = td.invalidateEventCache();

            % don't call update valid until both are updated
            for i = 1:td.nAlign
                td.alignInfoSet{i} = td.alignInfoSet{i}.applyToTrialData(td);
            end

            % will call update valid
            td = td.setConditionDescriptor(td.conditionInfo);
        end

         % synchronize valid between AlignInfo and ConditionINfo
         % shouldn't need to call this manually, but just in case
        function td = invalidateValid(td)
            td.warnIfNoArgOut(nargout);
            % this will flush .valid and .invalidCause
            td = invalidateValid@TrialData(td);

            if ~td.initialized
                % this gets called during initialization, so we return
                % early here
                return;
            end

            % reset the manual valid arrays of all condition info and align
            % info instances so that no attributes are overlooked
            td.conditionInfo = td.conditionInfo.setManualValidTo(td.manualValid);
            for iA = 1:td.nAlign
                td.alignInfoSet{iA} = td.alignInfoSet{iA}.setManualValidTo(td.manualValid);
            end

            % recompute valid using the conjunction condition info and
            % align info (this happens inside buildValid, which will be
            % called when we request .valid for the first time)
            valid = td.valid;

            % coordinate the valid arrays of all condition info and align
            % info instances
            td.conditionInfo = td.conditionInfo.setManualValidTo(valid);
            for iA = 1:td.nAlign
                td.alignInfoSet{iA} = td.alignInfoSet{iA}.setManualValidTo(valid);
            end

            % cause align summary to be recomputed
            td.alignSummarySet = [];
        end

        % print a short description
        function disp(td)
            if numel(td) > 1
                for i = 1:numel(td)
                    td(i).printDescriptionShort();
                end
                return
            end
            
            td.printDescriptionShort();

            td.conditionInfo.printDescription();
            if td.hasRandomizationSpecified
                tcprintf('inline', '  {bright blue}Randomization via: {none}%d samples, %s, seed %d\n', td.nRandomized, td.randomizationDescription, td.randomSeed);
            end
            for iA = 1:td.nAlign
                td.alignInfoSet{iA}.printDescription('active', td.alignInfoActiveIdx == iA);
            end

            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end

        % generates a short description of each trial
        function desc = getTrialDescriptions(td, varargin)
            p = inputParser();
            p.addParameter('includeAttributes', true, @islogical);
            p.addParameter('multiline', true, @islogical);
            p.parse(varargin{:});

            % builds a description of a given trial
            desc = cellvec(td.nTrials);
            valid = td.valid;
            cIdx = td.conditionIdx;
            cSubs = td.conditionSubs;
            cNames = cellvec(td.nTrials);
            cNames(valid) = td.conditionNames(cIdx(valid));

            if p.Results.multiline
                sep = newline; % newline
            else
                sep = ' ';
            end

            if p.Results.includeAttributes
                % include values of condition related params and
                % manually include trialDescriptionExtraParams
                attr = union(td.trialDescriptionExtraParams, td.attributeParams);
                valueStrings = td.getParamMultiAsStrings(attr, 'separator', sep, 'includeParamNames', true);
            else
                valueStrings = cellvec(td.nTrials);
            end

            for i = 1:td.nTrials
                if valid(i)
                    validStr = '';
                else
                    validStr = ' (invalid)';
                end
                lines = {};
                lines{1} = sprintf('Trial %d%s', i, validStr);
                if valid(i)
                    lines{2} = sprintf('Condition %d (%s): %s', cIdx(i), TrialDataUtilities.String.strjoin(cSubs(i, :), ','), cNames{i});
                    lines{3} = valueStrings{i};
                end

                desc{i} = TrialDataUtilities.String.strjoin(lines, sep);
            end
        end

        function descByGroup = getTrialDescriptionsGrouped(td, varargin)
            desc = td.getTrialDescriptions(varargin{:});
            descByGroup = td.groupElements(desc);
        end
    end

    methods % Low-level add and remove channels. Probably don't want to use addChannel directly
        function name = generateTemporaryChannelName(td, start)
            if nargin < 2 || isempty(start)
                start = 'temporary';
            end

            name = matlab.lang.makeUniqueStrings(...
                    matlab.lang.makeValidName(start, 'ReplacementStyle', 'delete'), ...
                    fieldnames(td.data));
        end


        function td = addChannel(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = addChannel@TrialData(td, varargin{:});

            % Don't do this anymore, we no longer auto-add attributes to
            % conditionInfo until they are needed

%             % detect whether any new condition info compatible params
%             % have been added
%             namesOld = td.listConditionInfoCompatibleParamChannels();
%             td = addChannel@TrialData(td, varargin{:});
%             names = td.listConditionInfoCompatibleParamChannels();
%
%             % if so, add them to the condition info with valueLists
%             % specified
%             newAttr = setdiff(names, namesOld);
%             for iA = 1:numel(newAttr)
%                 td.conditionInfo = td.conditionInfo.addAttribute(newAttr{iA}, ...
%                     'values', td.getParam(newAttr{iA}));
%             end
        end

        function td = dropChannels(td, names)
            names = string(names);

            % don't remove special channels
            names = setdiff(names, td.listSpecialChannels());

            % don't remove channels that don't exist
            names = intersect(names, td.listChannels('includeNamedSubChannels', true));
            if isempty(names)
                return;
            end

            % check whether any of the alignInfo events and error if so
            for iA = 1:td.nAlign
                alignEvents = td.alignInfoSet{iA}.getEventList();
                mask = ismember(alignEvents, names);
                if any(mask)
                    error('TrialData alignment depends on event %s', ...
                        TrialDataUtilities.String.strjoin(alignEvents(mask)));
                end
            end

            % check whether any of the events are in condition info
            conditionParams = td.conditionInfo.attributeNames;
            mask = ismember(conditionParams, names);
            if any(mask)
                error('TrialData conditioning depends on params %s', ...
                    TrialDataUtilities.String.strjoin(conditionParams(mask)));
            end

            td = dropChannels@TrialData(td, names);

            % in case we lost some event channels, update the alignInfo,
            % which has all of the events upfront. No need to update
            % conditionInfo, since it receives the param data as needed and
            % we've already checked that it isn't using any of the dropped
            % channels
            td = td.applyAlignInfoSet();
        end

        function td = dropNonConditionAlignChannelsExcept(td, names)
            % drop all channels except specified and param/events that could be
            % used for condition grouping and alignment

            if ~iscell(names)
                names = {names};
            end

            td.warnIfNoArgOut(nargout);
            td = td.dropChannelsExcept([td.listEventChannels(); td.listParamChannels(); makecol(names)]);
        end

        function td = setChannelUnitsPrimary(td, name, units)
            td.warnIfNoArgOut(nargout);
            td = setChannelUnitsPrimary@TrialData(td, name, units);

            % if the channel is an attribute channel, we need to tell
            % conditionInfo that its units have changed
            if ismember(name, td.conditionInfo.attributeNames)
                td.conditionInfo = td.conditionInfo.setAttributeUnits(name, units);
            end

            % we shouldn't need to call postUpdateConditionInfo, since this
            % should nominally only change the condition names, not
            % anything related to .valid or the condition structure
        end

        function td = trimAllChannelsToCurrentAlign(td)
            td.warnIfNoArgOut(nargout);
            [startTimes, stopTimes] = td.getTimeStartStopEachTrial();
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            td = td.trimAllChannelsRaw(startTimes + offsets, stopTimes + offsets);
        end

        function td = trimAllChannelsToMaxDuration(td, duration)
            td.warnIfNoArgOut(nargout);
            td = td.unalign().start('TrialStart').stop('TrialStart', duration);
            td = td.trimAllChannelsToCurrentAlign();
            td = td.unalign();
        end
    end

    % Copying settings from another tdca
    methods
        function td = setupLike(td, tdOther)
            td.warnIfNoArgOut(nargout);
            td = td.setConditionDescriptor(tdOther.conditionDescriptor);
            td = td.setAlignDescriptorSet(tdOther.alignInfoSet);
        end
    end

    % ConditionInfo control
    methods
        % parameters that are either scalar or strings
        function names = listConditionInfoCompatibleParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor') && ...
                (cd.isScalarByField(1) || cd.isStringByField(1)), ...
                channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end

        function paramStruct = getConditionInfoCompatibleParamStruct(td)
            names = td.listConditionInfoCompatibleParamChannels();
            paramStruct = keepfields(td.data, names);
        end

        function td = initializeConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            if isempty(td.conditionInfo)
                paramStruct = emptyStructArray(td.nTrials);
                td.conditionInfo = ConditionInfo.fromStruct(paramStruct);
            end
        end

        % given a cellvec or nmeric vector, group its elements
        function varargout = groupElements(td, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) TensorUtils.selectAlongDimension(data, 1, idx, false), ...
                    td.listByCondition, 'UniformOutput', false);
%                 varargout{i} = cellfun(@(idx) data(idx,:), td.listByCondition, ...
%                     'UniformOutput', false);
            end
        end

        % given a cellvec or nmeric vector, group its elements along
        % dimension 1 but without preserving the conditions tensor shape
        function varargout = groupElementsFlat(td, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) TensorUtils.selectAlongDimension(data, 1, idx, false), ...
                    td.listByCondition(:), 'UniformOutput', false);
            end
        end

        % given data with dimension 1 with size nTrials, group by condition
        % and map out{i} = fn(group{i})
        function out = mapByGroup(td, fn, varargin)
            dataByGroup = td.groupElements(varargin{:});
            out = cellfun(fn, dataByGroup, 'UniformOutput', false);
        end

%         function [means, sem, std, nTrials, totalWeights] = computeGroupMeansWeighted(td, data, minTrials, minTrialFraction)
%             deal(nan(td.nConditions, numel(tvec), nUnits));
%             for iC = 1:td.nConditions
%                 if ~isempty(data{iC})
%                     [means(iC, :, :), sem(iC, :, :), nTrials(iC, :, :), std(iC, :, :)] = ...
%                         TrialDataUtilities.Data.nanMeanSemMinCount(rateCell{iC}, 1, minTrials, minTrialFraction);
%                 end
%             end
%         end

        function td = setConditionDescriptor(td, cd)
            td.warnIfNoArgOut(nargout);

            % manually accept a condition descriptor instance
            if isa(cd, 'ConditionInfo')
                cd = ConditionDescriptor.fromConditionDescriptor(cd);
            end
            assert(isequal(class(cd), 'ConditionDescriptor'), 'Must be a ConditionDescriptor instance');

%             % grab the param data to feed to the condition descriptor
%             if isempty(cd.attributeNames)
%                 paramData = emptyStructArray(td.nTrials);
%             else
%                 paramData = td.getRawChannelDataAsStruct(cd.attributeNames);
%             end

            % build condition info from condition descriptor
            td.conditionInfo = ConditionInfo.fromConditionDescriptor(cd, td);
            td = td.postUpdateConditionInfo();
        end

        function td = setConditionAppearanceFn(td, fn)
            % Update the appearanceFn callback of conditionDescriptor
            % without invalidating any of the other cached info
            if nargin < 2
                fn = [];
            end
            assert(isempty(fn) || isa(fn, 'function_handle'));
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.appearanceFn = fn;
        end

        function appearances = testConditionAppearanceFn(td)
            appearances = td.conditionInfo.testConditionAppearanceFn();
        end

        function td = resetConditionAppearanceFn(td)
            % Update the appearanceFn callback of conditionDescriptor
            % without invalidating any of the other cached info
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.appearanceFn = [];
        end

        function td = freezeConditionAppearances(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.freezeAppearances();
        end

        function td = setConditionNameFn(td, fn)
            % Update the nameFn callback of conditionDescriptor
            % without invalidating any of the other cached info
            assert(isa(fn, 'function_handle'));
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.nameFn = fn;
        end
        
        function td = clearAppearanceModifications(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.clearAppearanceModifications(varargin{:});
        end

        function td = colorByAttributes(td, varargin)
            td.warnIfNoArgOut(nargout);

            % add any needed attributes to condition info
            td = td.addAttribute(varargin{1});

            td.conditionInfo = td.conditionInfo.colorByAttributes(varargin{:});
        end

        function td = colorByAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.colorByAxes(varargin{:});
        end

        function td = lineWidthByAttributes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.lineStyleByAttributes(varargin{:});
        end

        function td = lineWidthByAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.lineWidthByAxes(varargin{:});
        end

        function td = lineStyleByAttributes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.lineStyleByAttributes(varargin{:});
        end

        function td = lineStyleByAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.lineStyleByAxes(varargin{:});
        end

        function setAxisColorOrderToConditionColors(td, axh)
            if nargin < 2
                axh = gca;
            end
            set(axh, 'ColorOrder', td.conditionColors, 'ColorOrderIndex', 1);
        end

        function cellOfCategoricals = getConditionSubsAsCategorical(td)
           valLists = td.axisValueListsAsStringsShort;
           subs = td.conditionSubs;
           cellOfCategoricals = cellvec(td.nAxes);
           for iAx = 1:td.nAxes
               cellOfCategoricals{iAx} = categorical(subs(:, iAx), 1:numel(valLists{iAx}), valLists{iAx});
           end
        end

        function td = selectTrials(td, mask)
            % Apply a logical mask or index selection to the list of trials
            % within, appropriately notifying the condition descriptor and
            % align descriptor within
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.selectTrials(mask);
            if ~isempty(td.odc.conditionInfoRandomized)
                % only slice conditionInfoRandomized if it hasn't already
                % been built on the fly
                td.conditionInfoRandomized = td.conditionInfoRandomized.selectTrials(mask);
            end
            for i = 1:td.nAlign
                td.alignInfoSet{i} = td.alignInfoSet{i}.selectTrials(mask);
            end

            % select from cached event data as well
            if ~isempty(td.odc.eventData)
                c = td.odc.copy();

                flds = fieldnames(c.eventData);
                for iFld = 1:numel(flds)
                    fld = flds{iFld};
                    c.eventData.(fld) = c.eventData.(fld)(mask, :);
                    c.eventCounts.(fld) = c.eventCounts.(fld)(mask, :);
                end

                td.odc = c;
            end

            % update this after calling everythign above, since it will
            % call update valid and requires .valid lengths to be consistent
            % in the align descriptors and condition descriptor once that happens
            td = selectTrials@TrialData(td, mask);
        end

        function [td, namesModified] = addAttribute(td, names)
            % add attributes in names that aren't already in ConditionInfo
            % for non-param channels, the name of the attribute might be
            % modified to reflect the specific way in which the channel was
            % used to create the scalar attribute value
            %
            % for event channels, the first event occurrence time will be
            % used, relative to the current alignment zero.

            if isstringlike(names)
                wasChar = true;
                names = {names};
            else
                wasChar = false;
            end
            namesModified = names;

            for i = 1:numel(names)
                name = names{i};
                if ismember(name, td.conditionInfo.attributeNames)
                    continue;
                end

                cd = td.getChannelDescriptor(name);
                if isa(cd, 'EventChannelDescriptor')
                    % the event time will be referenced from the current
                    % zero event
                    nameMod = matlab.lang.makeValidName(sprintf("%s_from_%s", name, td.alignInfoActive.zeroLabel));
                    values = td.getEventFirst(name);

                elseif isa(cd, 'AnalogChannelDescriptor')
                    % take analog sample at zero time
                    nameMod = matlab.lang.makeValidName(sprintf("%s_at_%s", name, td.alignInfoActive.zeroLabel));
                    values = td.getAnalogSample(name);

                elseif isa(cd, 'ParamChannelDescriptor')
                    % params get used as is
                    nameMod = name;
                    values = td.getParamRaw(name);
                else
                    error('Unable to create attribute for channel type %s', class(cd));
                end

                units = td.getChannelUnitsPrimary(name);
                td.conditionInfo = td.conditionInfo.addAttribute(nameMod, 'values', values, 'units', units);
                namesModified{i} = nameMod;
            end

            if wasChar
                namesModified = namesModified{1};
            end

            td = td.postUpdateConditionInfo();
        end

        function td = groupBy(td, varargin)
            td.warnIfNoArgOut(nargout);

            % add any needed attributes to condition info
            for i = 1:numel(varargin)
                td = td.addAttribute(varargin{i});
            end

            td.conditionInfo = td.conditionInfo.groupBy(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = addAxis(td, attrList, varargin)
            td.warnIfNoArgOut(nargout);

            % add any needed attributes to condition info
            td = td.addAttribute(attrList);

            td.conditionInfo = td.conditionInfo.addAxis(attrList, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = removeAxis(td, attrList, varargin)
            td.warnIfNoArgOut(nargout);

            td.conditionInfo = td.conditionInfo.removeAxis(attrList, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = ungroup(td)
            % this only undoes the grouping axes, NOT the value list
            % filtering. use reset condition info for that
            td.warnIfNoArgOut(nargout);
            td = td.groupBy();
        end

        function td = ungroupMaintainValid(td)
            td.warnIfNoArgOut(nargout);
            td = td.markTrialsTemporarilyInvalid(~td.valid, 'Ungroup Fix Valid');
            td = td.groupBy();
        end

        function td = fixAllAxisValueLists(td)
            % this only undoes the grouping axes, NOT the value list
            % filtering. use reset condition info for that
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.fixAllAxisValueLists();
            td = td.postUpdateConditionInfo();
        end

        function td = reset(td)
            td.warnIfNoArgOut(nargout);
            td = reset@TrialData(td); % important to clear temporarily invalid
            td = td.resetConditionInfo();
            td = td.unalign();
        end

        function td = resetHard(td)
            % clears temporary and permanent invalid and resets to all trials being valid
            td.warnIfNoArgOut(nargout);
            td = resetHard@TrialData(td);
            td = td.resetConditionInfo();
            td = td.unalign();
        end

        function td = setConditionIncludeMask(td, m)
            td.warnIfNoArgOut(nargout);
            td = td.fixAllAxisValueLists();
            td.conditionInfo = td.conditionInfo.setConditionIncludeMask(m);
            td = td.postUpdateConditionInfo();
        end

        function td = clearConditionIncludeMask(td)
            % this will leave the axis value lists fixed
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.resetConditionIncludeMask();
            td = td.postUpdateConditionInfo();
        end

        function td = withConditions(td, m)
            td.warnIfNoArgOut(nargout);
            td = td.setConditionIncludeMask(m);
        end

        function td = withAllConditions(td)
            td.warnIfNoArgOut(nargout);
            td = td.clearConditionIncludeMask();
        end

        function td = withFirstNTrialsEachCondition(td, n)
            td.warnIfNoArgOut(nargout);
            lists = td.listByCondition;
            listsDrop = cell(size(lists));
            for iC = 1:td.nConditions
                if numel(lists{iC}) > n
                    listsDrop{iC} = lists{iC}(n+1:end);
                end
            end

            drop = cat(1, listsDrop{:});
            td = td.markTrialsTemporarilyInvalid(drop, 'withFirstNTrialsEachCondition');
        end
        
        function td = withFirstNValidTrials(td, n)
            td.warnIfNoArgOut(nargout);
            
            drop = true(td.nTrials, 1);
            inds = find(td.valid, n, 'first');
            drop(inds) = false;
            td = td.markTrialsTemporarilyInvalid(drop, 'withFirstNValidTrials');
        end

        function td = withTrials(td, mask, desc)
            td.warnIfNoArgOut(nargout);
            if nargin < 3
                desc = 'withTrials';
            end
            % temporarily keep only trial in mask
            keep = false(td.nTrials, 1);
            keep(mask) = true;
            td = td.markTrialsTemporarilyInvalid(~keep, desc);
        end

        function td = withRandomSubsetOfTrials(td, n)
            td.warnIfNoArgOut(nargout);
            % temporarily keep only n of the currently valid trials
            inds = find(td.valid);
            keep = randsample(inds, n, false);
            td = td.withTrials(keep, 'withRandomSubsetOfTrials');
        end

        function td = reshapeAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.reshapeAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = permuteAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.permuteAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = transposeAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.transposeAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = flattenAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.flattenAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = flattenAxesExceptFirst(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.flattenAxesExceptFirst(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = sortWithinConditionsBy(td, attrList, varargin)
            td.warnIfNoArgOut(nargout);

            if isstringlike(attrList)
                attrList = {attrList};
            end

            attrListModified = cellvec(numel(attrList));
            for i = 1:numel(attrList)
                if strncmp(attrList{i}, '-', 1)
                    [td, attrListModified{i}] = td.addAttribute(attrList{i}(2:end));
                else
                    [td, attrListModified{i}] = td.addAttribute(attrList{i});
                end
            end

            td.conditionInfo = td.conditionInfo.sortWithinConditionsBy(attrListModified, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function [td, maskC] = selectConditions(td, varargin)
            % select specific conditions by linear index or mask
            % and return a single-axis condition descriptor with just those
            % conditions selected
            td.warnIfNoArgOut(nargout);
            [td.conditionInfo, maskC] = td.conditionInfo.selectConditions(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function [td, maskC] = selectConditionsAlongAxis(td, varargin)
            td.warnIfNoArgOut(nargout);
            [td.conditionInfo, maskC] = td.conditionInfo.selectConditionsAlongAxis(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function [td, maskC] = matchSelectConditionsAlongAxis(td, varargin)
            td.warnIfNoArgOut(nargout);
            [td.conditionInfo, maskC] = td.conditionInfo.matchSelectConditionsAlongAxis(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = resetConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.setConditionDescriptor(ConditionDescriptor());
        end

        function td = setAllAttributeValueListsAuto(td)
            % will undo any filtering by attribute value lists and restore
            % auto-binning
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAllAttributeValueListsAuto();
            td = td.postUpdateConditionInfo();
        end

        function td = setAttributeValueListAuto(td, attr)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(attr);
            td.conditionInfo = td.conditionInfo.setAttributeValueListAuto(attr);
            td = td.postUpdateConditionInfo();
        end

        function td = binAttribute(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttribute(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = filterRange(td, name, range)
            td.warnIfNoArgOut(nargout);
            td = td.binAttribute(name, range);
        end

        function td = filterRangeExclude(td, attrName, rangeExclude)
            td.warnIfNoArgOut(nargout);
            td = td.binAttribute(attrName, [-Inf rangeExclude(1); rangeExclude(2) Inf]);
        end

        function td = binAttributeUniform(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttributeUniform(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = binAttributeQuantiles(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttributeQuantiles(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = setAttributeValueList(td, attrName, valueList, varargin)
            td.warnIfNoArgOut(nargout);
            attrName = char(attrName);
            td = td.addAttribute(attrName);
            if td.isChannelCategorical(attrName)
                valueList = string(valueList);
                valueList = categorical(valueList);
            end
            td.conditionInfo = td.conditionInfo.setAttributeValueList(attrName, valueList, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = setAttributeValueListDisplayAs(td, attrName, displayAs)
            td.warnIfNoArgOut(nargout);
            if isstringlike(attrName)
                td = td.addAttribute(attrName);
            end
            td.conditionInfo = td.conditionInfo.setAttributeValueListDisplayAs(attrName, displayAs);
            td = td.postUpdateConditionInfo();
        end

        function td = setAttributeDisplayAs(td, attrName, displayAs)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(attrName);
            td.conditionInfo = td.conditionInfo.setAttributeDisplayAs(attrName, displayAs);
            td = td.postUpdateConditionInfo();
        end

        function td = fixAttributeValueList(td, attrName)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.fixAttributeValueList(attrName);
            td = td.postUpdateConditionInfo();
        end

        function td = fixAllAttributeValueLists(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.fixAllAttributeValueLists();
            td = td.postUpdateConditionInfo();
        end

        function td = filter(td, attrName, valueList)
            % td = filter(td, attr, valueList)
            % shortcut for setAttributeValueList
            td.warnIfNoArgOut(nargout);
            attrName = char(attrName);
            td = td.setAttributeValueList(attrName, valueList);
        end

        function td = filterExclude(td, attrName, valueListExclude)
            td.warnIfNoArgOut(nargout);
            valueList = setdiff(td.getParamUnique(attrName), valueListExclude);
            td = td.setAttributeValueList(attrName, valueList);
        end

        function td = noFilter(td, attrName)
            td.warnIfNoArgOut(nargout);
            td = td.setAttributeValueListAuto(attrName);
        end

        function valueList = getAttributeValueList(td, attrName)
            td = td.addAttribute(attrName);
            valueList = td.conditionInfo.getAttributeValueList(attrName);
        end

        function td = setAxisValueList(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueList(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = setAxisValueListAutoAll(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueListAutoAll(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = setAxisValueListAutoOccupied(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueListAutoOccupied(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function valueList = getAxisValueList(td, varargin)
            valueList = td.conditionInfo.getAxisValueList(varargin{:});
        end

        function td = setAxisValueListDisplayAs(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueListAsStrings(varargin{:});
            % no need to update condition info since the conditions and
            % validity won't change
        end

        function td = setAxisValueListDisplayAsShort(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueListAsStringsShort(varargin{:});
            % no need to update condition info since the conditions and
            % validity won't change
        end

        function td = setConditionNames(td, varargin)
            % setConditionNames(names, [namesShort]);
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setConditionNames(varargin{:});
            % no need to update condition info since the conditions and
            % validity won't change
        end

        function td = postUpdateConditionInfo(td, clearRandomized)
            if nargin < 2
                clearRandomized = true; % sometimes don't want to clear conditionInfoRandomized when using .withRandomized( ) to access a specific set of data
            end
            td.warnIfNoArgOut(nargout);
            td = td.invalidateValid();

            if clearRandomized
                % flush conditionInfo randomized, needs to be recreated
                td.conditionInfoRandomized = [];

                % flush cached randomize lists if found
                td.odc = td.odc.copy();
                td.odc.flushRandomized();
            end
        end

        % filter trials that are valid based on ConditionInfo
        function td = selectValidTrialsConditionInfo(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.selectTrials(td.conditionInfo.valid);
        end

        % filter trials where a specific attribute matches a specific value
        % pass-thru to ConditionInfo
        function td = filteredByAttribute(td, varargin)
            % call on ConditionInfo, and filter based on its valid trials
            td.warnIfNoArgOut(nargout);
            % no need to copy, filteredBy already copies the conditionInfo
            td.conditionInfo = td.conditionInfo.filteredByAttribute(varargin{:});
            td = td.postUpdateConditionInfo();
%             td = td.filterValidTrialsConditionInfo();
        end

        % filter trials based on matching attribute values,
        % pass thru to ConditionInfo
        function td = filteredByAttributeStruct(td, varargin)
            % call on ConditionInfo, and filter based on its valid trials
            td.warnIfNoArgOut(nargout);
            % no need to copy, filteredBy already copies the conditionInfo
            td.conditionInfo = td.conditionInfo.filteredByAttributeStruct(varargin{:});
            td = td.postUpdateConditionInfo();
%             td = td.filterValidTrialsConditionInfo();
        end

        % manual setting of condition membership
        function td = setConditionMembershipManual(td, conditionMembership, reasonInvalid, varargin)
            td.warnIfNoArgOut(nargout);

            if nargin < 3
                reasonInvalid = {};
            end

            td.conditionInfo = td.conditionInfo.fixAllAxisValueLists();
            td.conditionInfo = td.conditionInfo.setConditionMembershipManual(conditionMembership, reasonInvalid, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = setConditionMembershipManualFn(td, conditionMembershipFn, varargin)
            td.warnIfNoArgOut(nargout);

            td.conditionInfo = td.conditionInfo.setConditionMembershipManualFn(conditionMembershipFn, varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function tdCell = getTrialDataByCondition(td)
            tdCell = arrayfun(@(condition) td.filteredByAttributeStruct(condition), ...
                td.conditions, 'UniformOutput', false);
        end

        function v = get.conditions(td)
            v = td.conditionInfo.conditions;
        end

        function v = get.conditionsAxisAttributesOnly(td)
            v = td.conditionInfo.conditionsAxisAttributesOnly;
        end

        function v = get.conditionsAsStrings(td)
            v = td.conditionInfo.conditionsAsStrings;
        end

        function tf = get.isGrouped(td)
            tf = td.conditionInfo.nAxes > 0;
        end

        function n = get.nConditions(td)
            n = td.conditionInfo.nConditions;
        end

        function sz = get.conditionsSize(td)
            sz = td.conditionInfo.conditionsSize;
        end

        function sz = get.conditionsSizeNoExpand(td)
            sz = td.conditionInfo.conditionsSizeNoExpand;
        end

        function m = get.conditionIncludeMask(td)
            m = td.conditionInfo.conditionIncludeMask;
        end

        function list = get.attributeParams(td)
            list = td.conditionInfo.attributeNames;
        end

        function n = get.nAxes(td)
            n = td.conditionInfo.nAxes;
        end

        function v = get.axisValueLists(td)
            v = td.conditionInfo.axisValueLists;
        end

        function v = get.axisValueListsAsStrings(td)
            v = td.conditionInfo.axisValueListsAsStrings;
        end

        function v = get.axisValueListsAsStringsShort(td)
            v = td.conditionInfo.axisValueListsAsStringsShort;
        end

        function v = get.axisNames(td)
            v = td.conditionInfo.axisNames;
        end

        function v = get.axisDescriptions(td)
            v = td.conditionInfo.axisDescriptions;
        end

        function v = get.conditionNames(td)
            v = td.conditionInfo.names;
        end

        function v = get.conditionNamesShort(td)
            v = td.conditionInfo.namesShort;
        end

        function v = get.conditionNamesMultiline(td)
            v = td.conditionInfo.namesMultiline;
        end

        function v = get.conditionNamesMultilineShort(td)
            v = td.conditionInfo.namesMultilineShort;
        end

        function v = get.conditionAppearances(td)
            v = td.conditionInfo.appearances;
        end

        function td = set.conditionAppearanceFn(td, v)
            td.conditionInfo.appearanceFn = v;
        end

        function v = get.conditionAppearanceFn(td)
            v = td.conditionInfo.appearanceFn;
        end

        function v = get.conditionColors(td)
            v = cat(1, td.conditionInfo.appearances.Color);
        end

        function v = get.conditionNameFn(td)
            v = td.conditionInfo.nameFn;
        end

        function v = get.listByCondition(td)
            v = td.conditionInfo.listByCondition;
        end

        function v = get.listByConditionWeights(td)
            v = td.conditionInfo.listByConditionWeights;
        end

        function v = get.nTrialsByCondition(td)
            v = cellfun(@numel, td.listByCondition);
        end

        function v = get.conditionIdx(td)
            v = td.conditionInfo.conditionIdx;
        end

        function v = get.conditionNamesEachTrial(td)
            idx = td.conditionIdx;
            namesFlat = td.conditionNames(:);
            namesValidTrials = namesFlat(idx(td.valid));
            v = TensorUtils.inflateMaskedTensor(namesValidTrials, 1, td.valid, '');
        end

        function v = get.conditionSubs(td)
            v = td.conditionInfo.conditionSubs;
        end

        function v = get.conditionMembership(td)
            v = td.conditionInfo.conditionMembership;
        end

        function minTimeDelta = get.minTimeDelta(td)
            minTimeDelta = td.alignDescriptor.minTimeDelta;
        end
    end

    methods(Hidden)
        function cdSub = expandConditionDescriptorForSubtraction(td, cdSub)
            if isa(cdSub, 'ConditionInfo')
                cdSub = ConditionInfo.fromConditionDescriptor(cdSub, td);
            end

            if ~isequal(td.conditionInfo.conditionsSizeNoExpand, cdSub.conditionsSizeNoExpand) || ...
               ~isequal(td.conditionInfo.axisAttributes, cdSub.axisAttributes)
                debug('Using ConditionDescriptor.expandToMatch to match condition descriptor sizes\n');
                [cdMatched, cdModifiedMask] = ConditionDescriptor.expandToMatch({td.conditionDescriptor, cdSub});
                if cdModifiedMask(1)
                    error('ConditionDescriptor.expandToMatch required modification of the current conditionDescriptor to match the provided subtractConditionDescriptor');
                end

                cdSub = cdMatched{2};
            end
        end
    end

    % AlignInfo control
    methods
        function td = initializeAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            if isempty(td.alignInfoSet)
                td.alignInfoSet = {AlignInfo()};
                td = td.applyAlignInfoSet();
            end
            if isempty(td.alignInfoActiveIdx)
                td.alignInfoActiveIdx = 1;
            end
        end

        function td = applyAlignInfoSet(td)
            td.warnIfNoArgOut(nargout);
            for i = 1:td.nAlign
                td.alignInfoSet{i} = td.alignInfoSet{i}.applyToTrialData(td);
            end
            
            td.interAlignGaps = [];
            td.manualAlignTimeOffsets = [];
            td = td.postUpdateAlignInfo();
        end

        function td = postUpdateAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.invalidateValid();
        end

        function td = align(td, varargin)
            td.warnIfNoArgOut(nargout);

            adSet = cell(numel(varargin), 1);
            for i = 1:numel(varargin)
                ad = varargin{i};

                if iscell(ad)
                    error('Please provide alignDescriptors as successive arguments');
                end
                if isstringlike(ad)
                    adSet{i} = AlignInfo(string(ad));
                else
                    if isa(ad, 'AlignDescriptor')
                        % convert to AlignInfo
                        adSet{i} = AlignInfo.fromAlignDescriptor(ad);
                    elseif isa(ad, 'AlignInfo')
                        adSet{i} = ad;
                    end
                end
            end

            td.alignInfoSet = adSet;
            td.alignInfoActiveIdx = 1;
            td = td.applyAlignInfoSet();
        end

        function td = setAlignDescriptor(td, ad)
            td.warnIfNoArgOut(nargout);
            td = td.align(ad);
        end

        function td = setAlignDescriptorSet(td, adSet)
            td.warnIfNoArgOut(nargout);
            td = td.align(adSet{1});
            for iA = 2:numel(adSet)
                td = td.addAlign(adSet{iA});
            end
            td = td.useAlign(1);
        end

        function td = addAlign(td, varargin)
            % whereas align replaces the existing align descriptor set,
            % this appends an align descriptor
            td.warnIfNoArgOut(nargout);

            adSet = cat(1, td.alignInfoSet, makecol(varargin));
            td = td.align(adSet{:});

            % and for convenience switch to make this new alignment active
            td = td.useAlign(td.nAlign);
        end

        function td = selectAlign(td, mask)
            td = td.align(td.alignInfoSet{mask});
        end

        function td = unalign(td)
            td.warnIfNoArgOut(nargout);
            td = td.align('TrialStart:TrialEnd @ TrialStart');
        end

        function td = setInterAlignGap(td, gaps)
            % set .interAlignGaps, which represent the time gaps between
            % successive alignments, mainly when plotting
            td.warnIfNoArgOut(nargout);

            if td.nAlign < 2
                error('Inter alignment gap not valid when only one align present');
            end
            if isscalar(gaps)
                gaps = repmat(gaps, td.nAlign - 1, 1);
            else
                assert(numel(gaps) == td.nAlign - 1, 'Gaps must be scalar or be length nAlign-1');
            end

            td.interAlignGaps = gaps;
        end
        
        function td = setAlignTimeOffsets(td, offsets)
            % set .manualAlignTimeOffsets, which represent the time gaps between
            % successive alignments, mainly when plotting
            td.warnIfNoArgOut(nargout);
            
            if ~isempty(offsets)
                assert(numel(offsets) == td.nAlign, 'Offsets must be length nAlign');
                td.manualAlignTimeOffsets = offsets;
            end
        end

        % the following methods pass-thru to alignInfo:

        function td = pad(td, window, expand)
            % add a padding window to the AlignInfo
            % may change which trials are valid
            % usage: pad([pre post]) or pad(pre, post)
            % pre > 0 means add padding before the start (typical case)
            td.warnIfNoArgOut(nargout);

            if nargin < 3
                expand = false;
            end

            updated = false;

            if expand
                for iA = 1:td.nAlign
                    if td.alignInfoSet{iA}.padPre < window(1) || td.alignInfoSet{iA}.padPost < window(2)
                        window = max(window, [td.alignInfoSet{iA}.padPre td.alignInfoSet{iA}.padPost]);
                        td.alignInfoSet{iA} = td.alignInfoSet{iA}.pad(window);
                        updated = true;
                    end
                end
            else
                for iA = 1:td.nAlign
                    if td.alignInfoSet{iA}.padPre ~= window(1) || td.alignInfoSet{iA}.padPost ~= window(2)
                        td.alignInfoSet{iA} = td.alignInfoSet{iA}.pad(window);
                        updated = true;
                    end
                end
            end

            if updated
                % to synchronize any changes in alignment validity
                td = td.postUpdateAlignInfo();
            end
        end

        function td = padForTimeBinning(td, timeDelta, binAlignmentMode, expand, includeExtraBin)
            % pad the edges of the alignment by the appropriate amount to
            % include data needed for time bin sampling
            td.warnIfNoArgOut(nargout);

            if nargin < 4
                expand = false;
            end
            if nargin < 5
                includeExtraBin = false;
            end

            switch binAlignmentMode
                case BinAlignmentMode.Acausal
                    window = [0 timeDelta];
                case BinAlignmentMode.Causal
                    window = [timeDelta 0];
                case BinAlignmentMode.Centered
                    window = [timeDelta/2 timeDelta/2];
                otherwise
                    error('Unknown binAlignmentMode');
            end

            if includeExtraBin
                window(1) = window(1) + timeDelta;
                window(2) = window(2) + timeDelta;
            end

            td = td.pad(window, expand);
        end

        function td = padForSpikeFilter(td, sf)
            % Pad trial data alignment for spike filter, including pre and
            % post windo of the filter, as well as the bin width used
            td.warnIfNoArgOut(nargout);
            td = td.pad(sf.padWindow);
        end

        function td = adjustStartStopAroundPadding(td)
            % replace current padding with actual shifts in start and stop times
            td.warnIfNoArgOut(nargout);
            ai = td.alignInfoActive;
            ai = ai.adjustStart('offset', ai.startOffset - ai.padWindow(1));
            ai = ai.adjustStop('offset', ai.stopOffset + ai.padWindow(2));
            ai = ai.pad(0, 0);
            td = td.setAlignDescriptor(ai);
        end

        function td = round(td, varargin)
            td.warnIfNoArgOut(nargout);
            for iAlign = 1:td.nAlign
                td.alignInfoSet{iAlign} = td.alignInfoSet{iAlign}.round(varargin{:});
            end
        end

        function td = noRound(td)
            td.warnIfNoArgOut(nargout);
            for iAlign = 1:td.nAlign
                td.alignInfoSet{iAlign} = td.alignInfoSet{iAlign}.noRound();
            end
        end

        function td = start(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.start(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = stop(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.stop(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = zero(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.zero(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = mark(td, eventStr, varargin)
            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParameter('index', [], @(x) isempty(x) || isstringlike(x) || isscalar(x));
            p.addParameter('as', AlignDescriptor.AUTO, @isstringlike);
            p.addParameter('color', [], @(x) isempty(x) || isstringlike(x) || isvector(x));
            p.addParameter('appear', [], @(x) isempty(x) || isa(x, 'AppearanceSpec'));
            p.addParameter('showOnData', true, @islogical);
            p.addParameter('showOnAxis', true, @islogical);
            p.addParameter('aggregateAllOccurrences', false, @islogical);
            p.parse(varargin{:});

            opts = p.Results;

            td.warnIfNoArgOut(nargout);

            eventName = td.alignInfoActive.parseEventOffsetString(eventStr, ...
                'mark', 'defaultIndex', ':');
            if td.hasAnalogChannelOrGroup(eventName)
                % make an event corresponding to this analog channel's
                % times so we can mark the sample times
                [td, eventField] = td.addEventFromAnalogTimes(eventName);
                eventStr = strrep(eventStr, eventName, eventField);

            elseif td.hasEventChannel(eventName)
                % an actual event, check for auto color
                cd = td.getChannelDescriptor(eventName);
                if isempty(opts.color)
                    opts.color = cd.color;
                end
            else
                error('Event %s not found', eventName);
            end

            % color arg must come last because of optional offset arg
            td.alignInfoActive = td.alignInfoActive.mark(eventStr, opts);
            td = td.postUpdateAlignInfo();
        end

        function td = removeMarkByString(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.removeMarkByString(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = removeMark(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.removeMark(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = interval(td, eventStart, eventStop, varargin)
            p = inputParser;
            p.addParameter('offsetStart', 0, @isscalar);
            p.addParameter('offsetStop', 0, @isscalar);
            p.addParameter('indexStart', ':', @(x) isstringlike(x) || isscalar(x));
            p.addParameter('indexStop', ':', @(x) isstringlike(x) || isscalar(x));
            p.addParameter('as', AlignDescriptor.AUTO, @isstringlike);
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.addParameter('showOnData', true, @islogical);
            p.addParameter('showOnAxis', true, @islogical);
            p.parse(varargin{:});
            opts = p.Results;

            td.warnIfNoArgOut(nargout);

            eventName = td.alignInfoActive.parseEventOffsetString(eventStart, 'interval start', 'defaultIndex', ':');
            if td.hasEventChannel(eventName)
                % an actual event, check for auto color
                cd = td.getChannelDescriptor(eventName);
                if isempty(opts.color)
                    opts.color = cd.color;
                end
            else
                error('Event %s not found', eventName);
            end

            td.alignInfoActive = td.alignInfoActive.interval(eventStart, eventStop, opts);
            td = td.postUpdateAlignInfo();
        end

        function td = intervalManual(td, intervalCell, as, varargin)
            % add an additional interval by adding the needed events with autogenerated names
            % intervalCell is nTrial x 1 cell of nOccur x 2 times, or is nOccur x 2 times matrix for all trials

            p = inputParser();
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.addParameter('showOnData', true, @islogical);
            p.addParameter('showOnAxis', true, @islogical);

            p.addParameter('eventStart', '', @isstringlike);
            p.addParameter('eventStop', '', @isstringlike);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            if isempty(p.Results.eventStart)
                eventStart = [as 'Start'];
            end
            eventStart = td.generateTemporaryChannelName(eventStart);
            if isempty(p.Results.eventStop)
                eventStop = [as 'Stop'];
            end
            eventStop = td.generateTemporaryChannelName(eventStop);

            % pull out start and stop times for each occurrence on each trial
            if ~iscell(intervalCell)
                assert(size(intervalCell, 2) == 2);
                starts = intervalCell(:, 1);
                stops = intervalCell(:, 2);
            else
                assert(numel(intervalCell) == td.nTrials);
                [starts, stops] = cellvec(td.nTrials);
                for iT = 1:td.nTrials
                    if ~isempty(intervalCell{iT})
                        assert(size(intervalCell{iT}, 2) == 2);

                        starts{iT} = intervalCell{iT}(:, 1);
                        stops{iT} = intervalCell{iT}(:, 2);
                    end
                end
            end

            td = td.addEvent(eventStart, starts);
            td = td.addEvent(eventStop, stops);

            extra = rmfield(p.Results, {'eventStart', 'eventStop'});
            td = td.interval(eventStart, eventStop, 'as', as, extra);
        end

        function td = removeInterval(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.removeInterval(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = setCurrentZero(td)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setCurrentZero();
            td = td.postUpdateAlignInfo();
        end
        
        function td = alignPreStart(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.buildPreStartAlign(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = alignPostStop(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.buildPostStopAlign(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = truncateBefore(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.truncateBefore(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = truncateAfter(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.truncateAfter(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = invalidateOverlap(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.invalidateOverlap(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = setOutsideOfTrialTruncate(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setOutsideOfTrialTruncate(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = setOutsideOfTrialInvalidate(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setOutsideOfTrialInvalidate(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = setOutsideOfTrialIgnore(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setOutsideOfTrialIgnore(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = alignRoundTimes(td, timeDelta)
            td.warnIfNoArgOut(nargout);
            if isempty(timeDelta)
                td.alignInfoActive = td.alignInfoActive.noRound();
            else
                td.alignInfoActive = td.alignInfoActive.round(timeDelta);
            end
            td = td.postUpdateAlignInfo();
        end

        % filter trials that are valid based on AlignInfo
        function td = filterValidTrialsAlignInfo(td, varargin)
            td.warnIfNoArgOut(nargout);
            avalid = truevec(td.nTrials);
            for iA = 1:td.nAlign
                avalid = avalid & td.alignInfoSet{iA}.computedValid;
            end
            td = td.selectTrials(avalid);
        end

        function td = abbrev(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.abbrev(varargin{:});
        end

        function td = lag(td, varargin)
            % shift so that aligned data is exctracted tDelay units earlier in time,
            % instead of start/stop at Event, start/stop at Event-tDelay. We must move
            % the zero so that the time reference for interpolation /
            % re-sampling remains fixed when we lag, otherwise we end up
            % with different numbers of samples. e.g. a signal extracted at
            % Move-100:+100 can be compared against a 50 ms lagged signal
            % which was taken from Move-150:+50.
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.lag(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function td = lagEachAlign(td, varargin)
            td.warnIfNoArgOut(nargout);
            for iA = 1:td.nAlign
                td.alignInfoSet{iA} = td.alignInfoSet{iA}.lag(varargin{:});
            end
            td = td.postUpdateAlignInfo();
        end

        function td = setStartAppearance(td, varargin)
            % updates the AppearanceSpec for start
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setStartAppearance(varargin{:});
        end

        function td = setStopAppearance(td, varargin)
            % updates the AppearanceSpec for stop
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setStopAppearance(varargin{:});
        end

        function td = setZeroAppearance(td, varargin)
            % updates the AppearanceSpec for zero
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setZeroAppearance(varargin{:});
        end

        function td = setMarkAppearance(td, varargin)
            % updates the AppearanceSpec for mark at index ind
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setMarkAppearance(varargin{:});
        end

        function td = setIntervalAppearance(td, varargin)
            % updates the AppearanceSpec for interval at index ind
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.setIntervalAppearance(varargin{:});
        end

        function td = clearMarksIntervals(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.clearMarksIntervals(varargin{:});
            td = td.postUpdateAlignInfo();
        end

        function ad = get.alignInfoActive(td)
            ad = td.alignInfoSet{td.alignInfoActiveIdx};
        end

        function td = set.alignInfoActive(td, ad)
            td.alignInfoSet{td.alignInfoActiveIdx} = ad;
        end

        function ad = get.alignSummaryActive(td)
            ad = td.alignSummarySet{td.alignInfoActiveIdx};
        end

        function td = set.alignSummaryActive(td, ad)
            td.alignSummarySet{td.alignInfoActiveIdx} = ad;
        end

        function td = useAlign(td, idx)
            td.warnIfNoArgOut(nargout);
            assert(isscalar(idx) && isnumeric(idx) && idx >= 1 && idx <= td.nAlign);
            td.alignInfoActiveIdx = idx;
        end

        function durations = getValidDurations(td)
            durations = td.alignInfoActive.getValidDurationByTrial();
            durations(~td.valid) = NaN;
        end

        function durations = getValidDurationsEachAlign(td)
            % get the time window for each trial
            % nTrials x nAlign
            durations = nan(td.nTrials, td.nAlign);
            for iA = 1:pset.nAlign
                durations(:, iA) = td.alignInfoSet{iA}.getValidDurationByTrial();
            end

            durations(~td.valid, :) = NaN;
        end

        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrial(td)
            % make relative to zero
            [tMinByTrial, tMaxByTrial] = td.alignInfoActive.getStartStopRelativeToZeroByTrial();
        end

        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrialWithPadding(td)
            % make relative to zero
            [tMinByTrial, tMaxByTrial] = td.alignInfoActive.getStartStopRelativeToZeroByTrialWithPadding();
        end

        function [tStartRelByTrial, tStopRelByTrial, tZeroRelByTrial] = getTimeStartStopZeroRelativeToTrialStartEachTrial(td)
            [tStartRelByTrial, tStopRelByTrial, tZeroRelByTrial] = td.alignInfoActive.getStartStopZeroRelativeToTrialStartByTrial();
            tStartRelByTrial(~td.valid) = NaN;
            tStopRelByTrial(~td.valid) = NaN;
            tZeroRelByTrial(~td.valid) = NaN;
        end

        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrialEachAlign(td)
            [tMinByTrial, tMaxByTrial] = deal(nan(td.nTrials, td.nAlign));

            for iA = 1:pset.nAlign
                [tMinByTrial(:, iA), tMaxByTrial(:, iA)] = ...
                    td.alignInfoSet{iA}.getStartStopRelativeToZeroByTrial();
            end

            tMinByTrial(~td.valid, :) = NaN;
            tMaxByTrial(~td.valid, :) = NaN;
        end

        function offsets = getTimeOffsetsFromZeroEachTrial(td, varargin)
            p = inputParser();
            p.addParameter('raw', false, @islogical);
            p.parse(varargin{:})

            offsets = td.alignInfoActive.getZeroByTrial();
            if ~p.Results.raw
                offsets(~td.valid) = NaN;
            end
        end

        function tf = alignIncludesFullTrial(td)
            tf = td.alignInfoActive.isFullTrial;
        end

        function offsets = getTimeOffsetsFromZeroEachTrialEachAlign(td)
            offsets = nan(td.nTrials, td.nAlign);
            for iA = 1:pset.nAlign
                offsets(:, iA) = td.alignInfoSet{iA}.getZeroByTrial();
            end
            offsets(~td.valid, :) = NaN;
        end

        function [dataRaw, timesRaw] = replaceDataWithinAlignWindow(td, dataRaw, timesRaw, dataAligned, timesAligned, varargin)
            % data/timesRaw are taken from the full and unaligned. data/timesAligned are aligned within the current alignment window
            % and should be spliced in over top of the align window in dataRaw to yield data/times (which will be unaligned)
            p = inputParser();
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('isAligned', true, @islogical); % timesAligned already relative to current zero
            p.parse(varargin{:});

            assert(iscell(timesRaw) && ismatrix(timesRaw) && size(timesRaw, 1) == td.nTrials);
            assert(iscell(dataRaw) && size(dataRaw, 1) == td.nTrials);
            assert(iscell(dataAligned) && size(dataAligned, 1) == td.nTrials);
            assert(iscell(timesAligned) && ismatrix(timesAligned) && size(timesAligned, 1) == td.nTrials);

            [~, indFirst, indLast] = td.alignInfoActive.getAlignedTimesMask(...
                timesRaw, 'includePadding', p.Results.includePadding);

            alignedTimesMask = td.alignInfoActive.getAlignedTimesMask(...
                timesAligned, 'includePadding', p.Results.includePadding, 'isAligned', p.Results.isAligned); % already zero-relative?
            if p.Results.isAligned
                offsets = td.getTimeOffsetsFromZeroEachTrial();
            else
                offsets = zeros(td.nTrials, 1);
            end
            offsets = repmat(offsets, 1, size(dataRaw, 2));

            for iT = 1:numel(dataRaw)
                if ~isnan(indFirst(iT))
                    dpre = dataRaw{iT}(1:indFirst(iT), :, :, :, :);
                    tpre = timesRaw{iT}(1:indFirst(iT), :, :, :, :);

                    dmid = dataAligned{iT}(alignedTimesMask{iT}, :, :, :, :, :);
                    tmid = timesAligned{iT}(alignedTimesMask{iT}, :, :, :, :, :) + offsets(iT);

                    dpost = dataRaw{iT}(indLast(iT)+1:end, :, :, :, :);
                    tpost = timesRaw{iT}(indLast(iT)+1:end, :, :, :, :);

                    dataRaw{iT} = cat(1, dpre, dmid, dpost);
                    timesRaw{iT} = cat(1, tpre, tmid, tpost);
                else
                    % select empty row so the result is cat1-able
                    dataRaw{iT} = dataRaw{iT}([], :, :, :, :);
                    timesRaw{iT} = timesRaw{iT}([], :, :, :, :);
                end
            end
        end

    end

    % Adjustment of adjacent trial-trial boundaries
    methods
        function [td, offsetTrialStart, offsetTrialEnd]  = adjustAdjacentTrialBoundaries(td, varargin)
            % this function assumes that any two successive trials with successive trial ids represent a continuous period
            % of time, such that the boundary between them could be adjusted:
            %  1. shift to next trial:
            %     If align start is TrialStart and align stop is Event + offset, removes data after <Event + offset> from each trial
            %     and moves it to the beginning of the next trial. If Event does not occur on a given trial, uses extendToIncludeRelTrialStart(next trial) < 0
            %     to determine where to start moving data to the beginning of the next trial.
            %  2. shift to prev trial:
            %     if align start is Event + offset and align stop is TrialEnd, removes data from the start of each trial and
            %     moves it to the end of the previous trial. If Event does not occur on a given trial, uses extendToIncludeRelTrialStart(prev trial) > 0 
            %     to determine where to start moving data to the end of the previous trial.
            %  3. If unaligned, (i.e. aligned TrialStart : TrialEnd), then extendToIncludeRelTrialStart will determine the behavior
            %     If all(extendToIncludeRelTrialStart > 0), will use shift to prev trial. if all(extendToIncludeRelTrialStart < 0), 
            %     will use shift to next trial. Otherwise an error is thrown
            %
            % copyDataAcrossBoundary is a critical parameter. If false (default), the trial boundaries will simply be shifted. If true,
            % data will be copied to the adjacent trial, so that two copies of the data will exist in the adjacent trials. This should only be done
            % at the very end and only once. 
            % 
            % if ignoreEvents it true, event markers will not be copied across trial boundaries, which can sometimes be problematic.

            p = inputParser();
            p.addParameter('showProgress', true, @islogical);
            p.addParameter('interTrialOffset', 1, @isscalar); % number of time units (typically ms) to insert between trials
            p.addParameter('trialSpliceEvent', 'TrialSplice', @TrialDataUtilities.String.isstringlike);
            p.addParameter('extendToIncludeRelTrialStart', [], @(x) isempty(x) || isvector(x)); 
            p.addParameter('copyDataAcrossBoundary', false, @isscalar);
            p.addParameter('ignoreEvents', false, @(x) islogical(x) || isstringlike(x));
            p.parse(varargin{:});
            interTrialOffset = p.Results.interTrialOffset;
            trialSpliceEvent = p.Results.trialSpliceEvent;
            copyDataAcrossBoundary = p.Results.copyDataAcrossBoundary;
            
            ignoreEvents = p.Results.ignoreEvents;
            if islogical(ignoreEvents)
                ignoreAllEvents = ignoreEvents;
                ignoreSpecificEvents = strings(0, 1);
            else
                ignoreAllEvents = false;
                ignoreSpecificEvents = string(ignoreEvents);
            end
                
            nTrials = td.nTrials;
            trialId = td.getParamRaw('trialId');
            is_adjacent_with_next = [diff(trialId) == 1; false];
            is_adjacent_with_prev = [false; diff(trialId) == 1];
            
            % gather all events up front
            ai = td.alignInfoActive;
%             align_valid = ai.computedValid;
            ev_trialStart = td.getEventRawFirst('TrialStart');
            ev_trialEnd = td.getEventRawFirst('TrialEnd');
            align_trialStart = ai.timeInfo.start;
            align_trialEnd = ai.timeInfo.stop;

            extendToIncludeRelTrialStart = p.Results.extendToIncludeRelTrialStart;
            
            respectStartStopBoundaries = true; % use extendToIncludeRelTrialStart only if the adjacent trial's start / stop is missing
            ai = td.alignInfoActive;
            if ai.isStartTrialStart
                if ai.isStopTrialEnd
                    assert(~isempty(extendToIncludeRelTrialStart) && ~all(isnan(extendToIncludeRelTrialStart)), 'Either start or stop alignment must not be TrialStart / TrialEnd, or extendToIncludeRelTrialStart must be specified');
                    % using extendToIncludeRelTrialStart to determine events, not start and stop
                    respectStartStopBoundaries = false;
                    mask_nonnan = ~isnan(extendToIncludeRelTrialStart);
                    if all(extendToIncludeRelTrialStart(mask_nonnan) < 0)
                        shiftToNext = true;
                    elseif all(extendToIncludeRelTrialStart(mask_nonnan) > 0)
                        shiftToNext = false;
                    else
                        error('All values in extendToIncludeRelTrialStart must have the same sign in order to determine which direction timepoints are shifting');
                    end
                else
                    shiftToNext = true;
                end

            elseif ai.isStopTrialEnd
                shiftToNext = false;

            else
                error('At least one of start and stop alignment must be TrialStart / TrialEnd in order to adjust trial boundaries');
            end
            
            if isempty(extendToIncludeRelTrialStart)
                extendToIncludeRelTrialStart = nan(nTrials, 1);
            end
            
            clamp = @(x, lo, hi) max(min(x, hi), lo);
            
            % shift trial start or trial end events and compute trial splice times
            ev_trialSplice = nan(nTrials, 1);
            if shiftToNext
                % shift timepoints from end of trial i to start of trial i+1 
                % shift the end event of trial i earlier to match align stop, 
                % and the start event of trial i+1 earlier to consume the extra samples at the end of trial i
                [new_trialStart, new_trialEnd] = deal(nan(nTrials, 1));
                shifted_from_prev = 0;
                for iR = 1:nTrials
                    if shifted_from_prev > 0
                        ev_trialSplice(iR) = ev_trialStart(iR);
                    end
                    new_trialStart(iR) = ev_trialStart(iR) - shifted_from_prev;
                    if is_adjacent_with_next(iR) 
                        if isnan(align_trialEnd(iR)) || ~respectStartStopBoundaries
                            % no trialend or we don't care, use the value in extendToIncludeRelTrialStart for the next trial to see where to shift samples over
                            if iR == nTrials || isnan(extendToIncludeRelTrialStart(iR+1))
                                % no value given, so we won't shift anything
                                new_trialEnd(iR) = ev_trialEnd(iR);
                            else
                                % figure out how far back in time the subsequent trial wants to proceed into (extendToIncludeRelTrialStart(iR+1) < 0)
                                % max is such that it can consume the entirety of this trial but no more. There used to be an interTrialOffset term here but no longer
                                new_trialEnd(iR) = clamp(ev_trialEnd(iR) + extendToIncludeRelTrialStart(iR+1), 0, ev_trialEnd(iR));
                            end
                        else
                            % shift everything after the current trial end
                            new_trialEnd(iR) = align_trialEnd(iR);
                        end
                        
                    else
                        % not adjacent or no end event to splice at, don't shift anything
                        new_trialEnd(iR) = ev_trialEnd(iR);
                    end
                    shifted_from_prev = ev_trialEnd(iR) - new_trialEnd(iR);
                end

            else
                % shift timepoints from beginning of trial i to the end of trial i-1
                % shift the start event of trial i later to match align start,
                % and the   end   event of trial i-1 later to consume the extra samples at the beginnning of trial i
                [new_trialStart, new_trialEnd] = deal(nan(nTrials, 1));
                shifted_from_next = 0;
                for iR = nTrials:-1:1
                    if shifted_from_next > 0
                        ev_trialSplice(iR) = ev_trialEnd(iR);
                    end
                    new_trialEnd(iR) = ev_trialEnd(iR) + shifted_from_next;
                    if is_adjacent_with_prev(iR)
                        if isnan(align_trialStart(iR)) || ~respectStartStopBoundaries
                            % no TrialStart or we don't care, use the value in extendToIncludeRelTrialStart for the previous trial to see where to shift samples over
                            if iR == 1 || isnan(extendToIncludeRelTrialStart(iR-1))
                                % no value given, so we won't shift anything
                                new_trialStart(iR) = ev_trialStart(iR);
                            else
                                % figure out how far forward in time the previous trial wants to take from the beginning of this one 
                                % (extendToIncludeRelTrialStart(iR-1 > 0))
                                % min is such that it can consume the entirety of this trial but no more
                                prev_trial_wants_beyond_end = extendToIncludeRelTrialStart(iR-1) - (ev_trialEnd(iR-1) - ev_trialStart(iR-1));
                                new_trialStart(iR) = clamp(ev_trialStart(iR) + prev_trial_wants_beyond_end, ev_trialStart(iR), ev_trialEnd(iR));
                            end
                        else
                            % shift everything before the current trial start
                            new_trialStart(iR) = align_trialStart(iR);
                        end
                    else
                        % not adjacent or no start event to splice at, don't shift anything
                        new_trialStart(iR) = ev_trialStart(iR);
                    end
                    shifted_from_next = new_trialStart(iR) - ev_trialStart(iR);
                end
            end
            
            function indLast = getReadjustedBoundariesShiftToNext(time_cell, timeScaling)
                % indLast is the last index that belongs in trial, indLast+1:end will go to next trial
                if ~iscell(time_cell)
                    time_cell = num2cell(time_cell);
                end
                nC = size(time_cell, 2);
                indLast = nan(nTrials, nC);
                for iF = 1:nTrials
                    for iC = 1:nC
                        mask = double(time_cell{iF, iC}) * timeScaling <= new_trialEnd(iF);
                        if any(mask)
                            indLast(iF, iC) = find(mask, 1, 'last');
                        else
                            indLast(iF, iC) = numel(mask);
                        end
                    end
                end
            end
            
            function indFirst = getReadjustedBoundariesShiftToPrev(time_cell, timeScaling)
                % indFirst is first valid index that belongs in trial, 1:indFirst-1 will go to previous trial
                if ~iscell(time_cell)
                    time_cell = num2cell(time_cell);
                end
                nC = size(time_cell, 2);
                indFirst = nan(nTrials, nC);
                for iF = 1:nTrials
                    for iC = 1:nC
                        mask = double(time_cell{iF, iC}) * timeScaling >= new_trialStart(iF);
                        if any(mask)
                            indFirst(iF, iC) = find(mask, 1, 'first');
                        else
                            indFirst(iF, iC) = 1;
                        end
                    end
                end
            end
            
            showProgress = p.Results.showProgress;
                
            % shift event channels
            if ~ignoreAllEvents
                channels = setdiff(td.listEventChannels(), ["TrialStart", "TrialEnd", "TimeZero"]);
                channels = setdiff(channels, ignoreSpecificEvents);
                
                if showProgress
                    prog = ProgressBar(numel(channels), 'Adjusting trial boundaries for event channels');
                end
                for c = 1:numel(channels)
                    ch = channels{c};
                    if showProgress, prog.update(c, 'Adjusting trial boundaries for event %s', ch); end
                    cd = td.channelDescriptorsByName.(ch);
                    adjust_boundary_internal(cd.dataFieldPrimary, [], false);
                end
                if showProgress, prog.finish(); end
            end

            % shift analog channels
            % organize access to analog fields based on their time field since time fields may be shared
            channels = cat(1, td.listAnalogChannels('includeNamedSubChannels', false, 'includeTransformChannels', false), ...
                                td.listAnalogChannelGroups('includeTransformChannels', false));
            timeFieldByChannel = arrayfun(@(ch) string(td.getAnalogTimeField(ch)), channels);
            [timeFields, ~, whichTimeField] = unique(timeFieldByChannel);
            if showProgress
                prog = ProgressBar(numel(channels), 'Adjusting trial boundaries for analog channels');
            end
            for c = 1:numel(timeFields)
                timeField = timeFields{c};
                associatedChannels = channels(whichTimeField==c);
                if showProgress, prog.update(c, 'Adjusting trial boundaries for %d analog channels with time field %s', numel(associatedChannels), timeField); end
                adjust_boundary_internal(timeField, associatedChannels, false);
            end
            if showProgress, prog.finish(); end

            % shift spike channels
            channels = cat(1, td.listSpikeChannels('includeArraySubChannels', false), ...
                                td.listExplicitSpikeArrays());
            if showProgress
                prog = ProgressBar(numel(channels), 'Adjusting trial boundaries for spike channels');
            end
            for c = 1:numel(channels)
                ch = channels{c};
                if showProgress, prog.update(c, 'Adjusting trial boundaries for spike %s', ch); end
                cd = td.channelDescriptorsByName.(ch);

                associated_fields = strings(0, 1);
                if cd.hasWaveforms, associated_fields(end+1) = string(cd.waveformsField); end %#ok<AGROW>
                if cd.hasBlankingRegions, associated_fields(end+1) = string(cd.blankingRegionsField); end %#ok<AGROW>
                if cd.hasSortQualityEachTrial, associated_fields(end+1) = string(cd.sortQualityEachTrialField); end %#ok<AGROW>
                isArray = isa(cd, 'SpikeArrayChannelDescriptor');
                adjust_boundary_internal(cd.dataFieldPrimary, associated_fields, isArray, cd.timeScaling);

                % update the data class to be signed
                if ~strcmp(cd.timeOriginalDataClass, "")
                    cd.timeOriginalDataClass = TrialDataUtilities.Data.getSignedDataType(cd.timeOriginalDataClass);
                    td.channelDescriptorsByName.(ch) = cd;
                end
            end
            if showProgress, prog.finish(); end

            if copyDataAcrossBoundary
                new_trialStart = min(ev_trialStart, new_trialStart);
                new_trialEnd = max(ev_trialEnd, new_trialEnd);
            end
            
            % extra outputs used to update any associated data that would be segmented into to trials to match TrialData
            offsetTrialStart = new_trialStart - ev_trialStart;
            offsetTrialEnd = new_trialEnd - ev_trialEnd;
            
            assert(~any(isnan(new_trialStart)) && ~any(isnan(new_trialEnd)));
            td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, ["TrialStart"; "TrialEnd"], [new_trialStart, new_trialEnd]);
            td = td.postDataChange(fieldnames(td.data));
            td = td.reset();
            if ~isempty(trialSpliceEvent)
                td = td.addEvent(trialSpliceEvent, ev_trialSplice, 'isAligned', false);
            end

            function adjust_boundary_internal(timeField, associatedFields, eachTrialIsCell, timeScaling)
                if nargin < 4
                    timeScaling = 1;
                end
                if eachTrialIsCell
                    time_data = cat(1, td.data.(timeField));
                else
                    time_data = {td.data.(timeField)}';
                end
                nativeClass = TrialDataUtilities.Data.getCellElementClass(time_data);
                signedClass = TrialDataUtilities.Data.getSignedDataType(nativeClass);
                toSignedClass = @(x) cast(x, signedClass);

                nA = numel(associatedFields);
                associated_data = cell(nA, 1);
                for iA = 1:nA
                    if eachTrialIsCell
                        associated_data{iA} = cat(1, td.data.(associatedFields{iA}));
                    else
                        associated_data{iA} = {td.data.(associatedFields{iA})}';
                    end
                end

                nDataCol = size(time_data(:, :), 2); % if eachTrialIsCell, data cell will have nTrial rows but have multiple columns

                if shiftToNext
                    % shift the end of this trial to the beginning of next trial
                    [time_next, time_buffer] = deal(cell(nDataCol, 1)); % end of previous trial to place at beginning of this one
                    [associated_buffers, associated_next] = deal(cell(nA, nDataCol));
                    trialEnd_previous = NaN; % need to keep track of offset from trial end

                    %[~, ~, indLastEachTrial] = ai.getAlignedTimesMask(time_data, 'raw', true, 'singleTimepointTolerance', 0, 'edgeTolerance', 0);
                    indLastEachTrial = getReadjustedBoundariesShiftToNext(time_data, timeScaling);
                    for iT = 1:td.nTrials-1
                        trialStart_this = ev_trialStart(iT);

                        for iC = 1:nDataCol
                            indLast = indLastEachTrial(iT, iC);

                            % cache end values for next trial, incorporating offset from trialEnd_previous into negative offset relative to this trialStart
                            associated_empty = cellfun(@(assoc) isempty(assoc{iT, iC}), associated_data); % sometimes time will be non-empty but the values will be empty
                            if ~any(associated_empty)
                                time_next{iC} = toSignedClass(time_data{iT, iC}(indLast+1:end, :));
                                for iA = 1:nA
                                    associated_next{iA, iC} = associated_data{iA}{iT, iC}(indLast+1:end, :, :, :, :, :, :, :, :, :);
                                end
                            end

                            if ~isempty(time_buffer{iC})
                                % prepend while slicing, including offsets
                                % we assume that something that occurred at the end of the previous trial now occurs at
                                % the start of this trial minus interTrialOffset

                                if copyDataAcrossBoundary
                                    retainTime = toSignedClass(time_data{iT, iC});
                                else
                                    retainTime = toSignedClass(time_data{iT, iC}(1:indLast, :));
                                end

%                                 time_data{iT, iC} = TrialDataUtilities.Data.catPromoteNumeric(1, time_buffer{iC} - trialEnd_previous + trialStart_this - interTrialOffset, retainTime);
                                time_data{iT, iC} = toSignedClass(cat(1, double(time_buffer{iC}) + double(-trialEnd_previous + trialStart_this - interTrialOffset)/timeScaling, double(retainTime)));

                                for iA = 1:nA
                                    if copyDataAcrossBoundary
                                        retainData = associated_data{iA}{iT, iC};
                                    else
                                        retainData = associated_data{iA}{iT, iC}(1:indLast, :, :, :, :, :, :, :, :, :);
                                    end
                                    associated_data{iA}{iT, iC} = TrialDataUtilities.Data.catPromoteNumeric(1, associated_buffers{iA, iC}, retainData);
                                end
                            elseif ~copyDataAcrossBoundary
                                % no need to prepend, just slice off end
                                time_data{iT, iC} = toSignedClass(time_data{iT, iC}(1:indLast, :));
                                for iA = 1:nA
                                    associated_data{iA}{iT, iC} = associated_data{iA}{iT, iC}(1:indLast, :, :, :, :, :, :, :, :, :);
                                end
                            end
                        end

                        % store next trial values in buffer
                        time_buffer = time_next;
                        associated_buffers = associated_next;
                        trialEnd_previous = ev_trialEnd(iT);
                    end

                else
                    % shift from the beginning of this trial to the end of the previous trial
                    % shift the end of this trial to the beginning of next trial
                    [time_prev, time_buffer] = deal(cell(nDataCol, 1)); % end of previous trial to place at beginning of this one
                    [associated_buffers, associated_prev] = deal(cell(nA, nDataCol));
                    trialStart_next = NaN; % need to keep track of offset from trial end

                    %[~, indFirstEachTrial, ~] = ai.getAlignedTimesMask(time_data, 'raw', true, 'singleTimepointTolerance', 0, 'edgeTolerance', 0);
                    indFirstEachTrial = getReadjustedBoundariesShiftToPrev(time_data, timeScaling);
                    for iT = td.nTrials:-1:1
                        trialEnd_this = ev_trialEnd(iT);

                        for iC = 1:nDataCol
                            indFirst = indFirstEachTrial(iT, iC);
                            % cache end values for next trial, incorporating offset from trialEnd_previous into negative offset relative to this trialStart
                            associated_empty = cellfun(@(assoc) isempty(assoc{iT, iC}), associated_data); % sometimes time will be non-empty but the values will be empty
                            if ~any(associated_empty)
                                time_prev{iC} = toSignedClass(time_data{iT, iC}(1:indFirst-1, :));
                                for iA = 1:nA
                                    associated_prev{iA, iC} = associated_data{iA}{iT, iC}(1:indFirst-1, :, :, :, :, :, :, :, :, :);
                                end
                            end

                            if ~isempty(time_buffer{iC})
                                % postpend while offsetting times, slice off beginning
                                % we assume that something that occurred at the start of the subsequent trial now occurs at
                                % the end of this trial plus interTrialOffset
                                if copyDataAcrossBoundary
                                    retainTime = toSignedClass(time_data{iT, iC});
                                else
                                    retainTime = toSignedClass(time_data{iT, iC}(indFirst:end, :));
                                end
                                
%                                 time_data{iT, iC} = TrialDataUtilities.Data.catPromoteNumeric(1, retainTime, time_buffer{iC} - trialStart_next + trialEnd_this + interTrialOffset);
                                time_data{iT, iC} = toSignedClass(cat(1, double(retainTime), double(time_buffer{iC}) + double(-trialStart_next + trialEnd_this + interTrialOffset)/timeScaling));
                                for iA = 1:nA
                                    if copyDataAcrossBoundary
                                        retainData = associated_data{iA}{iT, iC};
                                    else
                                        retainData = associated_data{iA}{iT, iC}(indFirst:end, :, :, :, :, :, :, :, :, :);
                                    end
                                    associated_data{iA}{iT, iC} = TrialDataUtilities.Data.catPromoteNumeric(1, retainData,  associated_buffers{iA, iC});
                                end
                                
                            elseif ~copyDataAcrossBoundary
                                % no need to postpend, just slice off beginning
                                time_data{iT, iC} = toSignedClass(time_data{iT, iC}(indFirst:end, :));
                                for iA = 1:nA
                                    associated_data{iA}{iT, iC} = associated_data{iA}{iT, iC}(indFirst:end, :, :, :, :, :, :, :, :, :);
                                end
                            end
                        end

                        % store next trial values in buffer
                        time_buffer = time_prev;
                        associated_buffers = associated_prev;
                        trialStart_next = ev_trialStart(iT);
                    end

                end

                if eachTrialIsCell % split back into nTrials x 1 cell of 1 x nDataCol cells
                    time_data = mat2cell(time_data, ones(td.nTrials, 1), nDataCol);
                    for iA = 1:nA
                        associated_data{iA} = mat2cell(associated_data{iA}, ones(td.nTrials, 1), nDataCol);
                    end
                end
                assign_data = [time_data, associated_data{:}];
                assign_fields = cat(1, string(timeField), string(associatedFields));
                td.data = TrialDataUtilities.Data.assignIntoStructArray(td.data, assign_fields, assign_data);
            end
        end
    end

    % Analog channel access
    methods
        function time = getAnalogTime(td, name, varargin)
            p = inputParser;
            p.addParameter('singleTimepointTolerance', Inf, @isscalar);
            p.parse(varargin{:});

            time = td.getAnalogTimeRaw(name);
            time = td.alignInfoActive.getAlignedTimesCell(time, false, p.Results); % no padding
        end

        % return aligned analog channel
        function [data, time] = getAnalog(td, name, varargin)
            p = inputParser;
            p.addParameter('subtractTrialBaseline', [], @(x) isempty(x) || isvector(x) || isscalar(x) || isstringlike(x) || isa(x, 'function_handle'))
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.addParameter('singleTimepointTolerance', Inf, @isscalar);
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('includeEdgeBins', false, @islogical);

            % if this is specified, the data for each trial will be sampled
            % at specific times, per-trial
            p.addParameter('sampleAtTimes', [], @(x) isvector(x) || iscell(x));
            p.addParameter('extrapolate', true, @islogical);

            % if these are specified, the data for each trial will be resampled
            p.addParameter('ensureUniformSampling', false, @islogical);
            p.addParameter('timeDelta', []);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);
            p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time vector will span time start stop even if no data is present
            p.parse(varargin{:});

            includePadding = p.Results.includePadding;

            [data, time] = getAnalog@TrialData(td, name);

            % resample if requested, don't use this when embedding in a
            % common matrix, better to figure out the common time vector
            % first
            if ~isempty(p.Results.sampleAtTimes)

                % do the sampling in the unaligned raw data to ensure good
                % resampling at the edges

%                 if p.Results.includePadding
%                     [tMin, tMax] = td.getTimeStartStopEachTrialWithPadding();
%                 else
%                     [tMin, tMax] = td.getTimeStartStopEachTrial();
%                 end

                sampleAt = p.Results.sampleAtTimes;
                if ~iscell(sampleAt)
                    sampleAt = repmat({sampleAt}, td.nTrials, 1);
                end
                % shift sampleAt so that it is unaligned
                zeroOffsets = td.getTimeOffsetsFromZeroEachTrial();
                sampleAtUnaligned = cellfun(@plus, sampleAt, num2cell(zeroOffsets), 'UniformOutput', false);

                % shift tMin and tMax as well
%                 tMin = tMin + zeroOffsets;
%                 tMax = tMax + zeroOffsets;

                [data, time] = TrialDataUtilities.Data.resampleDataCellAtSpecificTimes(data, time, sampleAtUnaligned, ...
                    'interpolateMethod', p.Results.interpolateMethod, ...
                    'extrapolate', p.Results.extrapolate, ...
                    'binAlignmentMode', p.Results.binAlignmentMode);

                % and then finally align everything to zero
                [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, includePadding, ...
                    'singleTimepointTolerance', p.Results.singleTimepointTolerance);

            elseif ~isempty(p.Results.timeDelta) || p.Results.ensureUniformSampling
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = td.getAnalogTimeDelta(name);
                end
                if p.Results.includePadding
                    [tMin, tMax] = td.getTimeStartStopEachTrialWithPadding();
                else
                    [tMin, tMax] = td.getTimeStartStopEachTrial();
                end

                % pad a bit forward or backwards depending on binning
                td = td.padForTimeBinning(timeDelta, p.Results.binAlignmentMode, p.Results.includePadding, p.Results.includeEdgeBins);

                [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, true, ...
                    'singleTimepointTolerance', p.Results.singleTimepointTolerance);

                time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, p.Results.timeReference);

                [data, time] = TrialDataUtilities.Data.resampleDataCellInTime(data, time, 'timeDelta', timeDelta, ...
                    'timeReference', p.Results.timeReference, 'binAlignmentMode', p.Results.binAlignmentMode, ...
                    'resampleMethod', p.Results.resampleMethod, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'tMinExcludingPadding', tMin, 'tMaxExcludingPadding', tMax, 'expandToTimeMinMax', p.Results.expandToTimeMinMax);
            else
                [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, includePadding, ...
                    'singleTimepointTolerance', p.Results.singleTimepointTolerance);
            end

            % subtract baseline on condition by condition
            if ~isempty(p.Results.subtractConditionBaselineAt)
                if strcmp(p.Results.subtractConditionBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(p.Results.subtractConditionBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baselineByCondition = tdBaseline.getAnalogMeanOverTimeEachTrialGroupMeans(name, 'singleTimepointTolerance', p.Results.singleTimepointTolerance);
                baselineForTrial = nanvec(tdBaseline.nTrials);
                mask = ~isnan(tdBaseline.conditionIdx);
                baselineForTrial(mask) = baselineByCondition(tdBaseline.conditionIdx(mask));
                baselineForTrial_splitTrials = mat2cell(baselineForTrial, ones(size(baselineForTrial, 1), 1), size(baselineForTrial, 2));
                data = cellfun(@minus, data, baselineForTrial_splitTrials, 'UniformOutput', false);
            end

            % subtract baseline on trial by trial basis
            if ~isempty(p.Results.subtractTrialBaselineAt)
                subtractTrialBaselineAt = p.Results.subtractTrialBaselineAt;
                if strcmp(subtractTrialBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(subtractTrialBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baseline = tdBaseline.getAnalogMeanOverTimeEachTrial(name, 'singleTimepointTolerance', p.Results.singleTimepointTolerance);
                baseline_splitTrials = mat2cell(baseline, ones(size(baselineForTrial, 1), 1), size(baselineForTrial, 2));
                data = cellfun(@minus, data, baseline_splitTrials, 'UniformOutput', false);
            end

            % subtract manual offset from each trial
            if ~isempty(p.Results.subtractTrialBaseline)
                sub = p.Results.subtractTrialBaseline;
                if isa(sub, 'function_handle')
                    % call it on each element of data
                    subValues = cellfun(sub, data);
                elseif isstringlike(sub)
                    % treat as parameter value
                    subValues = td.getParam(sub);
                elseif isscalar(sub)
                    subValues = repmat(sub, td.nTrials, 1);
                else
                    subValues = sub;
                end

                assert(numel(subValues) == td.nTrials, 'subtractTrialBaseline must be scalar of have nTrials entries');
                data = cellfun(@minus, data, num2cell(subValues), 'UniformOutput', false);
            end
        end

        function [data, time] = getAnalogEachAlign(td, name, varargin)
            % data, time are nTrials x nAlign cells
            [data, time] = deal(cell(td.nTrials, td.nAlign));
            for iA = 1:td.nAlign
                [data(:, iA), time(:, iA)] = td.useAlign(iA).getAnalog(name, varargin{:});
            end
        end

        function tvec = getAnalogCommonTimeVector(nameCell, varargin)
            p = inputParser();
            p.addParameter('timeDelta', [], @isscalar);
            p.parse(varargin{:});

            % infer timeDelta to use
            timeDelta = p.Results.timeDelta;
            if isempty(timeDelta)
                timeDelta = td.alignInfoActive.minTimeDelta;
                if isempty(timeDelta)
                    timeDelta = td.getAnalogTimeDelta(nameCell);
                    warning('timeDelta auto-computed from analog timestamps. Specify manually or call .round for consistent results');
                end
            end

            % build nTrials x nChannels cell of data/time vectors
            C = numel(nameCell);
            [timeCell, dataCell] = deal(cell(td.nTrials, C));
            for c = 1:C
                [dataCell(:, c), timeCell(:, c)] = td.getAnalog(nameCell{c}, 'includeEdgeBins', true, varargin{:});
            end

            tvec = TrialDataUtilies.Data.inferCommonTimeVectorForTimeseriesData(timeCell, dataCell, ...
                'timeDelta', timeDelta);
        end

        function [dataUnif, timeUnif, delta] = getAnalogUniformlySampled(td, name, varargin)
            [dataUnif, timeUnif, delta] = getAnalogUniformlySampled@TrialData(td, name, varargin{:});

            [dataUnif, timeUnif] = td.alignInfoActive.getAlignedTimeseries(dataUnif, timeUnif, false);
        end

        function [dCell, tCell] = getAnalogGrouped(td, name, varargin)
            [dataCell, timeCell] = td.getAnalog(name, varargin{:});
            [dCell, tCell] = td.groupElements(dataCell, timeCell);
        end

        function [dCell, tCell] = getAnalogGroupedEachAlign(td, name, varargin)
            [dataCell, timeCell] = td.getAnalogEachAlign(name, varargin{:});
            [dCell, tCell] = deal(cell(td.nConditions, td.nAlign));
            for iA = 1:td.nAlign
                [dCell(:, iA), tCell(:, iA)] = td.useAlign(iA).groupElementsFlat(dataCell(:, iA), timeCell(:, iA));
            end
        end

        function [means, tvec] = getAnalogMeanOverTimeEachTrial(td, name, varargin)
            [data, tvec] = td.getAnalog(name, varargin{:});
            %nCh = td.getAnalogChannelGroupSize(name);
            means = nanvec(td.nTrials);
            means(td.valid) = cellfun(@(x) mean(x, 'omitnan'), data(td.valid));
        end

        function [meansCell, tvec] = getAnalogMeanOverTimeEachTrialGrouped(td, name, varargin)
            [means, tvec] = td.getAnalogMeanOverTimeEachTrial(name, varargin{:});
            meansCell = td.groupElements(means);
        end

        function [meansCell, tvec] = getAnalogMeanOverTimeEachTrialGroupedRandomized(td, name, varargin)
            [means, tvec] = td.getAnalogMeanOverTimeEachTrial(name, varargin{:});
            meansCell = td.groupElementsRandomized(means);
        end

        function [means, tvec] = getAnalogMeanOverTimeEachTrialGroupMeansRandomized(td, name, varargin)
            [meansCell, tvec] = td.getAnalogMeanOverTimeEachTrialGroupedRandomized(name, varargin{:});
            means = cellfun(@(x) mean(x, 'omitnan'), meansCell);
        end
        
        function [meanMat, semMat, stdMat, nTrialsMat] = getAnalogMeanOverTimeEachTrialGroupMeans(td, name, varargin)
            % *Mat will be nConditions x T x ... matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            meansGrouped = td.getAnalogMeanOverTimeEachTrialGrouped(name, p.Unmatched);

            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan([td.nConditions, 1]));
            for iC = 1:td.nConditions
                if ~isempty(meansGrouped{iC})
                    [meanMat(iC, :, :, :, :, :, :), semMat(iC, :, :, :, :, :, :), ...
                        nTrialsMat(iC, :, :, :, :, :, :), stdMat(iC, :, :, :, :, :, :)] = ...
                        nanMeanSemMinCount(meansGrouped{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function [rms, ssqByTrial, countByTrial] = getAnalogRMSEachTrial(td, name, varargin)
            p = inputParser();
            p.addParameter('clip', [], @(x) numel(x) <= 2);
            p.addParameter('replace', NaN, @(x) numel(x) <= 2);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            data = td.getAnalog(name, p.Unmatched);
            if ~isempty(p.Results.clip)
                data = TensorUtils.clip(data, p.Results.clip, p.Results.replace);
            end

            ssqByTrial = cellfun(@(x) sum((x-mean(x, 'omitnan')).^2, 'omitnan'), data(td.valid));
            countByTrial = cellfun(@(x) nnz(~isnan(x)), data(td.valid));
            rms = sqrt(ssqByTrial ./ countByTrial);
            rms = TensorUtils.inflateMaskedTensor(rms, 1, td.valid, NaN);
            ssqByTrial = TensorUtils.inflateMaskedTensor(ssqByTrial, 1, td.valid, NaN);
            countByTrial = TensorUtils.inflateMaskedTensor(countByTrial, 1, td.valid, NaN);
        end

        function rms = getAnalogRMS(td, name, varargin)
            [~, ssqByTrial, countByTrial] = td.getAnalogRMSEachTrial(name, varargin{:});
            rms = sqrt(sum(ssqByTrial, 1, 'omitnan') ./ sum(countByTrial, 1, 'omitnan'))';
        end

        function [mat, tvec] = getAnalogAsMatrix(td, name, varargin)
            % return aligned analog channel, resampled and interpolated to
            % a uniformly spaced time vector around t=0 such that the
            % result can be embedded in a nTrials x nTime matrix. time will
            % be chosen to encapsulate the min / max timestamps across all
            % trials. Missing samples will be returned as NaN
            %
            % if name is a cellstr of multiple channels, mat will be nTrials x nTime x nChannels
            % tensor of values, but this will fail if the channels do not
            % share a common time vector to begin with. If this is not the
            % case, use getAnalogMultiAsTensor
            %
            % 'subtractTrialBaselineAt' : provide an alignment string that will
            %   optionally be used to estimate a 'baseline' value to subtract
            %   from each trace in the matrix. We'll align to this string,
            %   average the analog data over the whole interval, and subtract
            %   this off from each trial.
            %
            % 'subtractConditionBaselineAt' : same as above, but subtracts
            %   the average value for each condition, not each trial

            p = inputParser;
            p.addParameter('timeDelta', [], @(x) isscalar(x) || isempty(x));
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp

            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);

            p.addParameter('interpolateMethod', 'linear', @isstringlike); % see interp1 for details
            p.addParameter('assumeUniformSampling', false, @islogical);
            p.addParameter('minTrials', 0, @isscalar);
            p.addParameter('minTrialFraction', 0, @isscalar);
            
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time vector will span time start stop even if no data is present
            p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); % manually dictate time boundaries if specified, otherwise auto
            p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); 
            p.parse(varargin{:});

            % build nTrials cell of data/time vectors, and have getAnalog
            % do any resampling
            [dataCell, timeCell] = td.getAnalog(name, ...
                'subtractTrialBaseline', p.Results.subtractTrialBaseline, ...
                'subtractTrialBaselineAt', p.Results.subtractTrialBaselineAt, ...
                'subtractConditionBaselineAt', p.Results.subtractConditionBaselineAt, ...
                'includePadding', p.Results.includePadding, ...
                'includeEdgeBins', true, ...
                'ensureUniformSampling', true, ...
                'timeReference', 0, 'timeDelta', p.Results.timeDelta, ...
                'binAlignmentMode', p.Results.binAlignmentMode, ...
                'resampleMethod', p.Results.resampleMethod, ...
                'interpolateMethod', p.Results.interpolateMethod);

            tMin = p.Results.tMin;
            tMax = p.Results.tMax;
            if p.Results.expandToTimeMinMax
                if p.Results.includePadding
                    [tMin_, tMax_] = td.getTimeStartStopEachTrialWithPadding();
                else
                    [tMin_, tMax_] = td.getTimeStartStopEachTrial();
                end
                if isempty(tMin)
                    tMin = min(tMin_);
                end
                if isempty(tMax)
                    tMax = max(tMax_);
                end
            end

            % interpolate to common time vector
            % mat is nTrials x nTime
            [mat, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell, ...
                'assumeUniformSampling', true, ... % since getAnalog already ensured uniform sampling
                'minTrials', p.Results.minTrials, ...
                'fixNonmonotonicTimes', false, ... % already ensured by uniform sampling
                'minTrialFraction', p.Results.minTrialFraction, 'trialValid', td.valid, ...
                'tMin', tMin, 'tMax', tMax);
        end

        function [matCell, tvec] = getAnalogAsMatrixGrouped(td, name, varargin)
            [mat, tvec] = td.getAnalogAsMatrix(name, varargin{:});
            matCell = td.groupElements(mat);
        end

        function [mat, tvec, alignIdx] = getAnalogAsMatrixEachAlign(td, name, varargin)
            % similar to getAnalogAsMatrix, except each alignment will be
            % concatenated in time

            matCell = cellvec(td.nAlign);
            tvecCell = cellvec(td.nAlign);
            for iA = 1:td.nAlign
                [matCell{iA}, tvecCell{iA}] = td.useAlign(iA).getAnalogAsMatrix(name, varargin{:});
            end

            [mat, alignIdx] = TensorUtils.catWhich(2, matCell{:});
            tvec = cat(2, tvecCell{:});
        end

        function [matCell, tvec, alignIdx] = getAnalogAsMatrixGroupedEachAlign(td, name, varargin)
            [mat, tvec, alignIdx] = td.getAnalogAsMatrixEachAlign(name, varargin{:});
            matCell = td.groupElements(mat);
        end

        function [meanMat, semMat, tvec, stdMat, nTrialsMat] = getAnalogGroupMeans(td, name, varargin)
            % *Mat will be nConditions x T matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.addParameter('subtractConditionDescriptor', [], @(x) isempty(x) || isa(x, 'ConditionDescriptor'));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            [dCell, tvec] = td.getAnalogAsMatrixGrouped(name, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec)));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :), semMat(iC, :), nTrialsMat(iC, :), stdMat(iC, :)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end

            % subtract these conditions piecemeal from this one
            if ~isempty(p.Results.subtractConditionDescriptor)
                cdSub = td.expandConditionDescriptorForSubtraction(p.Results.subtractConditionDescriptor);

                % get the group means for this condition descriptor
                [sub_meanMat, sub_semMat, sub_tvec, sub_stdMat, ~] = ...
                    td.setConditionDescriptor(cdSub).getAnalogGroupMeans(name, ...
                    rmfield(p.Results, 'subtractConditionDescriptor'), p.Unmatched);

                % equalize the time vectors
                [out, tvec] = TrialDataUtilities.Data.equalizeTimeVectorsForTimeseries({...
                    meanMat, semMat, stdMat, sub_meanMat, sub_semMat, sub_stdMat}, ...
                    {tvec, tvec, tvec, sub_tvec, sub_tvec, sub_tvec}, 2);

                % do the subtraction
                meanMat = out{1} - out{4};

                % sd = sqrt(sd1^2 + sd2^2)
                % sem = sqrt(sem1^2 + sem2^2)
                semMat = sqrt(out{2}.^2 + out{5}.^2);
                stdMat = sqrt(out{3}.^2 + out{6}.^2);
            end

        end

        function [meanMat, semMat, tvec, stdMat, nTrialsMat] = getAnalogGroupMeansRandomized(td, name, varargin)
            % *Mat will be nConditions x T x nRandomized matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            % *Mat will be nConditions x T x nRandomized matrices
            [mat, tvec] = td.getAnalogAsMatrix(name, p.Unmatched);

             % now do grouping and mean computation in loop for each iRandom
            % so as to save memory
            prog = ProgressBar(td.nRandomized, 'Computing mean over randomizations');
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), td.nRandomized));

            for iR = 1:td.nRandomized
                prog.update(iR);

                % group for this randomization
                dCell = td.groupElementsRandomizedSingle(iR, mat);

                for iC = 1:td.nConditions
                    if ~isempty(dCell{iC})
                        [meanMat(iC, :, iR), semMat(iC, :, iR), nTrialsMat(iC, :, iR), stdMat(iC, :, iR)] = ...
                            nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                    end
                end
            end
            prog.finish();
        end

        function quantileMat = getAnalogGroupMeansRandomizedQuantiles(td, name, varargin)
            p = inputParser();
            p.addParameter('quantiles', [0.025 0.975], @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            meanMat = td.getAnalogGroupMeansRandomized(name, p.Unmatched);

            quantileMat = quantile(meanMat, p.Results.quantiles, 3); % nRandom along 3rd dimension for this
        end

        function [meanMat, semMat, tvec, stdMat, nTrialsMat, alignIdx] = getAnalogGroupMeansEachAlign(td, name, varargin)
            % *Mat will be nConditions x T x nChannels matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            if iscell(name)
                nC = numel(name);
            else
                nC = 1;
            end
            [dCell, tvec, alignIdx] = td.getAnalogAsTensorGroupedEachAlign(name, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nC));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), stdMat(iC, :, :)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function [meanMat, semMat, tvec, nTrialsMat, stdMat] = ...
                 getAnalogMultiGroupMeans(td, nameCell, varargin)
             % *Mat will be nConditions x T x nChannels tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.addParameter('subtractConditionDescriptor', [], @(x) isempty(x) || isa(x, 'ConditionDescriptor'));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            nameCell = cellstr(nameCell);
            assert(iscellstr(nameCell), 'Must be cellstr of channel names');
            nChannels = numel(nameCell);

            % dCell is size(conditions) with nTrials x T x nChannels inside
            [dCell, tvec] = td.getAnalogMultiAsTensorGrouped(nameCell, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nChannels));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), ...
                        stdMat(iC, :, :)] = nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end

            % subtract these conditions piecemeal from this one
            if ~isempty(p.Results.subtractConditionDescriptor)
                cdSub = td.expandConditionDescriptorForSubtraction(p.Results.subtractConditionDescriptor);

                % get the group means for this condition descriptor
                [sub_meanMat, sub_semMat, sub_tvec, ~, sub_stdMat] = ...
                    td.setConditionDescriptor(cdSub).getAnalogMultiGroupMeans(nameCell, ...
                    rmfield(p.Results, 'subtractConditionDescriptor'), p.Unmatched);

                % equalize the time vectors
                [out, tvec] = TrialDataUtilities.Data.equalizeTimeVectorsForTimeseries({...
                    meanMat, semMat, stdMat, sub_meanMat, sub_semMat, sub_stdMat}, ...
                    {tvec, tvec, tvec, sub_tvec, sub_tvec, sub_tvec}, 2);

                % do the subtraction
                meanMat = out{1} - out{4};

                % sd = sqrt(sd1^2 + sd2^2)
                % sem = sqrt(sem1^2 + sem2^2)
                semMat = sqrt(out{2}.^2 + out{5}.^2);
                stdMat = sqrt(out{3}.^2 + out{6}.^2);
            end
         end

         function [meanMat, semMat, tvec, nTrialsMat, stdMat, alignIdx] = ...
                 getAnalogMultiGroupMeansEachAlign(td, nameCell, varargin)
             % *Mat will be nConditions x T x nChannels tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            nameCell = cellstr(nameCell);
%             assert(iscellstr(nameCell), 'Must be cellstr of channel names');
            nChannels = numel(nameCell);

            % dCell is size(conditions) with nTrials x T x nChannels inside
            [dCell, tvec, alignIdx] = td.getAnalogMultiAsTensorGroupedEachAlign(nameCell, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nChannels));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), ...
                        stdMat(iC, :, :)] = nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
         end
    end

    methods % Analog first sample as a scalar, useful for aligning to single timepoint
         function [dataVec, timeVec] = getAnalogSample(td, name, varargin)
             % same as in TrialData, except issues warning if alignment has
             % more than one sample. Pulls a single time point from the t=0
             % aligned time for each trial. Guaranteed to be scalar per
             % trial, i.e. vector over trials
             p = inputParser;
             p.addParameter('subtractTrialBaseline', [], @(x) isempty(x) || isvector(x) || isstringlike(x) || isa(x, 'function_handle'))
             p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
             p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
             p.addParameter('singleTimepointTolerance', Inf, @isscalar);
             p.parse(varargin{:});
             if ~td.alignInfoActive.isStartStopEqual
                 warning('Getting analog sample when alignment is not for single timepoint');
             end
             [dataVec, timeVec] = getAnalogSample@TrialData(td, name, varargin{:});

             % subtract baseline on condition by condition
            if ~isempty(p.Results.subtractConditionBaselineAt)
                if strcmp(p.Results.subtractConditionBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(p.Results.subtractConditionBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(~td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baselineByCondition = tdBaseline.getAnalogMeanOverTimeEachTrialGroupMeans(name, 'singleTimepointTolerance', p.Results.singleTimepointTolerance);
                baselineForTrial = nanvec(tdBaseline.nTrials);
                mask = ~isnan(tdBaseline.conditionIdx);
                baselineForTrial(mask) = baselineByCondition(tdBaseline.conditionIdx(mask));
                dataVec = dataVec - baselineForTrial;
            end

            % subtract baseline on trial by trial basis
            if ~isempty(p.Results.subtractTrialBaselineAt)
                if strcmp(p.Results.subtractTrialBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(p.Results.subtractTrialBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baseline = tdBaseline.getAnalogMeanOverTimeEachTrial(name, 'singleTimepointTolerance', p.Results.singleTimepointTolerance);
                dataVec = dataVec - baseline;
            end

            % subtract manual offset from each trial
            if ~isempty(p.Results.subtractTrialBaseline)
                sub = p.Results.subtractTrialBaseline;
                if isa(sub, 'function_handle')
                    % call it on each element of data
                    subValues = arrayfun(sub, dataVec);
                elseif isstringlike(sub)
                    % treat as parameter value
                    subValues = td.getParam(sub);
                elseif isscalar(sub)
                    subValues = repmat(sub, td.nTrials, 1);
                else
                    subValues = sub;
                end

                assert(numel(subValues) == td.nTrials, 'subtractTrialBaseline must be scalar of have nTrials entries');
                dataVec = dataVec - subValues;
            end
        end

        function [dataCell, timeCell] = getAnalogSampleGrouped(td, name, varargin)
            [dataVec, timeVec] = td.getAnalogSample(name, varargin{:});
            [dataCell, timeCell] = td.groupElements(dataVec, timeVec);
        end

        function [dataCell, timeCell] = getAnalogSampleGroupedRandomized(td, name, varargin)
            [dataVec, timeVec] = td.getAnalogSample(name, varargin{:});
            [dataCell, timeCell] = td.groupElementsRandomized(dataVec, timeVec);
        end

        function [meanMat, semMat, stdMat, nTrialsMat] = getAnalogSampleGroupMeans(td, name, varargin)
            % get averaged analog sample value within each group
            % *Mat will be size(conditions) tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            dCell = td.getAnalogSampleGrouped(name, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.conditionsSize));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC), semMat(iC), nTrialsMat(iC), stdMat(iC)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function [meanMat, semMat, stdMat, nTrialsMat] = getAnalogSampleGroupMeansRandomized(td, name, varargin)
            % get averaged analog sample value within each group
            % *Mat will be size(conditions) x nRandomized tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            dCell = td.getAnalogSampleGroupedRandomized(name, p.Unmatched);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan([td.conditionsSize td.nRandomized]));
            for iC = 1:numel(dCell)
                if ~isempty(dCell{iC})
                    [meanMat(iC), semMat(iC), nTrialsMat(iC), stdMat(iC)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end
    end

    methods % Analog channel modification
        % differentiation
        function [diffData, time, diffUnits] = differentiateAnalogChannel(td, name, varargin)
            % data will be in units / second, checking .timeUnitsPerSecond
            % in order to normalize appropriately
            p = inputParser;
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('interpolateMethod', 'linear', @isstringlike);
            p.addParameter('resampleMethod', 'interp', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('smoothing', NaN, @(x) isnan(x) || (isscalar(x) && mod(x, 2) == 1));
            p.addParameter('smoothingMs', NaN, @isscalar); %
            p.addParameter('differentiationOrder', 1, @isscalar);
            p.addParameter('polynomialOrder', 2, @isscalar);
            p.addParameter('progress', false, @islogical);
            p.parse(varargin{:});
            progress = p.Results.progress;

            name = string(name);
            if numel(name) == 1 
                if td.hasAnalogChannelGroup(name)
                    % fetch as group, rather than single channel
                    [data, time, delta] = td.getAnalogChannelGroupUniformlySampled(name, ...
                        'timeDelta', p.Results.timeDelta, 'interpolateMethod', p.Results.interpolateMethod, ...
                        'resampleMethod', p.Results.resampleMethod, ...
                        'ensureUniformSampling', true);
                else
                    [data, time, delta] = td.getAnalogUniformlySampled(name, ...
                    'timeDelta', p.Results.timeDelta, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'resampleMethod', p.Results.resampleMethod, ...
                    'ensureUniformSampling', true);
                end
            elseif numel(name) > 1
                [data, time, delta] = td.getAnalogMultiCommonTimeUniformlySampled(name, ...
                    'timeDelta', p.Results.timeDelta, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'resampleMethod', p.Results.resampleMethod, ...
                    'ensureUniformSampling', true);
            end

            diffData = cellvec(td.nTrials);

            deltaMs = delta / td.timeUnitsPerMs;
            minSmoothingSamples = p.Results.polynomialOrder + 1;
            minSmoothingMs = minSmoothingSamples * deltaMs;

            % figure out smoothing in samples if specified in ms
            if ~isnan(p.Results.smoothing)
                smoothing = p.Results.smoothing;
            elseif ~isnan(p.Results.smoothingMs)
                smoothingMs = p.Results.smoothingMs;
                % convert from ms to samples using ms * samples/ms
                smoothing = round(smoothingMs / deltaMs);
                if mod(smoothing, 2) == 0
                    smoothing = smoothing + 1;
                end
            else
                % default to minimum allowed
                smoothing = minSmoothingSamples;
            end

            % check smoothing is large enough for sampling rate
            assert(smoothing >= minSmoothingSamples, 'Smoothing width is too small, must be at least %d samples or %g ms', ...
                minSmoothingSamples, minSmoothingMs);

%             w = -1 / (delta / td.timeUnitsPerSecond) ^ p.Results.order;
            if iscell(name)
                nameStr = 'channels';
            else
                nameStr = name;
            end
            if progress, prog = ProgressBar(td.nTrials, 'Smoothing/Differentiating %s', nameStr); end
            for iT = 1:td.nTrials
                if isempty(data{iT})
                    continue;
                end

                if any(sum(~isnan(data{iT}), 1) < smoothing)
                    % too few samples
                    diffData{iT} = nan(size(data{iT}));
                else
                    diffData{iT} = TrialDataUtilities.Data.savitzkyGolayFilt( ...
                        data{iT}, 'polynomialOrder', p.Results.polynomialOrder, 'differentiationOrder', p.Results.differentiationOrder, ...
                        'frameSize', smoothing, 'samplingIntervalMs', delta / td.timeUnitsPerMs);
                end
                if progress, prog.update(iT); end
            end
            if progress, prog.finish(); end

            diffUnits = sprintf('%s/sec', td.getChannelUnitsPrimary(name));
        end

        function td = addDifferentiatedAnalogChannel(td, name, diffName, varargin)
            % drops alignment and grouping before differentiating since
            % this is typically what the user wants

            td.warnIfNoArgOut(nargout);
            [diffData, time, diffUnits] = td.reset().differentiateAnalogChannel(name, varargin{:});
            td = td.addAnalog(diffName, diffData, time, 'units', diffUnits, 'isAligned', false);
        end

        function td = addAnalogChannelModifiedFromExisting(td, oldName, newName, data)
            % if we do a reset here we'll mess up the valid array the user
            % was expecting
%             td = td.reset();
            td.warnIfNoArgOut(nargout);
            td = td.copyChannel(oldName, newName);
            td = td.setAnalog(newName, data, 'updateValidOnly', true, 'clearForInvalid', true);
        end

        function td = addScaledAnalogChannel(td, name, scaledName, multiplyBy, scaledUnits)
            td.warnIfNoArgOut(nargout);

            [data, time] = td.getAnalogRaw(name);
            data = cellfun(@(d) d*multiplyBy, data, 'UniformOutput', false);
            td = td.addAnalog(scaledName, data, time, 'units', scaledUnits, 'isAligned', false);
        end

        function td = addResampledAnalogChannel(td, name, newName, timeDelta, varargin)
            td.warnIfNoArgOut(nargout);

            [data, time] = td.reset().getAnalog(name, 'timeDelta', timeDelta, varargin{:});

            td = td.copyChannel(name, newName);
            td = td.setAnalog(newName, data, time, 'isAligned', false, 'updateValidOnly', true, 'clearForInvalid', true);
        end

        function [data, time] = getAnalogFiltered(td, name, B, A, varargin)
            p = inputParser;
            p.addParameter('filtfilt', false, @islogical);
            p.addParameter('subtractFirstSample', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if td.hasAnalogChannel(name)
                [data, time] = td.getAnalog(name, p.Unmatched);
            elseif td.hasAnalogChannelGroup(name)
                [data, time] = td.getAnalogChannelGroup(name, p.Unmatched);
            end
            data = TrialDataUtilities.Data.filterIgnoreLeadingTrailingNaNs(B, A, data, ...
                'filtfilt', p.Results.filtfilt, 'subtractFirstSample', p.Results.subtractFirstSample);
        end

        function td = addFilteredAnalogChannel(td, name, filtName, B, A, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            data = td.getAnalogFiltered(name, B, A, varargin{:});
            td = td.addAnalogChannelModifiedFromExisting(name, filtName, data);
        end

        function [data, time] = getAnalogLowPassFiltered(td, name, order, cornerHz, varargin)
            Fs = td.getAnalogSamplingRateHz(name);
            cornerNormalized = cornerHz / (Fs/2);
            [B, A] = butter(order, cornerNormalized, 'low');
            [data, time] = td.getAnalogFiltered(name, B, A, varargin{:});
        end

        function td = addLowPassFilteredAnalogChannel(td, name, filtName, order, cornerHz, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            data = td.getAnalogLowPassFiltered(name, order, cornerHz, varargin{:});
            td = td.addAnalogChannelModifiedFromExisting(name, filtName, data);
        end

        function [data, time] = getAnalogHighPassFiltered(td, name, order, cornerHz, varargin)
            Fs = td.getAnalogSamplingRateHz(name);
            cornerNormalized = cornerHz / (Fs/2);
            [B, A] = butter(order, cornerNormalized, 'high');
            [data, time] = td.getAnalogFiltered(name, B, A, varargin{:});
        end

        function td = addHighPassFilteredAnalogChannel(td, name, filtName, order, cornerHz, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            data = td.getAnalogHighPassFiltered(name, order, cornerHz, varargin{:});
            td = td.addAnalogChannelModifiedFromExisting(name, filtName, data);
        end

        function [data, time] = getAnalogBandPassFiltered(td, name, order, cornerHz, varargin)
            Fs = td.getAnalogSamplingRateHz(name);
            cornerNormalized = cornerHz / (Fs/2);
            [B, A] = butter(order, cornerNormalized);
            [data, time] = td.getAnalogFiltered(name, B, A, varargin{:});
        end

        function td = addBandPassFilteredAnalogChannel(td, name, filtName, order, cornerHz, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            data = td.getAnalogBandPassFiltered(name, order, cornerHz, varargin{:});
            td = td.addAnalogChannelModifiedFromExisting(name, filtName, data);
        end

        function [data, time] = getAnalogBandStopFiltered(td, name, order, cornerHz, varargin)
            Fs = td.getAnalogSamplingRateHz(name);
            cornerNormalized = cornerHz / (Fs/2);
            [B, A] = butter(order, cornerNormalized, 'stop');
            [data, time] = td.getAnalogFiltered(name, B, A, varargin{:});
        end

        function td = addBandStopFilteredAnalogChannel(td, name, filtName, order, cornerHz, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            data = td.getAnalogBandStopFiltered(name, order, cornerHz, varargin{:});
            td = td.addAnalogChannelModifiedFromExisting(name, filtName, data);
        end

        function td = setAnalogWithinAlignWindow(td, name, values, varargin)
            % replace the analog data within the curr ent align window with
            % kthe data in values
            p = inputParser();
            p.addOptional('times', [], @(x) iscell(x) || isnumeric(x));
            p.addParameter('preserveContinuity', false, @islogical);
            p.addParameter('clearOverlappingTimesOutsideAlignWindow', false, @islogical);
            p.addParameter('keepScaling', false, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            td.assertHasChannel(name);

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

            times = p.Results.times;
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

                % check that times have same length as data
                assert(numel(times) == numel(values), 'Times and values must have same number of trials');
                nTimes = cellfun(@numel, times);
                nValues = cellfun(@numel, values);
                assert(all(nTimes == nValues), 'Mismatch between number of times and values. If the number of times has changed be sure to specify times parameter');

                updateTimes = true;
            else
                % here we're only talking about the aligned region anyway,
                % so if nothing is specified, those are the times we're
                % using anyway
                updateTimes = false;
                times = td.getAnalogTime(name);

            end

            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);

            cd = td.channelDescriptorsByName.(name);

            % for shared column channels where the scaling or time vectors
            % change, it needs to be separated from the shared column form
            if isa(cd, 'AnalogChannelDescriptor') && cd.isColumnOfSharedMatrix && ((cd.hasScaling && ~p.Results.keepScaling) || updateTimes)
                error('Setting analog channel while ''keepScaling'' is false or when specifying new sample times requires this channel to be separated from its analog channel group. Use separateAnalogChannelFromGroup if you want to do this. Or use setAnalogChannelGroupWithinAlignWindow to set all channels in the group at once.');
%                 td = td.separateAnalogChannelFromGroup(name, false);
            end

            % data being passed in is now in original units
            % so change scaling factors
            if ~p.Results.keepScaling
                td = td.convertAnalogChannelToNoScaling(name);
            end

            cd = td.channelDescriptorsByName.(name); %#ok<NASGU>

            % figure out where to splice in the data into the full timeseries
            [fullData, fullTimes] = td.getAnalogRaw(name);
            [~, timesMask] = td.alignInfoActive.getAlignedTimesCell(fullTimes, false, 'singleTimepointTolerance', Inf); % no padding

            % splice each trial's data in
            prog = ProgressBar(td.nTrials, 'Splicing in analog data');
            for iT = 1:td.nTrials
                prog.update(iT);
                mask = timesMask{iT}; % indicates where in fullTimes that currentTimes lives
                first = find(mask, 1);
                last = find(mask, 1, 'last');

                if updateTimes
                    % check whether the times we're writing don't exceed the
                    % window we're aligned to

                    preTimes = fullTimes{iT}(1:first-1);
                    preData = fullData{iT}(1:first-1);
                    postTimes = fullTimes{iT}(last+1:end);
                    postData = fullData{iT}(last+1:end);

                    % check time overlap with times outside the window
                    maxTimePre = preTimes(end);
                    minTimePost = postTimes(1);
                    if p.Results.clearOverlappingTimesOutsideAlignWindow
                        % remove overlap from the old data
                        maxTimeNew = max(times{iT}, [], 'omitnan');
                        minTimeNew = min(times{iT}, [], 'omitnan');

                        maskRemove = falsevec(numel(times{iT}));
                        maskRemove(1:first-1) = fullTimes{iT}(1:first-1) >= minTimeNew;
                        maskRemove(last+1:end) = fullTimes{iT}(1:first-1) <= maxTimeNew;

                        preTimes = preTimes(~maskRemove(1:first-1));
                        preData = preData(~maskRemove(1:first-1));
                        postTimes = postTimes(~maskRemove(last+1:end));
                        postData = postData(~maskRemove(last+1:end));
                    else
                        % issue an error if overlap
                        newTimeOutside = times{iT} >= minTimePost | times{iT} <= maxTimePre;
                        if any(newTimeOutside)
                            error('Trial %d has time values that overlap with times outside of the alignment window. Set parameter ''clearOverlappingTimesOutsideAlignWindow'' true to clear these automatically');
                        end
                    end

                else
                    preData = fullData{iT}(1:first-1);
                    postData = fullData{iT}(last+1:end);
                end

                if p.Results.preserveContinuity
                    % offset values{iT} to match last value of preData
                    % insert the data
                    values{iT} = values{iT} - (values{iT}(1) - preData(end));

                    % offset postData to match last value of values{iT}
                    postData = postData - (postData(1) - values{iT}(1));
                end

                fullData{iT} = cat(1, preData, values{iT}, postData);
                if updateTimes
                    fullTimes{iT} = cat(1, preTimes, times{iT}, postTimes);
                end
            end
            prog.finish();

            % setChannelData will call repairData which will update
            % memoryDataClassByField{1} to reflect the type of values
            if updateTimes
                td = td.setAnalog(name, fullData, fullTimes, 'isAligned', false, 'keepScaling', p.Results.keepScaling);
            else
                td = td.setAnalog(name, fullData, 'timesMatchFullTrialRaw', true, ...
                    'isAligned', false, 'keepScaling', p.Results.keepScaling);
            end
        end

        function td = trimAnalogChannelToCurrentAlign(td, names)
            td.warnIfNoArgOut(nargout);

            [startTimes, stopTimes] = td.getTimeStartStopEachTrial();

            names = TrialDataUtilities.Data.wrapCell(names);
            timeFields = unique(cellfun(@(name) td.getAnalogTimeField(name), names, 'UniformOutput', false));

            for i = 1:numel(timeFields)
                td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeFields{i}, startTimes, stopTimes);
            end
        end
    end

    methods % Analog channel group methods
        % return aligned analog channel
        function [data, time] = getAnalogChannelGroup(td, groupName, varargin)
            p = inputParser();
            p.addParameter('singleTimepointTolerance', Inf, @isscalar);
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('includeEdgeBins', false, @islogical);
            p.addParameter('subtractTrialBaseline', [], @(x) isempty(x) || isvector(x) || isstringlike(x) || isa(x, 'function_handle'))
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);

            % if these are specified, the data for each trial will be resampled
            p.addParameter('ensureUniformSampling', false, @islogical);
            p.addParameter('timeDelta', []);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);

            p.addParameter('sliceSelectSubChannels', {}, @(x) true);
            p.addParameter('slice', {}, @(x) true); % subscript args to slice the data from each sample
            p.addParameter('averageOverSlice', false, @islogical); % average within each slice

            p.addParameter('linearCombinationWeights', [], @(x) true); % alternatively, take a weighted combination over samples in the slice, size should be [size of analog channel, number of weighted combinations]

			% these apply to the weighted combination
			p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            p.addParameter('keepNaNIfAllNaNs', false, @islogical); % when replaceNaNWithZero is true, keep the result as NaN if every entry being combined is NaN
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical); % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis this enables nanmean like computations

            p.addParameter('applyScaling', true, @islogical);
            
            p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); % manually dictate time boundaries if specified, otherwise auto
            p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time vector will span time start stop even if no data is present
            p.parse(varargin{:});

            [data, time] = getAnalogChannelGroup@TrialData(td, groupName, ...
                'slice', p.Results.slice, 'sliceSelectSubChannels',  p.Results.sliceSelectSubChannels, 'averageOverSlice', p.Results.averageOverSlice, ...
                'linearCombinationWeights', p.Results.linearCombinationWeights, ...
                'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                'keepNaNIfAllNaNs', p.Results.keepNaNIfAllNaNs, ...
                'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                'applyScaling', p.Results.applyScaling);

            includePadding = p.Results.includePadding;

            % resample if requested, don't use this when embedding in a
            % common matrix, better to figure out the common time vector
            % first
            if ~isempty(p.Results.timeDelta) || p.Results.ensureUniformSampling
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = td.getAnalogChannelGroupTimeDelta(groupName);
                end
                if p.Results.includePadding
                    [tMin, tMax] = td.getTimeStartStopEachTrialWithPadding();
                else
                    [tMin, tMax] = td.getTimeStartStopEachTrial(); 
                end
                % override with auto
                if ~isempty(p.Results.tMin)
                    tMin = p.Results.tMin;
                end
                if ~isempty(p.Results.tMax)
                    tMax = p.Results.tMax;
                end

                % pad a bit forward or backwards depending on binning
                td = td.padForTimeBinning(timeDelta, p.Results.binAlignmentMode, p.Results.includePadding, p.Results.includeEdgeBins);

                [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, true, ...
                    'singleTimepointTolerance', p.Results.singleTimepointTolerance);
                [data, time] = TrialDataUtilities.Data.resampleDataCellInTime(data, time, 'timeDelta', timeDelta, ...
                    'timeReference', p.Results.timeReference, 'binAlignmentMode', p.Results.binAlignmentMode, ...
                    'resampleMethod', p.Results.resampleMethod, 'interpolateMethod', p.Results.interpolateMethod,  ...
                    'tMinExcludingPadding', tMin, 'tMaxExcludingPadding', tMax, ...
                    'expandToTimeMinMax', p.Results.expandToTimeMinMax);
            else
                [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, includePadding, ...
                    'singleTimepointTolerance', p.Results.singleTimepointTolerance);
            end

            % subtract baseline on condition by condition
            if ~isempty(p.Results.subtractConditionBaselineAt)
                if strcmp(p.Results.subtractConditionBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(p.Results.subtractConditionBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baselineByCondition = tdBaseline.getAnalogChannelGroupMeanOverTimeEachTrialGroupMeans(groupName, 'singleTimepointTolerance', p.Results.singleTimepointTolerance, ...
                    'slice', p.Results.slice, ...
                    'linearCombinationWeights', p.Results.linearCombinationWeights, ...
                    'applyScaling', p.Results.applyScaling);
                sz = size(baselineByCondition);
                sz(1) = tdBaseline.nTrials;
                baselineForTrial = nan(sz);
                mask = tdBaseline.valid;
                baselineForTrial(mask, :, :, :, :, :, :, :, :) = baselineByCondition(tdBaseline.conditionIdx(mask), :, :, :, :, :, :, :, :);
                baselineForTrial_splitTrials = mat2cell(baselineForTrial, ones(size(baselineForTrial, 1), 1), size(baselineForTrial, 2));
                data = cellfun(@minus, data, baselineForTrial_splitTrials, 'UniformOutput', false);
            end

            % subtract baseline on trial by trial basis
            if ~isempty(p.Results.subtractTrialBaselineAt)
                if strcmp(p.Results.subtractTrialBaselineAt, '*')
                    tdBaseline = td;
                else
                    tdBaseline = td.align(p.Results.subtractTrialBaselineAt);
                    tdBaseline = tdBaseline.setManualValidTo(td.valid); % this shouldn't matter since the samples will be nans anyway, but just in case
                end
                baseline = tdBaseline.getAnalogChannelGroupMeanOverTimeEachTrial(groupName, 'singleTimepointTolerance', p.Results.singleTimepointTolerance, ...
                    'slice', p.Results.slice, ...
                    'linearCombinationWeights', p.Results.linearCombinationWeights, ...
                    'applyScaling', p.Results.applyScaling);
                data = cellfun(@minus, data, TensorUtils.splitAlongDimension(baseline, 1), 'UniformOutput', false);
            end

            % subtract manual offset from each trial
            if ~isempty(p.Results.subtractTrialBaseline)
                sub = p.Results.subtractTrialBaseline;
                if isa(sub, 'function_handle')
                    % call it on each element of data
                    subValues = cellfun(sub, data);
                elseif ischar(sub) || isStringScalar(sub)
                    % treat as parameter value
                    subValues = td.getParam(sub);
                elseif isscalar(sub)
                    subValues = repmat(sub, td.nTrials, 1);
                elseif iscellstr(sub) || isstring(sub)
                    subValues = td.getParamMultiAsCell(sub);
                    subValues = cell2mat(subValues);
                else
                    subValues = sub;
                end

                assert(size(subValues, 1) == td.nTrials, 'subtractTrialBaseline must be scalar of have nTrials entries');
                data = cellfun(@(d, s) d - s, data, num2cell(subValues, 2), 'UniformOutput', false);
            end
        end

        function [dataCell, timeCell] = getAnalogChannelGroupGrouped(td, groupName, varargin)
            % dataCell will be size(td.conditions)
            % contents will be nTrials x T x nChannels
            tdValid = td.selectValidTrials();
            [data, time] = tdValid.getAnalogChannelGroup(groupName, varargin{:});
            [dataCell, timeCell] = tdValid.groupElements(data, time);
        end

        function time = getAnalogChannelGroupTime(td, groupName)
            time = getAnalogChannelGroupTime@TrialData(td, groupName);
            time = td.alignInfoActive.getAlignedTimesCell(time, false, 'singleTimepointTolerance', Inf); % no padding
        end
        
        function [rms, ssqByTrial, countByTrial] = getAnalogChannelGroupRMSEachTrial(td, name, varargin)
            p = inputParser();
            p.addParameter('clip', [], @(x) numel(x) <= 2);
            p.addParameter('replace', NaN, @(x) numel(x) <= 2);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            data = td.getAnalogChannelGroup(name, p.Unmatched);
            if ~isempty(p.Results.clip)
                data = TensorUtils.clip(data, p.Results.clip, p.Results.replace);
            end

            temp = cellfun(@(x) sum((x-mean(x, 1, 'omitnan')).^2, 1, 'omitnan'), data(td.valid), 'UniformOutput', false);
            ssqByTrial = cat(1, temp{:}); % nTrials x nChannels
            temp = cellfun(@(x) sum(~isnan(x), 1, 'omitnan'), data(td.valid), 'UniformOutput', false);
            countByTrial = cat(1, temp{:}); % nTrials x nChannels
            rms = sqrt(ssqByTrial ./ countByTrial);
            rms = TensorUtils.inflateMaskedTensor(rms, 1, td.valid, NaN);
            ssqByTrial = TensorUtils.inflateMaskedTensor(ssqByTrial, 1, td.valid, NaN);
            countByTrial = TensorUtils.inflateMaskedTensor(countByTrial, 1, td.valid, NaN);
        end

        function meanVal = getAnalogChannelGroupGlobalMeanOverTime(td, name, varargin)
            data = td.getAnalogChannelGroup(name, varargin{:});
            sums = cellfun(@(x) sum(x, 1, 'omitnan'), data, 'UniformOutput', false);
            totals = cellfun(@(x) sum(~isnan(x), 1), data, 'UniformOutput', false);

            meanVal = sum(cat(1, sums{:}), 1) ./ sum(cat(1, totals{:}), 1);
        end

        function rms = getAnalogChannelGroupRMS(td, name, varargin)
            [~, ssqByTrial, countByTrial] = td.getAnalogChannelGroupRMSEachTrial(name, varargin{:});
            rms = sqrt(sum(ssqByTrial, 1, 'omitnan') ./ sum(countByTrial, 1, 'omitnan'))';
        end

        function [dataCell, timeCell] = getAnalogMulti(td, name, varargin)
            % [data, time] = getAnalogMulti(td, chNames, varargin)
            % data and time are cell(nTrials, nChannels)
            % this overrides a method in TrialData
            p = inputParser();
            p.addParameter('raw', false, @islogical);
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            name = string(name);
            nChannels = numel(name);

            % apply subtractTrialBaseline, if its a cell, feed individually for each channel
            if ~isempty(p.Results.subtractTrialBaseline)
                subBase = p.Results.subtractTrialBaseline;
                if ~iscell(subBase)
                    subBase = repmat({subBase}, nChannels, 1);
                else
                    assert(numel(subBase) == nChannels, 'Size of subtractTrialBaseline cell must match nChannels');
                end
            else
                subBase = cellvec(nChannels);
            end

            % build nTrials x nChannels cell of data/time vectors
            C = numel(name);
            [dataCell, timeCell] = deal(cell(td.nTrials, C));
%             prog = ProgressBar(C, 'Fetching analog channels');
            for c = 1:C
%                 prog.update(c);
                if p.Results.raw
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalogRaw(name{c}, p.Unmatched);
                else
                    [dataCell(:, c), timeCell(:, c)] = td.getAnalog(name{c}, 'subtractTrialBaseline', subBase{c}, p.Unmatched);
                end
            end
%             prog.finish();
        end

        function [dataCell, timeCell] = getAnalogMultiCommonTime(td, names, varargin)
            % dataCell is cell(nTrials, 1) containing nTime x nChannels mat
            % timeCell is cell(nTrials, 1) containing nTime x 1 vectors
            % All data vectors will be interpolated to a
            % common time vector independently on each trial. Use
            % getAnalogMultiAsMatrix to register to a common time vector
            % across all trials.

            p = inputParser;
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('interpolateMethod', 'linear', @isstringlike); % see interp1 for details
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
                colIdx = td.getAnalogChannelColumnIdxInGroup(names);
%                 dataCell = cellvec(td.nTrials);
                [dataCell, timeCell] = td.getAnalogChannelGroup(groupName, 'timeDelta', p.Results.timeDelta, 'slice', colIdx, ...
                    'includeEdgeBins', true, 'timeReference', p.Results.timeReference, 'interpolateMethod', p.Results.interpolateMethod, ...
                    'binAlignmentMode', p.Results.binAlignmentMode, 'resampleMethod', p.Results.resampleMethod, p.Unmatched);

                % then go and grab the correct columns
%                 colIdx = td.getAnalogChannelColumnIdxInGroup(names);
%                 prog = ProgressBar(td.nTrials, 'Slicing columns from analog channel group');
%                 for iT = 1:td.nTrials
%                     prog.update(iT)
%                     dataCell{iT} = matCell{iT}(:, colIdx);
%                 end
%                 prog.finish();

            else
%                 sameTime = td.checkAnalogChannelsShareTimeField(names);

                % pick common sampling rate up front
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    % use common sampling rate and upsample
                    timeDelta = td.getAnalogTimeDelta(names); % will choose smallest sampling interval
                    if isempty(timeDelta)
                        timeDelta = td.alignInfoActive.minTimeDelta;
                    end
                end

                [dataCellRaw, timeCellRaw] = td.getAnalogMulti(names, 'timeDelta', timeDelta, ...
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

        function [dataUnif, timeUnif, delta] = getAnalogMultiCommonTimeUniformlySampled(td, name, varargin)
            [dataUnif, timeUnif] = td.getAnalogMultiCommonTime(name, varargin{:}, 'ensureUniformSampling', true);
            delta = TrialData.computeDeltaFromTimes(timeUnif);
        end

        function [dataByCondition, timeByCondition] = getAnalogMultiCommonTimeGrouped(td, names, varargin)
            [data, time] = td.getAnalogMultiCommonTime(names, varargin{:});
            [dataByCondition, timeByCondition] = td.groupElements(data, time);
        end

        function [dataTensor, tvec] = getAnalogMultiAsTensor(td, names, varargin)
            % data is nTrials x nTime x nChannels
            % tvec is nTime x 1 time vector
            %
            % parameters:
            %   assumeUniformScaling [false]: % assumes that all trials use consistent sampling
            %     relative to alignment 0, setting this true avoids an interpolation step
            %
            p = inputParser;

            % for resampling
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);
            p.addParameter('assumeUniformSampling', false, @islogical);

            p.addParameter('includePadding', false, @islogical);
            p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time vector will span time start stop even if no data is present
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            assert(isstring(names) || iscell(names));

            % pick common sampling rate up front
            timeDelta = p.Results.timeDelta;
            if isempty(timeDelta)
                % use common sampling rate and upsample
                timeDelta = td.getAnalogTimeDelta(names); % will choose smallest sampling interval
                if isempty(timeDelta)
                    timeDelta = td.alignInfoActive.minTimeDelta;
                end
            end

            % build nTrials x nTime x nChannels cell of data/time vectors
            [dataCell, timeCell] = td.getAnalogMulti(names, ...
                'includeEdgeBins', true, ...
                'timeDelta', timeDelta, ...
                'timeReference', p.Results.timeReference, 'interpolateMethod', p.Results.interpolateMethod, ...
                'binAlignmentMode', p.Results.binAlignmentMode, 'resampleMethod', p.Results.resampleMethod, ...
                'includePadding', p.Results.includePadding, ...
                'expandToTimeMinMax', p.Results.expandToTimeMinMax, p.Unmatched);
            
            tMin = p.Results.tMin;
            tMax = p.Results.tMax;
            if p.Results.expandToTimeMinMax
                if p.Results.includePadding
                    [tMin_, tMax_] = td.getTimeStartStopEachTrialWithPadding();
                else
                    [tMin_, tMax_] = td.getTimeStartStopEachTrial();
                end
                if isempty(tMin)
                    tMin = tMin_;
                end
                if isempty(tMax)
                    tMax = tMax_;
                end

                extra_args = {'tMin', min(tMin), 'tMax', max(tMax)};
            else
                extra_args = {};
            end
            
            % interpolate to common time vector, we can assume uniform
            % sampling since they've already been interpolated to timeDelta
            [dataTensor, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell, ...
                'assumeUniformSampling', true, extra_args{:});
        end

        function [dataTensor, tvec] = getAnalogChannelGroupAsTensor(td, groupName, varargin)
            p = inputParser;
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);

            p.addParameter('minTrials', 0, @isscalar);
            p.addParameter('minTrialFraction', 0, @isscalar);
            
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time vector will span time start stop even if no data is present
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [dataCell, timeCell] = td.getAnalogChannelGroup(groupName, ...
                'ensureUniformSampling', true, ...
                'timeDelta', p.Results.timeDelta, ...
                'includeEdgeBins', true, ...
                'binAlignmentMode', p.Results.binAlignmentMode, ...
                'resampleMethod', p.Results.resampleMethod, ...
                'interpolateMethod', p.Results.interpolateMethod, ...
                'includePadding', p.Results.includePadding, p.Unmatched);

            tMin = p.Results.tMin;
            tMax = p.Results.tMax;
            if p.Results.expandToTimeMinMax
                if p.Results.includePadding
                    [tMin_, tMax_] = td.getTimeStartStopEachTrialWithPadding();
                else
                    [tMin_, tMax_] = td.getTimeStartStopEachTrial();
                end
                if isempty(tMin)
                    tMin = min(tMin_);
                end
                if isempty(tMax)
                    tMax = max(tMax_);
                end
            end
            % interpolate to common time vector
            [dataTensor, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell, ...
                'assumeUniformSampling', true, ...
                'minTrials', p.Results.minTrials, ...
                'minTrialFraction', p.Results.minTrialFraction, 'trialValid', td.valid, ...
                'tMin', min(tMin), 'tMax', max(tMax));
        end

        function [dataCell, tvec] = getAnalogChannelGroupAsTensorGrouped(td, nameCell, varargin)
            % dataCell will be size(td.conditions)
            % contents will be nTrials x T x nChannels
            tdValid = td.selectValidTrials(); % saves memory to slice trials first to avoid nan-fills on large tensors
            [data, tvec] = tdValid.getAnalogChannelGroupAsTensor(nameCell, varargin{:});
            dataCell = tdValid.groupElements(data);
        end

        function [meanMat, semMat, tvec, stdMat, nTrialsMat] = getAnalogChannelGroupGroupMeans(td, groupName, varargin)
            % *Mat will be nConditions x T x ... matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            [dCell, tvec] = td.getAnalogChannelGroupAsTensorGrouped(groupName, p.Unmatched);
            ne = find(~cellfun(@isempty, dCell), 1);
            sampleSize = size(dCell{ne});
            sampleSize = sampleSize(3:end);

            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan([td.nConditions, numel(tvec), sampleSize]));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :, :, :, :, :), semMat(iC, :, :, :, :, :, :), ...
                        nTrialsMat(iC, :, :, :, :, :, :), stdMat(iC, :, :, :, :, :, :)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end
        
        function [snr, range, noise]  = getAnalogChannelGroupSNR(td, groupName, varargin)
            % take range = max - min (or e.g. 0.0.25, 0.975 quantile) and noise = max (or 0.95 quantile) s.e.m.
            p = inputParser();
            p.addParameter('noiseQuantile', 0.95, @isscalar);
            p.addParameter('rangeQuantile', 0.95, @(x) isscalar(x) || numel(x) == 2); % either width of centered quantile band (e.g. 0.95 --> 0.025 to 0.975) or [low high]
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            nq = p.Results.noiseQuantile;
            rq = p.Results.rangeQuantile;
            if isscalar(rq)
                rq = [(1-rq)/2 1-(1-rq)/2];
            end

            % options same as getAnalogChannelGroupGroupMeans
            [meanMat, semMat] = td.getAnalogChannelGroupGroupMeans(groupName, p.Unmatched);

            % single trial estimated average values don't count
            meanMat(semMat == 0) = NaN;
            semMat(semMat == 0) = NaN;

            noise = squeeze(TensorUtils.quantileMultiDim(semMat, nq, [1 2]));
            range = squeeze(diff(TensorUtils.quantileMultiDim(meanMat, rq, [1 2]), 1, 1));
            snr = range ./ noise;
        end

        function [snr, range, noise]  = getAnalogChannelGroupSNROverConditions(td, groupName, varargin)
            % take range = max - min over conditions at each time and noise = max s.e.m. at each time
            % then take quantile of these snr values over time
            p = inputParser();
            p.addParameter('quantile', 0.95, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            q = p.Results.quantile;
            
            % options same as getAnalogChannelGroupGroupMeans
            [meanMat, semMat] = td.getAnalogChannelGroupGroupMeans(groupName, p.Unmatched);
            range = max(meanMat,[],1, 'omitnan') - min(meanMat,[],1, 'omitnan');
            noise = max(semMat, [], 1, 'omitnan');
            snr = range ./ noise;
            snr = squeeze(quantile(snr, q, 2));
        end
        
        function [snr_c, range_c, noise_c] = getAnalogChannelGroupSNREachCondition(td, groupName, varargin)
            p = inputParser();
            p.addParameter('noiseQuantile', 0.95, @isscalar);
            p.addParameter('rangeQuantile', 0.95, @(x) isscalar(x) || numel(x) == 2); % either width of centered quantile band (e.g. 0.95 --> 0.025 to 0.975) or [low high]
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            nq = p.Results.noiseQuantile;
            rq = p.Results.rangeQuantile;
            if isscalar(rq)
                rq = [(1-rq)/2 1-(1-rq)/2];
            end

            % options same as getAnalogChannelGroupGroupMeans
            [meanMat, semMat] = td.getAnalogChannelGroupGroupMeans(groupName, p.Unmatched);

            % single trial estimated average values don't count
            meanMat(semMat == 0) = NaN;
            semMat(semMat == 0) = NaN;
  
            noise_c = squeeze(TensorUtils.quantileMultiDim(semMat, nq, 2));
            range_c = squeeze(diff(TensorUtils.quantileMultiDim(meanMat, rq, 2), 1, 2));
            snr_c = range_c ./ noise_c;
        end
            
        function tbl = getAnalogChannelGroupSNRTable(td, groupName, varargin)
            names = td.listAnalogChannelsInGroupByColumn(groupName);
            nCols = numel(names);
            index = (1:nCols)';

            [snr, range, noise] = td.getAnalogChannelGroupSNR(groupName, varargin{:});

            s = struct('col_index', num2cell(index), 'snr', num2cell(snr), 'range', num2cell(range), 'noise', num2cell(noise));
            tbl = struct2table(s, 'RowNames', string(names), 'AsArray', true);
        end

        function [mat, tvec, alignIdx] = getAnalogMultiAsTensorEachAlign(td, names, varargin)
            % similar to getAnalogMultiAsMatrix, except each alignment will be
            % concatenated in time

            matCell = cellvec(td.nAlign);
            tvecCell = cellvec(td.nAlign);
            for iA = 1:td.nAlign
                [matCell{iA}, tvecCell{iA}] = td.useAlign(iA).getAnalogMultiAsTensor(names, varargin{:});
            end

            [mat, alignIdx] = TensorUtils.catWhich(2, matCell{:});
            tvec = cat(1, tvecCell{:});
        end

        function [dataCell, tvec] = getAnalogMultiAsTensorGrouped(td, nameCell, varargin)
            % dataCell will be size(td.conditions)
            % contents will be nTrials x T x nChannels
            [data, tvec] = td.getAnalogMultiAsTensor(nameCell, varargin{:});
            dataCell = td.groupElements(data);
        end

        function [dCell, tvec, alignIdx] = getAnalogMultiAsTensorGroupedEachAlign(td, nameCell, varargin)
            [mat, tvec, alignIdx] = td.getAnalogMultiAsTensorEachAlign(nameCell, varargin{:});
            dCell = td.groupElements(mat);
        end

%         function [dataUnif, timeUnif, delta] = getAnalogChannelGroupUniformlySampled(td, name, varargin)
%             [dataUnif, timeUnif, delta] = getAnalogChannelGroupUniformlySampled@TrialData(td, name, varargin{:});
% 
%             THIS IS WRONG, getAnalogChannelGroupUniformlySampled@TrialData will returned aligned data because it calls through to a TDCA overriden method
%             %[dataUnif, timeUnif] = td.alignInfoActive.getAlignedTimeseries(dataUnif, timeUnif, false);
%         end

        function means = getAnalogChannelGroupMeanOverTimeEachTrial(td, name, varargin)
            [data, ~] = td.getAnalogChannelGroup(name, varargin{:});
            meansValid = cellfun(@(x) mean(x, 1, 'omitnan'), data(td.valid), 'UniformOutput', false);
            means = TensorUtils.inflateMaskedTensor(cat(1, meansValid{:}), 1, td.valid, NaN);
        end

        function meansCell = getAnalogChannelGroupMeanOverTimeEachTrialGrouped(td, name, varargin)
            means = td.getAnalogChannelGroupMeanOverTimeEachTrial(name, varargin{:});
            meansCell = td.groupElements(means);
        end
        
%         function means = getAnalogChanneMeanOverTimeGroupMeans(td, name, varargin)
%             meansCell = td.getAnalogChannelGroupMeanOverTimeEachTrialGrouped(name, varargin{:});
%             temp = cellfun(@(x) mean(x, 1, 'omitnan'), meansCell, 'UniformOutput', false);
%             means = cat(1, temp{:});
%         end

        function [meanMat, semMat, stdMat, nTrialsMat] = getAnalogChannelGroupMeanOverTimeEachTrialGroupMeans(td, name, varargin)
            % *Mat will be nConditions x T x ... matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            meansGrouped = td.getAnalogChannelGroupMeanOverTimeEachTrialGrouped(name, p.Unmatched);

            ne = find(~cellfun(@isempty, meansGrouped), 1);
            sampleSize = size(meansGrouped{ne});
            sampleSize = sampleSize(2:end);

            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan([td.nConditions, sampleSize]));
            for iC = 1:td.nConditions
                if ~isempty(meansGrouped{iC})
                    [meanMat(iC, :, :, :, :, :, :), semMat(iC, :, :, :, :, :, :), ...
                        nTrialsMat(iC, :, :, :, :, :, :), stdMat(iC, :, :, :, :, :, :)] = ...
                        nanMeanSemMinCount(meansGrouped{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function td = setAnalogChannelGroupWithinAlignWindow(td, groupName, values, varargin)
            % replace the analog data within the current align window with
            % the data in values.
            % values is nTrials x time x channels
            % (nTrials can be nTrialsValid if parameter
            % ''dataSpansValidTrialsOnly'' is true
            p = inputParser();
            p.addOptional('times', [], @(x) ~isstringlike(x) && (iscell(x) || ismatrix(x)));
            p.addParameter('preserveContinuity', false, @islogical);
            p.addParameter('clearOverlappingTimesOutsideAlignWindow', false, @islogical);
            p.addParameter('keepScaling', false, @islogical);
            p.addParameter('removeNaNSamples', true, @islogical);
            p.addParameter('dataSpansValidTrialsOnly', false, @islogical); % values and times have size(., 1) == nTrialsValid and should be inflated
            p.parse(varargin{:});
            times = makecol(p.Results.times);

            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelGroup(groupName);

            chList = td.listAnalogChannelsInGroup(groupName);
            nCh = numel(chList);

            if p.Results.dataSpansValidTrialsOnly
                sizeDim1Expected = td.nTrialsValid;
            else
                sizeDim1Expected = td.nTrials;
            end

            % check the values and convert to nTrials cellvec
            if isnumeric(values) %&& isnumeric(values)
                % values must be nTrials x nTimes x nChannels
                assert(size(values, 1) == sizeDim1Expected, 'Values as matrix must be nTrials (or nTrialsValid if ''dataSpansValidTrialsOnly'' set) along dimension 1');

                % convert to cellvec
                tensor = values;
                values = cellvec(sizeDim1Expected);
                for r = 1:sizeDim1Expected
                    values{r} = TensorUtils.squeezeDims(tensor(r, :, :), 1); % TODO
                end
                clear tensor;

            elseif iscell(values)
                assert(numel(values) == sizeDim1Expected, 'Values as cell must have length nTrials (or nTrialsValid if ''dataSpansValidTrialsOnly'' set)');
                nCol = cellfun(@(x) size(x, 2), values);
                assert(all(nCol == nCh), 'All elements of values must have same number of columns as channels (%d)', nCh);
                values = makecol(values);
            else
                error('Values must be numeric tensor or cell array');
            end

            if ~isempty(times)
                if iscell(times)
                    assert(numel(times) == sizeDim1Expected, 'numel(times) must match nTrials (or nTrialsValid if ''dataSpansValidTrialsOnly'' set)');
                elseif ismatrix(times)
                    if isvector(times)
                        times = repmat(makerow(times), sizeDim1Expected, 1);
                    end
                    assert(size(times,1) == sizeDim1Expected, 'size(times, 1) must match nTrials(or nTrialsValid if ''dataSpansValidTrialsOnly'' set)');
                    times = mat2cell(times', size(times, 2), onesvec(sizeDim1Expected))';
                else
                    error('Times must be numeric matrix or cell array');
                end
                updateTimes = true;
            else
                % here we're only talking about the aligned region anyway,
                % so if nothing is specified, those are the times we're
                % using anyway
                updateTimes = false;

                % pass along the current times since the data is coming in with the
                % existing alignment
                times = td.getAnalogChannelGroupTime(groupName);

                if p.Results.dataSpansValidTrialsOnly
                    times = times(td.valid);
                end
            end

            % check that times have same length as data
            assert(numel(times) == numel(values), 'Times and values must have same number of trials');
            nTimes = cellfun(@numel, times);
            nValues = cellfun(@(x) size(x, 1), values);
            assert(all(nTimes == nValues), 'Mismatch between number of times and values. If the number of times has changed be sure to specify times parameter');

            % now slice off NaN samples of each trial's data
            % this is useful in general but also specifically allows this
            % function to accept data returned from
            % getAnalogChannelGroupAsTensor
            if p.Results.removeNaNSamples
                for iT = 1:sizeDim1Expected
                    keep = any(~isnan(values{iT}), 2);
                    if ~all(keep), updateTimes = true; end
                    values{iT} = values{iT}(keep, :);
                    times{iT} = times{iT}(keep, :);
                end
            end

            % now that values and times are cells, we can inflate to
            % nTrials without wasting memory
            if p.Results.dataSpansValidTrialsOnly
                values = TensorUtils.inflateMaskedTensor(values, 1, td.valid);
                times = TensorUtils.inflateMaskedTensor(times, 1, td.valid);
            end

            % add the zero offset to the time vector for each trial
            % this is mostly for TDCA, so that alignments info is
            % preserved
            offsets = td.getTimeOffsetsFromZeroEachTrial();
            times = cellfun(@plus, times, num2cell(offsets), 'UniformOutput', false);

            % if we're updating the times field, it's possible for other
            % fields outside this group to reference our time field, so we
            % request that copies be made where necessary. Since we're not
            % going through setChannelData, we need to handle this
            % ourselves
            if updateTimes
                % time field is at index 2
                td = td.copyRenameSharedChannelFields(chList, 2);
            end

            timeField = td.channelDescriptorsByName.(groupName).timeField;

            if ~p.Results.keepScaling
                % data being passed in is now in original units
                % so change scaling factors of the channel and the data
                td = td.convertAnalogChannelGroupToNoScaling(groupName);
                % and convert back to memory anyway in case the data class
                % has changed, although this shouldn't do any scaling
                values = td.channelDescriptorsByName.(groupName).convertAccessDataCellToMemory(1, values);
            else
                % take new data back into scaled values to match the
                % existing
%                 assert(td.checkAnalogChannelGroupHasUniformScaling(groupName), 'Analog channel group must have uniform scaling');
                cd = td.channelDescriptorsByName.(groupName);
                values = cd.convertAccessDataCellToMemory(1, values);
            end

            fullTimes = td.getAnalogChannelGroupTimeRaw(groupName);
            [~, timesMask] = td.alignInfoActive.getAlignedTimesCell(fullTimes, false, 'singleTimepointTolerance', Inf); % no padding

            % splice each trial's data in
            prog = ProgressBar(td.nTrials, 'Splicing in analog data');
            for iT = 1:td.nTrials
                if ~td.valid(iT), continue; end
                prog.update(iT);
                mask = timesMask{iT}; % indicates where in fullTimes that currentTimes lives
                first = find(mask, 1);
                last = find(mask, 1, 'last');

                fullDataThis = td.data(iT).(groupName);

                if updateTimes
                    % check whether the times we're writing don't exceed the
                    % window we're aligned to

                    preTimes = fullTimes{iT}(1:first-1);
                    preData = fullDataThis(1:first-1, :);
                    postTimes = fullTimes{iT}(last+1:end);
                    postData = fullDataThis(last+1:end, :);

                    % check time overlap with times outside the window
                    maxTimePre = preTimes(end);
                    minTimePost = postTimes(1);
                    if p.Results.clearOverlappingTimesOutsideAlignWindow
                        % remove overlap from the old data
                        maxTimeNew = max(times{iT}, [], 'omitnan');
                        minTimeNew = min(times{iT}, [], 'omitnan');

                        maskRemove = falsevec(numel(times{iT}));
                        maskRemove(1:first-1) = fullTimes{iT}(1:first-1) >= minTimeNew;
                        maskRemove(last+1:end) = fullTimes{iT}(1:first-1) <= maxTimeNew;

                        preTimes = preTimes(~maskRemove(1:first-1));
                        preData = preData(~maskRemove(1:first-1), :);
                        postTimes = postTimes(~maskRemove(last+1:end));
                        postData = postData(~maskRemove(last+1:end), :);
                    else
                        % issue an error if overlap
                        newTimeOutside = times{iT} >= minTimePost | times{iT} <= maxTimePre;
                        if any(newTimeOutside)
                            error('Trial %d has time values that overlap with times outside of the alignment window. Set parameter ''clearOverlappingTimesOutsideAlignWindow'' true to clear these automatically');
                        end
                    end
                else
                    preData = fullDataThis(1:first-1, :);
                    postData = fullDataThis(last+1:end, :);
                end

                if p.Results.preserveContinuity
                    % offset values{iT} to match last value of preData
                    % insert the data
                    values{iT} = bsxfun(@minus, values{iT}, values{iT}(1, :) - preData(end, :));

                    % offset postData to match last value of values{iT}
                    postData = bsxfun(@minus, postData, postData(1, :) - values{iT}(end, :));
                end

                fullDataThis = cat(1, preData, values{iT}, postData);
                td.data(iT).(groupName) = fullDataThis;

                if updateTimes
                    fullTimesThis = cat(1, preTimes, times{iT}, postTimes);
                    td.data(iT).(timeField) = fullTimesThis;
                end

                assert(size(fullDataThis, 1) == numel(td.data(iT).(timeField)), 'Mismatch introduced between data and time vectors');
            end
            prog.finish();

            if updateTimes
                td = td.postDataChange({groupName, timeField});
            else
                td = td.postDataChange({groupName});
            end
        end

        function td = addAnalogChannelGroupModifiedFromExisting(td, oldName, newName, data, varargin)
            % if we do a reset here we'll mess up the valid array the user
            % was expecting
%             td = td.reset();
            td.warnIfNoArgOut(nargout);
            td = td.copyChannel(oldName, newName);
            td = td.setAnalogChannelGroup(newName, data, varargin{:});
        end

        function td = addDifferentiatedAnalogChannelGroup(td, groupName, diffName, varargin)
            % drops alignment and grouping before differentiating since
            % this is typically what the user wants

            p = inputParser();
            p.addParameter('channelNames', {}, @iscellstr);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            [diffData, diffTime, diffUnits] = td.reset().differentiateAnalogChannel(groupName, p.Unmatched);
            td = td.addAnalogChannelGroupModifiedFromExisting(groupName, diffName, diffData, ...
                'times', diffTime, 'units', diffUnits, 'isAligned', false, p.Results);
        end

        function td = addAnalogChannelGroupViaSelection(td, groupName, newName, varargin)
            p = inputParser();
            p.addParameter('chNames', {}, @iscellstr);
            p.addParameter('slice', [], @(x) true); % this is used to index specifically into each sample
            p.addParameter('linearCombinationWeights', [], @(x) true); % alternatively, take a weighted combination over samples in the slice, size should be [size of analog channel, number of weighted combinations]
            p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            p.addParameter('keepNaNIfAllNaNs', false, @islogical); % when replaceNaNWithZero is true, keep the result as NaN if every entry being combined is NaN
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical); % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis this enables nanmean like computations
            p.addParameter('averageOverSlice', false, @islogical); % average within each slice

            p.parse(varargin{:});
            td.warnIfNoArgOut(nargout);
            data = td.getAnalogChannelGroup(groupName, 'slice', p.Results.slice, ...
                'linearCombinationWeights', p.Results.linearCombinationWeights, ...
                'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                'keepNaNIfAllNaNs', p.Results.keepNaNIfAllNaNs, ...
                'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                'averageOverSlice', p.Results.averageOverSlice);

            td = td.addAnalogChannelGroupModifiedFromExisting(groupName, newName, data, 'channelNames', p.Results.chNames);

%             if ~isempty(p.Results.chNames)
%                 td = td.setAnalogChannelGroupSubChannelNames(newName, p.Results.chNames);
%             end
        end

        function td = transformAnalogChannelGroupInPlace(td, groupName, transformFn, varargin)
            % dataMat = tranformFn(dataMat, tvec, trialInd)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical); % hand transformFn scaled data?
            p.addParameter('preserveContinuity', false, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            td.assertHasAnalogChannelGroup(groupName);

            %assert(td.checkAnalogChannelGroupHasUniformScaling(groupName), 'Analog channel group must be uniformly scaled');

            cd = td.getChannelDescriptor(groupName);
            timeField = cd.timeField;

            fullTrial = td.alignIncludesFullTrial();
            if fullTrial
                % we do this so that the timestamps match up with what
                % we're about to replace them with
                td = td.trimAnalogChannelGroupToTrialStartEnd(groupName);
            end

            % figure out where to splice in the data into the full timeseries
            fullTimes = td.getAnalogChannelGroupTimeRaw(groupName);
            [~, timesMask] = td.alignInfoActive.getAlignedTimesCell(fullTimes, false, 'singleTimepointTolerance', Inf); % no padding

            prog = ProgressBar(td.nTrials, 'Tranforming %s in place', groupName);
            for iT = 1:td.nTrials
                if ~td.valid(iT), continue; end
                prog.update(iT);

                dataFull = td.data(iT).(groupName);
                timeFull = td.data(iT).(timeField);

                dataMat = dataFull(timesMask{iT}, :);
                tvec = timeFull(timesMask{iT});

                if p.Results.applyScaling
                    dataMat = cd.convertDataSingleOnAccess(1, dataMat);
                end

                % call the supplied transform function
                dataMatReplace = transformFn(dataMat, tvec, iT);

                assert(isequal(size(dataMat), size(dataMatReplace)), ...
                    'TransformFn must return data the same size as input data matrix, trial %d', iT);

                % undo scaling if we applied it in the first place
                if p.Results.applyScaling
                    dataMatReplace = cd.convertAccessDataSingleToMemory(dataMatReplace);
                else
                    assert(isequal(class(dataMatReplace), class(dataMat)), 'TransformFn must preserve class of received data');
                end

                if fullTrial
                    % just set the new value
                    oldVal = td.data(iT).(groupName);
                    assert(isequal(size(oldVal), size(dataMatReplace)), ...
                    'Internal error with data trimming on trial %d', iT);

                    td.data(iT).(groupName) = dataMatReplace;
                else
                    % splice in the new value
                    if p.Results.preserveContinuity
                        first = find(timesMask{iT}, 1, 'first');
                        last = find(timesMask{iT}, 1, 'last');
                        dataMatReplace = bsxfun(@minus, dataMatReplace, dataMatReplace(1, :) - dataFull(first-1, :));

                        dataFull(timesMask{iT}, :) = dataMatReplace;

                        dataFull(last+1, :) = bsxfun(@minus, dataFull(last+1, :), dataFull(last+1, :) - dataFull(last, :));
                    else
                        dataFull(timesMask{iT}, :) = dataMatReplace;
                    end
                    td.data(iT).(groupName) = dataFull;
                end
            end
            prog.finish();
        end

        function td = trimAnalogChannelGroupToCurrentAlign(td, groupNames)
            td.warnIfNoArgOut(nargout);

            [startTimes, stopTimes] = td.getTimeStartStopEachTrial();

            groupNames = TrialDataUtilities.Data.wrapCell(groupNames);
            timeFields = unique(cellfun(@(group) td.getAnalogChannelGroupTimeField(group), groupNames, 'UniformOutput', false));

            for i = 1:numel(timeFields)
                td = td.trimAnalogChannelTimeFieldAndReferencingChannelsRaw(timeFields{i}, startTimes, stopTimes);
            end
        end

        function td = resampleAnalogChannelGroup(td, groupName, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.reset();
            [data, time] = td.getAnalogChannelGroup(groupName, varargin{:});

            td = td.setAnalogChannelGroup(groupName, data, time);
        end
    end

    % Event channel access
    methods
        % return aligned event times
        function timesCell = getEvent(td, name)
            timesCell = getEvent@TrialData(td, name);
            timesCell = td.alignInfoActive.getAlignedTimesCell(timesCell, false, 'singleTimepointTolerance', 0, 'edgeTolerance', 0);
        end

        function dCell = getEventGrouped(td, name)
            dCell = td.groupElements(td.getEvent(name));
        end

        function firstCell = getEventFirstGrouped(td, name)
            firstCell = td.groupElements(td.getEventFirst(name));
        end

        function lastCell = getEventLastGrouped(td, name)
            lastCell = td.groupElements(td.getEventLast(name));
        end

        function td = addEvent(td, eventName, varargin)
            td.warnIfNoArgOut(nargout);
            td = addEvent@TrialData(td, eventName, varargin{:});

            % force .eventData and .eventCounts to be updated
            td = td.updateEventData(eventName);
        end

        function td = addEventByThresholdingAnalogChannel(td, eventChName, analogChName, level, varargin)
            td.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addParameter('lockoutPeriod', 0, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [data, time] = td.getAnalog(analogChName);
            thresholdTimes = TrialDataUtilities.Data.findThresholdCrossings(data, time, level, p.Results.lockoutPeriod);

            td = td.addEvent(eventChName, thresholdTimes, p.Unmatched);
        end

        function td = trimEventToCurrentAlign(td, names)
            % Timepoints that lie outside of TrialStart and TrialStop will
            % never be accessible via getTimes since they will be filtered
            % out by the AlignInfo

            td.warnIfNoArgOut(nargout);
            % default is TrialStart and TrialEnd, so just pass it along
            [startTimes, stopTimes] = td.getTimeStartStopEachTrial();
            td = td.trimSpikeChannel(names, startTimes, stopTimes);
        end
    end

    % Param channel access
    methods
        function dCell = getParamGrouped(td, name)
            % get parameter values grouped by condition
            % dCell will be size(conditions) cell tensor with value arrays
            % within
            dCell = td.groupElements(td.getParam(name));
        end

        function dCell = getParamGroupedRandomized(td, name)
            dCell = td.groupElementsRandomized(td.getParam(name));
        end

        function values = getParamUniqueGrouped(td, name)
            vCell = td.getParamGrouped(name);
            values = cellfun(@getUnique, vCell, 'UniformOutput', false);

            function values = getUnique(vals)
                if ~iscell(vals) && ~isstring(vals)
                    vals = removenan(vals);
                end
                values = unique(vals);
            end
        end

        function [meanMat, semMat, stdMat, nTrialsMat] = getParamGroupMeans(td, name, varargin)
            % get averaged parameter value within each group
            % *Mat will be size(conditions) tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            dCell = td.getParamGrouped(name);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.conditionsSize));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC), semMat(iC), nTrialsMat(iC), stdMat(iC)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function [meanMat, semMat, stdMat, nTrialsMat] = getParamGroupMeansRandomized(td, name, varargin)
            % get averaged parameter value within each group
            % *Mat will be size(conditions) x nRandomized tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;

            dCell = td.getParamGroupedRandomized(name);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan([td.conditionsSize td.nRandomized]));
            for iC = 1:numel(dCell)
                if ~isempty(dCell{iC})
                    [meanMat(iC), semMat(iC), nTrialsMat(iC), stdMat(iC)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials, p.Results.minTrialFraction);
                end
            end
        end

        function prepareAxesForChannels(td, nameX, nameY, varargin)
            p = inputParser();
            p.addParameter('axh', gca, @isscalar);
            p.parse(varargin{:});
            axh = p.Results.axh;

            cdX = td.channelDescriptorsByName.(nameX);
            cdY = td.channelDescriptorsByName.(nameY);
            TrialDataUtilities.Plotting.setupAxisForChannel(cdX, ...
                'which', 'x', 'style', 'tickBridge', 'axh', axh);
            TrialDataUtilities.Plotting.setupAxisForChannel(cdY, ...
                'which', 'y', 'style', 'tickBridge', 'axh', axh);
        end

        function plotParamScatter(td, nameX, nameY, varargin)
            dataX = td.getParamGrouped(nameX);
            dataY = td.getParamGrouped(nameY);
            cdX = td.getParamChannelDescriptor(nameX);
            cdY = td.getParamChannelDescriptor(nameY);
            td.plotProvidedGroupedScatterData(dataX, dataY, ...
                'axisInfoX', cdX, 'axisInfoY', cdY, varargin{:});
        end

        function plotParamScatter3(td, nameX, nameY, nameZ, varargin)
            dataX = td.getParamGrouped(nameX);
            dataY = td.getParamGrouped(nameY);
            dataZ = td.getParamGrouped(nameZ);
            cdX = td.getParamChannelDescriptor(nameX);
            cdY = td.getParamChannelDescriptor(nameY);
            cdZ = td.getParamChannelDescriptor(nameZ);
            td.plotProvidedGroupedScatterData(dataX, dataY, dataZ, ...
                'axisInfoX', cdX, 'axisInfoY', cdY, 'axisInfoZ', cdZ, varargin{:});
        end


        function plotParamScatterMatrixPCA(td, name, varargin)
            p = inputParser();
            p.addParameter('dims', 3, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            dims = p.Results.dims;
            data = td.getParam(name);
            [~, dataPCA] = pca(data(td.valid, :), 'NumComponents', dims);
            dataPCA = TensorUtils.inflateMaskedTensor(dataPCA, 1, td.valid);

            dataX = td.groupElements(dataPCA(:, 1));
            dataY = td.groupElements(dataPCA(:, 2));

            if dims == 2
                td.plotProvidedGroupedScatterData(dataX, dataY, ...
                    'axisInfoX', 'PC 1', 'axisInfoY', 'PC 2', p.Unmatched);
            else
                dataZ = td.groupElements(dataPCA(:, 3));
                td.plotProvidedGroupedScatterData(dataX, dataY, dataZ, ...
                    'axisInfoX', 'PC 1', 'axisInfoY', 'PC 2', 'axisInfoZ', 'PC 3', p.Unmatched);
            end
        end


        function plotProvidedGroupedScatterData(td, varargin)
            % common utility function for drawing data points grouped by
            % condition
            %
            % time is either:
            %     vector, for nAlign == 1
            %     nAlign x 1 cell of ime vectors for each align
            %
            %  data, dataErrorY is either
            %     nConditions x T x D matrix
            %     nAlign x 1 cell of nConditions x T x D matrices
            %
            p = inputParser();
            p.addRequired('dataX', @(x) isnumeric(x) || iscell(x));
            p.addRequired('dataY', @(x) isnumeric(x) || iscell(x));
            p.addOptional('dataZ', [], @(x) isnumeric(x) || iscell(x));

            p.addParameter('axisInfoX', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoZ', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));

            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @isstringlike);
            p.addParameter('axisStyleY', 'tickBridge', @isstringlike);

            p.addParameter('conditionIdx', 1:td.nConditions, @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('alpha', 0.7, @isscalar);
            p.addParameter('edgeAlpha', [], @isscalar);
            p.addParameter('markerSize', 6, @isscalar);
            p.addParameter('useThreeVector', true, @islogical);

            p.addParameter('xJitter', 0, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

%             alpha = p.Results.alpha;
%             if isempty(p.Results.edgeAlpha)
%                 edgeAlpha = alpha;
%             else
%                 edgeAlpha = p.Results.edgeAlpha;
%             end

            axh = td.getRequestedPlotAxis(p.Unmatched);

            if isempty(p.Results.dataZ)
                dims = 2;
            else
                dims = 3;
            end

            conditionIdx = p.Results.conditionIdx; %#ok<*PROPLC>
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            app = td.conditionAppearances(conditionIdx);

            dataX = p.Results.dataX;
            dataY = p.Results.dataY;
            if isempty(dataX) || isempty(dataY)
                error('Must provide dataX and dataY');
            end

            if dims == 2
                if ~ismember('scaleBars', p.UsingDefaults) && p.Results.scaleBars
                    axisStyleX = 'scaleBar';
                    axisStyleY = 'scaleBar';
                else
                    axisStyleX = p.Results.axisStyleX;
                    axisStyleY = p.Results.axisStyleY;
                end

                h = TrialDataUtilities.Plotting.allocateGraphicsHandleVector(nConditionsUsed);
                for iC = 1:nConditionsUsed
                    idxC = conditionIdx(iC);
                    args = app(idxC).getMarkerPlotArgs(false);
                    if isempty(dataX{iC}) || isempty(dataY{iC})
                        continue;
                    end

                    if p.Results.xJitter > 0
                        X = dataX{iC} + (rand(size(dataX{iC}))-0.5) * p.Results.xJitter;
                    else
                        X = dataX{iC};
                    end

    %                 h(iC) = plot(dataX{iC}, dataY{iC}, 'o', 'MarkerSize', p.Results.markerSize, ...
    %                     args{:}, p.Results.plotOptions{:});
                    h(iC) = scatter(X, dataY{iC}, p.Results.markerSize, ...
                        args{:}, 'MarkerFaceAlpha', p.Results.alpha, p.Results.plotOptions{:});
                    hold on;

                    TrialDataUtilities.Plotting.showInLegend(h(iC), td.conditionNamesShort{idxC});
                end

                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                        'which', 'x', 'style', axisStyleX);
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                        'which', 'y', 'style', axisStyleY);
                end

                hold(axh, 'off');
                box(axh, 'off');
                axis(axh, 'tight');

                AutoAxis.updateIfInstalled(axh);

            else
                % 3D plot
                dataZ = p.Results.dataZ;
                if isempty(dataX) || isempty(dataY)
                    error('Must provide dataX and dataY');
                end

                h = TrialDataUtilities.Plotting.allocateGraphicsHandleVector(nConditionsUsed);
                for iC = 1:nConditionsUsed
                    idxC = conditionIdx(iC);
                    args = app(idxC).getMarkerPlotArgs(false);
                    if isempty(dataX{iC}) || isempty(dataY{iC}) || isempty(dataZ{iC})
                        continue;
                    end

                    if p.Results.xJitter > 0
                        X = dataX{iC} + (rand(size(dataX{iC}))-0.5) * p.Results.xJitter;
                    else
                        X = dataX{iC};
                    end

                    h(iC) = scatter3(X, dataY{iC}, dataZ{iC}, p.Results.markerSize, ...
                        args{:}, 'MarkerFaceAlpha', p.Results.alpha, p.Results.plotOptions{:});
                    hold on;

                    TrialDataUtilities.Plotting.showInLegend(h(iC), td.conditionNamesShort{idxC});
                end

                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                        'which', 'x', 'style', 'label');
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                        'which', 'y', 'style', 'label');
                end
                if ~isempty(p.Results.axisInfoZ)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoZ, ...
                        'which', 'z', 'style', 'label');
                end

                hold(axh, 'off');
                box(axh, 'off');
                axis(axh, 'tight');
                axis(axh, 'vis3d');

                AutoAxis.updateIfInstalled(axh);
            end

%             for iC = 1:numel(h)
%                 % not sure why this needs to go last
%                 SaveFigure.setMarkerOpacity(h(iC), alpha, edgeAlpha);
%             end
        end

        function plotParamVsParamMeanLinePlot(td, nameX, nameY, varargin)
            % for each bin or grouping of nameX, plot the mean of paramY
            % as a line plot with vertical error bars conveying either sem
            % or std (according to 'errorType' parameter). One line plot
            % for each group defined in td (before adding nameX to the
            % list).

            p = inputParser();
            p.addParameter('errorType', 'sem', @isstringlike);
            p.addParameter('LineWidth', 1, @isscalar);
            p.addParameter('connectMissing', false, @islogical);
            p.addParameter('xJitter', 0, @isscalar);
            p.KeepUnmatched = true;
            p.CaseSensitive = false;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);
            hold(axh, 'on');

            % cache before messing with axis shapes
            conditionNames = td.conditionNamesShort;

            % add as the last axis
            td = td.addAxis(nameX);
            td = td.setAxisValueListAutoAll(nameX); % important to match attribute values in the axis list

            % put nameX as axis 1, combine other axes as axis 2
            if td.nAxes > 1
                td = td.reshapeAxes(td.nAxes, 1:(td.nAxes-1));
            end

            xValues = td.conditionInfo.getAttributeValueList(nameX);
            if td.conditionInfo.getIsAttributeBinned(nameX)
                % bins are in cell array of 1 x 2 vectors
                xValues = mean(cat(1, xValues{:}), 2);
            else
                % leave as is
                if numel(xValues) > 50
                    warning('Grouping parameter %s has >50 unique values. You may wish to bin or group this value before plotting means against it', nameX);
                end
            end

            [meanY, semY, stdY] = td.getParamGroupMeans(nameY);
            if strcmp(p.Results.errorType, 'sem')
                errorY = semY;
            elseif strcmp(p.Results.errorType, 'std')
                errorY = stdY;
            else
                error('Unknown errorType %s', p.Results.errorType);
            end

            app = td.conditionAppearances;

            nCondOther = td.conditionsSize(2);
            h = gobjects(nCondOther, 1);
            for iC = 1:nCondOther
                if p.Results.connectMissing
                    mask = ~isnan(makecol(xValues)) & ~isnan(meanY(:, iC)) & ~isnan(errorY(:, iC));
                    if ~any(mask)
                        continue;
                    end
                else
                    mask = true(size(xValues));
                end
                x = xValues(mask) + iC*p.Results.xJitter;
                y = meanY(mask, iC);
                ey = errorY(mask, iC);

                h(iC) = TrialDataUtilities.Plotting.errorline(x, y, ey, ...
                    'Color', app(1, iC).Color, 'axh', axh, 'LineWidth', p.Results.LineWidth);
                TrialDataUtilities.Plotting.showInLegend(h(iC), conditionNames{iC});
            end

            td.prepareAxesForChannels(nameX, nameY, 'axh', axh);

            AutoAxis.updateIfInstalled(axh);
        end

        function plotParamVsParamCovariance(td, nameX, nameY, varargin)
            p = inputParser();
            p.addParameter('LineWidth', 1, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            dataX = td.getParamGrouped(nameX);
            dataY = td.getParamGrouped(nameY);

            axh = td.getRequestedPlotAxis(p.Unmatched);
            hold(axh, 'on');

            h = nan(td.nConditions, 1);
            for iC = 1:td.nConditions
                h(iC) = TrialDataUtilities.Plotting.plotCovarianceEllipse(dataX{iC}, dataY{iC}, ...
                    'Color', td.conditionAppearances(iC).Color, 'axh', axh, 'LineWidth', p.Results.LineWidth);
                TrialDataUtilities.Plotting.showInLegend(h(iC), td.conditionNames{iC});
            end

            td.prepareAxesForChannels(nameX, nameY, 'axh', axh);

            AutoAxis.updateIfInstalled(axh);
        end
    end

    % Spike data related
    methods
        % return aligned unit spike times
        function [timesCell, timeScaling] = getSpikeTimes(td, unitNames, varargin)
            % nTrials x nUnits
            p = inputParser();
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('combine', false, @islogical);
            p.addParameter('slice', [], @(x) true); % for subselecting units from array
            p.addParameter('combineAligns', false, @islogical);
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('applyScaling', true, @islogical);
            p.parse(varargin{:});

            [timesCell, timeScaling] = getSpikeTimes@TrialData(td, unitNames, 'combine', p.Results.combine, ...
                'slice', p.Results.slice, 'applyScaling', p.Results.applyScaling);
            if p.Results.applyScaling
                % convert everything to double and in ms
                preserveNativeScaling = 0;
                singleTimepointTolerance = 0;
                edgeTolerance = 0;
            else
                % leave as is, requires multiply by timeScaling to get ms
                preserveNativeScaling = timeScaling;
                % with preserve native scaling, we can be a bit more precise about the edgeTolerance here
                nativeClass = TrialDataUtilities.Data.getCellElementClass(timesCell);
                edgeTolerance = double(ones(1, nativeClass)) * 1e-3;
                singleTimepointTolerance = edgeTolerance;
            end

            if ~p.Results.combineAligns
                timesCell = td.alignInfoActive.getAlignedTimesCell(timesCell, p.Results.includePadding, ...
                    'singleTimepointTolerance', singleTimepointTolerance, 'edgeTolerance', edgeTolerance, ...
                    'preserveNativeScaling', preserveNativeScaling);
            else
                alignIdx = p.Results.alignIdx;
                nAlign = numel(alignIdx);
                timesCellAlign = cell([size(timesCell), nAlign]);
                
                alignOffsets = td.getAlignPlottingTimeOffsets([], 'alignIdx', alignIdx); % use full windows of valid trials, not set by spikes
                for iA = 1:nAlign
                    timesCellAlign(:, :, iA) = td.alignInfoSet{iA}.getAlignedTimesCell(timesCell, p.Results.includePadding, ...
                    'singleTimepointTolerance', singleTimepointTolerance, 'edgeTolerance', edgeTolerance, ...
                    'preserveNativeScaling', preserveNativeScaling);
                end
                
                nTrials = size(timesCell, 1);
                nUnits = size(timesCell, 2);
                for iT = 1:nTrials
                    for iU = 1:nUnits
                        for iA = 1:nAlign
                            timesCellAlign{iT, iU, iA} = timesCellAlign{iT, iU, iA} + alignOffsets(iA);
                        end
                        timesCell{iT, iU} = cat(1, timesCellAlign{iT, iU, :});
                    end
                end
            end
                
        end

        function timesCell = getSpikeTimesEachAlign(td, unitNames, varargin)
            % nTrials x nUnits x nAlign
            p = inputParser();
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            timesCellUnaligned = td.getSpikeTimesUnaligned(unitNames, 'combine', p.Results.combine);

            timesCell = cell([size(timesCellUnaligned) td.nAlign]);
            for iA = 1:td.nAlign
                timesCell(:, :, iA) = td.alignInfoSet{iA}.getAlignedTimesCell(timesCellUnaligned, p.Results.includePadding, 'singleTimepointTolerance', 0);
            end
        end

        function [timesMask, indFirst, indLast] = getSpikeTimesMask(td, unitNames, varargin)
            p = inputParser();
            p.addParameter('includeInvalid', true, @islogical); % if true, include a mask on the invalid trials that rejects everything
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            if p.Results.includeInvalid
                timesRaw = td.getRawSpikeTimes(unitNames, 'combine', p.Results.combine);
            else
                timesRaw = td.getSpikeTimesUnaligned(unitNames, 'combine', p.Results.combine);
            end

            [timesMask, indFirst, indLast] = td.alignInfoActive.getAlignedTimesMask(...
                timesRaw, 'includePadding', p.Results.includePadding, ...
                'singleTimepointTolerance', 0, 'edgeTolerance', 0);
        end

        %%%%%
        % Spike blank intervals
        %%%%%

        function intervalCell = getSpikeBlankingRegions(td, unitNames, varargin)
            p = inputParser();
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            intervalCell = getSpikeBlankingRegions@TrialData(td, unitNames, 'combine', p.Results.combine);

            % align the intervals to the current align info
            if ~isempty(intervalCell)
                intervalCell = td.alignInfoActive.getAlignedIntervalCell(intervalCell, p.Results.includePadding);
            end
        end

        function td = blankSpikesWithinAlignWindow(td, spikeCh)
            td.warnIfNoArgOut(nargout);
            % build interval cell with current alignment window

            [start, stop] = td.alignInfoActive.getStartStopRelativeToZeroByTrial();
            intervalCell = num2cell([start, stop], 2);

            td = td.blankSpikesInTimeIntervals(spikeCh, intervalCell);
        end

        %%%%%%%
        % Binned spike counts
        %%%%%%%

        function [counts, tvec, hasSpikes, tBinEdges, td] = getSpikeBinnedCounts(td, unitName, varargin)
            % if raggedCell is false (default): counts is nTrials x T x nUnits, tvec is T x 1
            % if raggedCell is true (default): counts is nTrials { T x nUnits }, tvec is nTrials {T x 1}
            % is nTrials x nUnits
            p = inputParser;
            p.addParameter('binWidthMs', 1, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Causal, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('combine', false, @islogical);
            p.addParameter('slice', [], @(x) true); % for subselecting units from array
            p.addParameter('raggedCell', false, @islogical); % false, returns countsTensor
            p.parse(varargin{:});
            binWidth = p.Results.binWidthMs;
            binAlignmentMode = p.Results.binAlignmentMode;
            raggedCell = p.Results.raggedCell;

            % pad a bit forward or backwards depending on binning
            td = td.padForTimeBinning(binWidth, binAlignmentMode, false, true);
            spikeCell = td.getSpikeTimes(unitName, 'includePadding', true, 'combine', p.Results.combine, 'slice', p.Results.slice);

            % provide an indication as to which trials have spikes
            hasSpikes = ~cellfun(@isempty, spikeCell);
            valid = td.valid;

            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            % use padded times to incorporate the padding for the binning
            % we used
            timeInfo = td.alignInfoActive.timeInfo;
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];

            [tvec_common, tbinsForHistc, tbinsValidMat] = binAlignmentMode.generateCommonBinnedTimeVector(...
                tMinByTrial, tMaxByTrial, binWidth);

            % get spike blanking regions
            %blankIntervals = td.getSpikeBlankingRegions(unitName);
            blankIntervals = {};

            % use hist c to bin the spike counts
            nUnits = size(spikeCell, 2);
            nTrials = size(spikeCell, 1);
            
            if ~raggedCell
                tvec = tvec_common;
                counts = nan(nTrials, numel(tvec), nUnits);
                for iTrial = 1:nTrials
                    if valid(iTrial)
                        for iU = 1:nUnits
                            if ~isempty(spikeCell{iTrial, iU})
                                counts(iTrial, :, iU) = histcounts(spikeCell{iTrial, iU}, tbinsForHistc);
                            else
                                counts(iTrial, :, iU) = zeros(1, numel(tvec));
                            end

                            % update tbinsValidMat to reflect spike blanking intervals
                            tvecValidMask = makecol(tbinsValidMat(iTrial, :));

                            if ~isempty(blankIntervals)
                                bi = blankIntervals{iTrial};
                                for iInt = 1:size(bi, 1)
                                    % keep regions where the time bin starts after the
                                    % interval or ends before the interval
                                    tvecValidMask = tvecValidMask & makecol(tbinsForHistc(1:end-1) >= bi(iInt, 2) | tbinsForHistc(2:end) <= bi(iInt, 1));
                                end
                            end

                            % mark NaN bins not valid on that trial
                            counts(iTrial, ~tvecValidMask, iU) = NaN;
                        end
                    end
                end
            else
                [counts, tvec] = deal(cell(nTrials, 1));
                for iTrial = 1:nTrials
                    if valid(iTrial)
                        % tbinsValidMat is nTrials x nTime, find first and last time
                        tindFirst = find(tbinsValidMat(iTrial, :), 1, 'first');
                        tindLast = find(tbinsValidMat(iTrial, :), 1, 'last');
                        
                        nTimeThis = tindLast - tindFirst + 1;
                        tvec{iTrial} = tvec_common(tindFirst:tindLast);
                        counts{iTrial} = zeros(nTimeThis, nUnits);
                        
                        for iU = 1:nUnits
                            if ~isempty(spikeCell{iTrial, iU})
                                counts{iTrial}(:, iU) = histcounts(spikeCell{iTrial, iU}, tbinsForHistc(tindFirst:tindLast+1))';
                            end

%                             % update tbinsValidMat to reflect spike blanking intervals
%                             tvecValidMask = makecol(tbinsValidMat(iTrial, :));
% 
%                             if ~isempty(blankIntervals)
%                                 bi = blankIntervals{iTrial};
%                                 for iInt = 1:size(bi, 1)
%                                     % keep regions where the time bin starts after the
%                                     % interval or ends before the interval
%                                     tvecValidMask = tvecValidMask & makecol(tbinsForHistc(1:end-1) >= bi(iInt, 2) | tbinsForHistc(2:end) <= bi(iInt, 1));
%                                 end
%                             end
% 
%                             % mark NaN bins not valid on that trial
%                             counts(iTrial, ~tvecValidMask, iU) = NaN;
                        end
                    else
                        counts{iTrial} = zeros(0, nUnits);
                    end
                end
                 
            end

            tBinEdges = tbinsForHistc;
        end

        function [countsMat, tvec, hasSpikes, tBinEdgesCell, alignVec] = getSpikeBinnedCountsEachAlign(td, unitName, varargin)
            countsCell = cellvec(td.nAlign);
            tvecCell = cellvec(td.nAlign);
            hasSpikesCell = cellvec(td.nAlign);
            tBinEdgesCell = cellvec(td.nAlign);
            for iA = 1:td.nAlign
                 [countsCell{iA}, tvecCell{iA}, hasSpikesCell{iA}, tBinEdgesCell{iA}] =  td.useAlign(iA).getSpikeBinnedCounts(unitName, varargin{:});
            end

            [countsMat, alignVec] = TensorUtils.catWhich(2, countsCell{:});
            tvec = cat(1, tvecCell{:});
            hasSpikesMat = cat(2, hasSpikesCell{:});
            hasSpikes = any(hasSpikesMat, 2);
        end

        function [countsGrouped, tvec, hasSpikesGrouped, tBinEdges] = getSpikeBinnedCountsGrouped(td, unitName, varargin)
            [countsMat, tvec, hasSpikes, tBinEdges, td] = td.getSpikeBinnedCounts(unitName, varargin{:}); % we take back td so that we keep the validity changed during padding
            countsGrouped = td.groupElements(countsMat);
            hasSpikesGrouped = td.groupElements(hasSpikes);
        end

        function [countsGrouped, tvec, hasSpikesGrouped, tBinEdges] = getSpikeBinnedCountsGroupedRandomized(td, unitName, varargin)
            [countsMat, tvec, hasSpikes, tBinEdges, td] = td.getSpikeBinnedCounts(unitName, varargin{:});
            countsGrouped = td.groupElementsRandomized(countsMat);
            hasSpikesGrouped = td.groupElementsRandomized(hasSpikes);
        end

        function [countsGrouped, tvec, hasSpikesGrouped, tBinEdges, alignVec] = getSpikeBinnedCountsGroupedEachAlign(td, unitName, varargin)
            [countsMat, tvec, hasSpikes, tBinEdges, alignVec] = td.getSpikeBinnedCountsEachAlign(unitName, varargin{:});
            countsGrouped = td.groupElements(countsMat);
            hasSpikesGrouped = td.groupElements(hasSpikes);
        end

        function [countsGrouped, tvec, hasSpikesGrouped, tBinEdges, alignVec] = getSpikeBinnedCountsGroupedEachAlignRandomized(td, unitName, varargin)
            [countsMat, tvec, hasSpikes, tBinEdges, alignVec] = td.getSpikeBinnedCountsEachAlign(unitName, varargin{:});
            countsGrouped = td.groupElementsRandomized(countsMat);
            hasSpikesGrouped = td.groupElementsRandomized(hasSpikes);
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat, tBinEdges] = ...
                getSpikeBinnedCountsGroupMeans(td, varargin)
            % *Mat will be nConditions x T matrices
            p = inputParser();
            p.addRequired('unitName', @(x) isstringlike(x) || iscellstr(x) || isstring(x));
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('removeZeroSpikeTrials', false, @islogical);
            p.addParameter('assumePoissonStatistics', false, @islogical); % comutes the sem and std assuming Poisson statistics (i.e. var == mean)
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unitName = p.Results.unitName;

            if isempty(p.Results.minTrials)
                minTrials = 1;
            else
                minTrials = p.Results.minTrials;
            end

            if isempty(p.Results.minTrialFraction)
                minTrialFraction = 0;
            else
                minTrialFraction = p.Results.minTrialFraction;
            end

            [countsGrouped, tvec, hasSpikesGrouped, tBinEdges] = td.getSpikeBinnedCountsGrouped(unitName, p.Unmatched);

            nUnits = size(countsGrouped{1}, 3);

            % remove trials from each group that have no spikes by setting
            % the whole trial to NaN on per-unit basis
            if p.Results.removeZeroSpikeTrials
                for iC = 1:td.nConditions
                    for iU = 1:nUnits
                        countsGrouped{iC}(hasSpikesGrouped{iC}(:, iU), :, iU) = NaN;
                    end
                end
            end

            % comute the means
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec), nUnits));
            for iC = 1:td.nConditions
                if ~isempty(countsGrouped{iC})
                    [psthMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), stdMat(iC, :, :)] = ...
                        TrialDataUtilities.Data.nanMeanSemMinCount(countsGrouped{iC}, 1, minTrials, minTrialFraction, ...
                        'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', 1); % poisson multiplier is 1 since these are raw counts
                end
            end
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat] = ...
                getSpikeBinnedCountsGroupMeansEachAlign(td, varargin)
            % *Mat will be nConditions x T matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addRequired('unitName', @isstringlike);
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('removeZeroSpikeTrials', false, @islogical);
            p.addParameter('assumePoissonStatistics', false, @islogical); % comutes the sem and std assuming Poisson statistics (i.e. var == mean)
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unitName = p.Results.unitName;

            if isempty(p.Results.minTrials)
                minTrials = 1;
            else
                minTrials = p.Results.minTrials;
            end

            if isempty(p.Results.minTrialFraction)
                minTrialFraction = 0;
            else
                minTrialFraction = p.Results.minTrialFraction;
            end

            [countsGrouped, tvec, hasSpikesGrouped] = td.getSpikeBinnedCountsGroupedEachAlign(unitName, p.Unmatched);
            nUnits = size(countsGrouped{1}, 3);

            % remove trials from each group that have no spikes by setting
            % the whole trial to NaN on per-unit basis
            if p.Results.removeZeroSpikeTrials
                for iC = 1:td.nConditions
                    for iU = 1:nUnits
                        countsGrouped{iC}(hasSpikesGrouped{iC}(:, iU), :, iU) = NaN;
                    end
                end
            end

            % comute the means
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec), nUnits));
            for iC = 1:td.nConditions
                if ~isempty(countsGrouped{iC})
                    [psthMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), stdMat(iC, :, :)] = ...
                        TrialDataUtilities.Data.nanMeanSemMinCount(countsGrouped{iC}, 1, minTrials, minTrialFraction, ...
                        'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', 1); % poisson multiplier is 1 since these are raw counts
                end
            end
        end

        function plotSpikeBinBoundaries(td, unitName, varargin)
            % utility function for showing raster with superimposed spike
            % bin boundaries for debugging bin widths and bin alignment
            % mode
            p = inputParser;
            p.addParameter('binWidthMs', 1, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Causal, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});
            binWidth = p.Results.binWidthMs;
            binAlignmentMode = p.Results.binAlignmentMode;

            [~, ~, ~, tBinEdges] = td.getSpikeBinnedCounts(unitName, ...
                'binWidthMs', binWidth, 'binAlignmentMode', binAlignmentMode, ...
                'combine', p.Results.combine);

            td.plotRaster(unitName, 'combine', p.Results.combine);
            hold on;
            TrialDataUtilities.Plotting.verticalLine(tBinEdges, 'Color', grey(0.5));
            hold off;
        end

    end

    methods % Filtered spike rates
        function [rateCell, timeCell, hasSpikes, poissonCountMultiplier] = getSpikeRateFiltered(td, unitName, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            sf = p.Results.spikeFilter;

            % Pad trial data alignment for spike filter
            td = td.padForSpikeFilter(sf);

            spikeCell = td.getSpikeTimes(unitName, 'includePadding', true, 'combine', p.Results.combine);
            timeInfo = td.alignInfoActive.timeInfo;

            % provide an indication as to which trials have spikes
            hasSpikes = ~cellfun(@isempty, spikeCell);

            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rateCell, timeCell, poissonCountMultiplier] = sf.filterSpikeTrainsWindowByTrial(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond);

            % now we need to nan out the regions affected by blanking
            blankingIntervals = td.getSpikeBlankingRegions(unitName);
            [rateCell, timeCell] = TrialDataUtilities.SpikeData.markNanBlankedIntervals(...
                blankingIntervals, rateCell, timeCell, 'padding', [sf.preWindow sf.postWindow]);
        end

        function [rates, tvec, hasSpikes, poissonCountMultiplier] = getSpikeRateFilteredAsMatrix(td, unitNames, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('combine', false, @islogical);
            p.addParameter('showProgress', true, @islogical);
            p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); % manually dictate time boundaries if specified, otherwise auto
            p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); 
            p.addParameter('useNativeScaling', false, @islogical);
            p.parse(varargin{:});

            sf = p.Results.spikeFilter;
            if isempty(sf)
                sf = SpikeFilter.getDefaultFilter();
            end

            useNativeScaling = p.Results.useNativeScaling;

            % Pad trial data alignment for spike filter
            td = td.padForSpikeFilter(sf);

            if isnumeric(unitNames)
                unitNames = td.lookupSpikeChannelByIndex(unitNames);
            end

            % critical to include the spike times in the padded window
            [spikeCell, timeScaling] = td.getSpikeTimes(unitNames, 'includePadding', true, 'combine', p.Results.combine, 'applyScaling', ~useNativeScaling);
            if useNativeScaling
                preserveNativeScaling = timeScaling;
            else
                preserveNativeScaling = 0;
            end
            [tMinByTrialExcludingPadding, tMaxByTrialExcludingPadding] = td.alignInfoActive.getStartStopRelativeToZeroByTrial();
            [tMinByTrialWithPadding, tMaxByTrialWithPadding] = td.alignInfoActive.getStartStopRelativeToZeroByTrialWithPadding();
            
            % provide an indication as to which trials have spikes
            hasSpikes = ~cellfun(@isempty, spikeCell);

            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            [rates, tvec, poissonCountMultiplier] = sf.filterSpikeTrainsWindowByTrialAsMatrix(spikeCell, ...
                tMinByTrialWithPadding, tMaxByTrialWithPadding, td.timeUnitsPerSecond, ...
                'showProgress', p.Results.showProgress, ...
                'tMinByTrialExcludingPadding', tMinByTrialExcludingPadding, ...
                'tMaxByTrialExcludingPadding', tMaxByTrialExcludingPadding, ...
                'tMin', p.Results.tMin, ...
                'tMax', p.Results.tMax, ...
                'preserveNativeScaling', preserveNativeScaling);
            tvec = makecol(tvec);

            % now we need to nan out the regions affected by blanking
%             blankingIntervals = td.getSpikeBlankingRegions(unitNames, 'combine', p.Results.combine);
%             if ~isempty(blankingIntervals)
%                 [rates, tvec] = TrialDataUtilities.SpikeData.markNanBlankedIntervals(...
%                     blankingIntervals, rates, tvec, 'padding', [sf.preWindow sf.postWindow]);
%             end
        end

        function [rates, tvec, hasSpikes, alignVec, poissonCountMultiplier] = getSpikeRateFilteredAsMatrixEachAlign(td, unitName, varargin)
            p = inputParser;
            p.addParameter('combine', false, @islogical);
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('tMin', [], @(x) isempty(x) || isvector(x)); % manually dictate time boundaries if specified, otherwise auto
            p.addParameter('tMax', [], @(x) isempty(x) || isvector(x)); 
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alignIdx = p.Results.alignIdx;

            if isnumeric(unitName)
                unitName = td.lookupSpikeChannelByIndex(unitNames);
            end
            unitName = string(unitName);
            if p.Results.combine
                nUnits = 1;
            else
                nUnits = td.computeSpikeUnitCountFromName(unitName);
            end

            nAlign = numel(alignIdx);
            ratesCell = cell(nAlign, 1);
            tvecCell = cell(nAlign, 1);
            hasSpikesMat = nan(td.nTrials, nAlign, nUnits);
            for iA = 1:nAlign
                if isempty(p.Results.tMin)
                    tMinThis = [];
                else
                    tMinThis = p.Results.tMin(iA);
                end
                if isempty(p.Results.tMax)
                    tMaxThis = [];
                else
                    tMaxThis = p.Results.tMax(iA);
                end
                [ratesCell{iA}, tvecCell{iA}, hasSpikesMat(:, iA, :), poissonCountMultiplier] = ...
                    td.useAlign(alignIdx(iA)).getSpikeRateFilteredAsMatrix(unitName, 'combine', p.Results.combine, ...
                    'tMin', tMinThis, 'tMax', tMaxThis, ...
                    p.Unmatched);
            end

            [tvec, alignVec] = td.catTimeVecOverAlignWithSeparator(tvecCell, 'alignIdx', p.Results.alignIdx);
            rates = td.catDataOverAlignWithSeparator(2, ratesCell);
            hasSpikes = any(hasSpikesMat, 2);
        end

        function [rateCell, timeCell, hasSpikesGrouped] = getSpikeRateFilteredGrouped(td, unitName, varargin)
            [rateCell, timeCell, hasSpikes] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElementsFlat(rateCell);
            timeCell = td.groupElementsFlat(timeCell);
            hasSpikesGrouped = td.groupElementsFlat(hasSpikes);
        end

        function [rateCell, timeCell, hasSpikesGrouped, poissonCountMultiplier] = getSpikeRateFilteredGroupedRandomized(td, unitName, varargin)
            [rateCell, timeCell, poissonCountMultiplier] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElementsFlatRandomized(rateCell);
            timeCell = td.groupElementsFlatRandomized(timeCell);
            hasSpikesGrouped = td.groupElementsFlatRandomized(hasSpikes);
        end

        function [rateCell, tvec, hasSpikesGrouped, poissonCountMultiplier] = getSpikeRateFilteredAsMatrixGrouped(td, unitNames, varargin)
            [rates, tvec, hasSpikes, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrix(unitNames, varargin{:});
            rateCell = td.groupElementsFlat(rates);
            hasSpikesGrouped = td.groupElementsFlat(hasSpikes);
        end

        function [rateCell, tvec, hasSpikesGrouped, poissonCountMultiplier] = getSpikeRateFilteredAsMatrixGroupedRandomized(td, unitNames, varargin)
            [rates, tvec, hasSpikes, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrix(unitNames, varargin{:});
            rateCell = td.groupElementsFlatRandomized(rates);
            hasSpikesGrouped = td.groupElementsFlatRandomized(hasSpikes);
        end

        function [rateCell, tvec, hasSpikesGrouped, alignVec, poissonCountMultiplier] = getSpikeRateFilteredAsMatrixGroupedEachAlign(td, unitName, varargin)
            [rates, tvec, hasSpikes, alignVec, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrixEachAlign(unitName, varargin{:});
            [rateCell, hasSpikesGrouped] = td.groupElementsFlat(rates, hasSpikes);
        end

        function [rateCell, tvec, hasSpikesGrouped, alignVec, poissonCountMultiplier] = getSpikeRateFilteredAsMatrixGroupedRandomizedEachAlign(td, unitName, varargin)
            [rates, tvec, hasSpikes, alignVec, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrixEachAlign(unitName, varargin{:});
            [rateCell, hasSpikesGrouped] = td.groupElementsFlatRandomized(rates, hasSpikes);
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                getSpikeRateFilteredGroupMeans(td, unitNames, varargin)
            % *Mat will be nConditions x T matrices
            % if randomized, will be nConditions x T x nRandomized
            p = inputParser();
            p.addRequired('unitNames', @(x) isnumeric(x) || ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('subtractConditionDescriptor', [], @(x) isempty(x) || isa(x, 'ConditionDescriptor'));

            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('removeZeroSpikeTrials', false, @islogical);
            p.addParameter('assumePoissonStatistics', false, @islogical);

            p.addParameter('combine', false, @islogical);
            p.addParameter('useNativeScaling', false, @islogical); % do computations within the native scaling of the spike times for precision

            % these are to enable this function to support multiple roles
            p.addParameter('eachAlign', false, @islogical);
            
            % manually dictate time boundaries if specified, otherwise auto
            p.addParameter('tMin', [], @(x) isempty(x) || isvector(x)); 
            p.addParameter('tMax', [], @(x) isempty(x) || isvector(x)); 
            p.parse(unitNames, varargin{:});
            unitNames = p.Results.unitNames;
            
            if isnumeric(unitNames)
                unitNames = td.lookupSpikeChannelByIndex(unitNames);
            end
            
            if isempty(p.Results.minTrials)
                minTrials = 1;
            else
                minTrials = p.Results.minTrials;
            end

            if isempty(p.Results.minTrialFraction)
                minTrialFraction = 0;
            else
                minTrialFraction = p.Results.minTrialFraction;
            end

            % pick the right function to call
            if p.Results.eachAlign
                 [rateCell, tvec, hasSpikesGrouped, whichAlign, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrixGroupedEachAlign(unitNames, ...
                     'spikeFilter', p.Results.spikeFilter, 'combine', p.Results.combine, 'tMin', p.Results.tMin, 'tMax', p.Results.tMax, 'useNativeScaling', p.Results.useNativeScaling);
            else
                [rateCell, tvec, hasSpikesGrouped, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrixGrouped(unitNames, ...
                    'spikeFilter', p.Results.spikeFilter, 'combine', p.Results.combine, 'tMin', p.Results.tMin, 'tMax', p.Results.tMax, 'useNativeScaling', p.Results.useNativeScaling);
                whichAlign = td.alignInfoActiveIdx * ones(size(tvec));
            end

%             if p.Results.combine || ischar(unitNames) || isStringScalar(unitNames)
%                 cd = td.getChannelDescriptor(unitNames);
%                 if isa(cd, 'SpikeArrayChannelDescriptor')
%                     nUnits = cd.nChannels;
%                 else
%                     nUnits = 1;
%                 end
%             else
%                 nUnits = numel(unitNames);
%             end

            nUnits = size(rateCell{1}, 3);

            % remove trials from each group that have no spikes
            if p.Results.removeZeroSpikeTrials
                for iC = 1:td.nConditions
                    for iU = 1:nUnits
                        rateCell{iC}(~hasSpikesGrouped{iC}(:, iU), :, iU) = NaN;
                    end
                end
            end

            % comute the means
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec), nUnits));
            for iC = 1:td.nConditions
                if ~isempty(rateCell{iC})
                    [psthMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), stdMat(iC, :, :)] = ...
                        TrialDataUtilities.Data.nanMeanSemMinCount(rateCell{iC}, 1, minTrials, minTrialFraction, ...
                        'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', poissonCountMultiplier);
                end
            end

            % subtract these conditions piecemeal from this one
            if ~isempty(p.Results.subtractConditionDescriptor)
                cdSub = td.expandConditionDescriptorForSubtraction(p.Results.subtractConditionDescriptor);

                % get the group means for this condition descriptor
                if p.Results.eachAlign
                    [sub_psthMat, sub_tvec, sub_semMat, sub_stdMat, ~, sub_whichAlign] = ...
                        td.setConditionDescriptor(cdSub).getSpikeRateFilteredGroupMeansEachAlign(unitNames, ...
                        rmfield(p.Results, 'subtractConditionDescriptor'), p.Unmatched);

                    % split by align, equalize time vectors for each align,
                    % and recombine
                    psthMat = TensorUtils.splitAlongDimensionByIndex(psthMat, 2, whichAlign);
                    sub_psthMat = TensorUtils.splitAlongDimensionByIndex(sub_psthMat, 2, sub_whichAlign);
                    tvec = TensorUtils.splitAlongDimensionByIndex(tvec, 1, whichAlign);
                    sub_tvec = TensorUtils.splitAlongDimensionByIndex(sub_tvec, 1, sub_whichAlign);
                    semMat = TensorUtils.splitAlongDimensionByIndex(semMat, 2, whichAlign);
                    sub_semMat = TensorUtils.splitAlongDimensionByIndex(sub_semMat, 2, sub_whichAlign);
                    stdMat = TensorUtils.splitAlongDimensionByIndex(stdMat, 2, whichAlign);
                    sub_stdMat = TensorUtils.splitAlongDimensionByIndex(sub_stdMat, 2, sub_whichAlign);

                    V = 6;
                    nA = numel(psthMat);
                    outByAlign = cell(V, nA);
                    tvecByAlign = cell(nA, 1);
                    for iA = 1:numel(psthMat)
                        [outByAlign(:, iA), tvecByAlign{iA}] = TrialDataUtilities.Data.equalizeTimeVectorsForTimeseries({...
                            psthMat{iA}, semMat{iA}, stdMat{iA}, sub_psthMat{iA}, sub_semMat{iA}, sub_stdMat{iA}}, ...
                            {tvec{iA}, tvec{iA}, tvec{iA}, sub_tvec{iA}, sub_tvec{iA}, sub_tvec{iA}}, 2);
                    end
                    out = cell(V, 1);
                    for i = 1:Vdbquit
                        out{i} = cat(2, outByAlign{i, :});
                    end
                    tvec = cat(1, tvecByAlign{:});
                else
                    [sub_psthMat, sub_tvec, sub_semMat, sub_stdMat] = ...
                        td.setConditionDescriptor(cdSub).getSpikeRateFilteredGroupMeans(unitNames, ...
                        rmfield(p.Results, 'subtractConditionDescriptor'), p.Unmatched);

                    % equalize the time vectors
                    [out, tvec] = TrialDataUtilities.Data.equalizeTimeVectorsForTimeseries({...
                        psthMat, semMat, stdMat, sub_psthMat, sub_semMat, sub_stdMat}, ...
                        {tvec, tvec, tvec, sub_tvec, sub_tvec, sub_tvec}, 2);
                end

                % do the subtraction
                psthMat = out{1} - out{4};

                % sd = sqrt(sd1^2 + sd2^2)
                % sem = sqrt(sem1^2 + sem2^2)
                semMat = sqrt(out{2}.^2 + out{5}.^2);
                stdMat = sqrt(out{3}.^2 + out{6}.^2);
            end
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                getSpikeRateFilteredGroupMeansEachAlign(td, unitName, varargin)
            [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                td.getSpikeRateFilteredGroupMeans(unitName, 'eachAlign', true, varargin{:});
        end

        function [snr, range, noise]  = getSpikeChannelSNR(td, unitName, varargin)
            % take range = max - min (or e.g. 0.0.25, 0.975 quantile) and noise = max (or 0.95 quantile) s.e.m.
            p = inputParser();
            p.addParameter('noiseQuantile', 0.95, @isscalar);
            p.addParameter('rangeQuantile', 0.95, @(x) isscalar(x) || numel(x) == 2); % either width of centered quantile band (e.g. 0.95 --> 0.025 to 0.975) or [low high]
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            nq = p.Results.noiseQuantile;
            rq = p.Results.rangeQuantile;
            if isscalar(rq)
                rq = [(1-rq)/2 1-(1-rq)/2];
            end

            % options same as getSpikeRateFilteredGroupMeansEachAlign
            [psthMat, ~, semMat] = td.getSpikeRateFilteredGroupMeansEachAlign(unitName, p.Unmatched);

            % single trial estimated average values don't count
            psthMat(semMat == 0) = NaN;
            semMat(semMat == 0) = NaN;

            noise = squeeze(TensorUtils.quantileMultiDim(semMat, nq, [1 2]));
            range = squeeze(diff(TensorUtils.quantileMultiDim(psthMat, rq, [1 2]), 1, 1));
            snr = range ./ noise;
        end

        function [snr, range, noise]  = getSpikeChannelSNROverConditions(td, unitName, varargin)
            % take range = max - min over conditions at each time and noise = max s.e.m. at each time
            % then take quantile of these snr values over time
            p = inputParser();
            p.addParameter('quantile', 0.95, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            q = p.Results.quantile;
            [psthMat, ~, semMat] = td.getSpikeRateFilteredGroupMeansEachAlign(unitName, p.Unmatched);
            range = max(psthMat,[],1, 'omitnan') - min(psthMat,[],1, 'omitnan');
            noise = max(semMat, [], 1, 'omitnan');
            snr = range ./ noise;
            snr = quantile(snr, q);
        end
        
        function [snr_c, range_c, noise_c] = getSpikeChannelSNREachCondition(td, unitName, varargin)
            p = inputParser();
            p.addParameter('noiseQuantile', 0.95, @isscalar);
            p.addParameter('rangeQuantile', 0.95, @(x) isscalar(x) || numel(x) == 2); % either width of centered quantile band (e.g. 0.95 --> 0.025 to 0.975) or [low high]
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            nq = p.Results.noiseQuantile;
            rq = p.Results.rangeQuantile;
            if isscalar(rq)
                rq = [(1-rq)/2 1-(1-rq)/2];
            end

            % options same as getSpikeRateFilteredGroupMeansEachAlign
            [psthMat, ~, semMat] = td.getSpikeRateFilteredGroupMeansEachAlign(unitName, p.Unmatched);

            % single trial estimated average values don't count
            psthMat(semMat == 0) = NaN;
            semMat(semMat == 0) = NaN;
  
            noise_c = squeeze(TensorUtils.quantileMultiDim(semMat, nq, 2));
            range_c = squeeze(diff(TensorUtils.quantileMultiDim(psthMat, rq, 2), 1, 2));
            snr_c = range_c ./ noise_c;
        end
            

        function tbl = getSpikeChannelSNRTable(td, varargin)
            units = td.listSpikeChannels();

            [snr, range, noise] = td.getSpikeChannelSNR(units, varargin{:}, 'combine', false);

            s = struct('snr', num2cell(snr), 'range', num2cell(range), 'noise', num2cell(noise));
            tbl = struct2table(s, 'RowNames', units, 'AsArray', true);
        end

        %   overConditions: at each time, take max - min over conditions at each time, and noise = max s.e.m at each time, then ratio, and take


        function [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                getSpikeRateFilteredGroupMeansRandomized(td, unitNames, varargin)
            % *Mat will be nConditions x T x nRandomized
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addRequired('unitNames', @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('assumePoissonStatistics', false, @islogical);
                
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('removeZeroSpikeTrials', false, @islogical);

            % these are to enable this function to support multiple roles
            p.addParameter('eachAlign', false, @islogical);
            p.parse(unitNames, varargin{:});
            unitNames = p.Results.unitNames;

            if isempty(p.Results.minTrials)
                minTrials = 1;
            else
                minTrials = p.Results.minTrials;
            end

            if isempty(p.Results.minTrialFraction)
                minTrialFraction = 0;
            else
                minTrialFraction = p.Results.minTrialFraction;
            end

            % grab the data ungrouped
            if ~p.Results.eachAlign
                [rates, tvec, hasSpikes, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrix(unitNames,...
                    'spikeFilter', p.Results.spikeFilter, 'timeDelta', p.Results.timeDelta);
                whichAlign = td.alignInfoActiveIdx * ones(size(tvec));
            else
                [rates, tvec, hasSpikes, whichAlign, poissonCountMultiplier] = td.getSpikeRateFilteredAsMatrixEachAlign(unitNames, ...
                    'spikeFilter', p.Results.spikeFilter, 'timeDelta', p.Results.timeDelta);
            end

            % now do grouping and mean computation in loop for each iRandom
            % so as to save memory
            prog = ProgressBar(td.nRandomized, 'Computing mean spike rate over randomizations');
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec), td.nRandomized));

            for iR = 1:td.nRandomized
                prog.update(iR);

                % group for this randomization
                [rateCell, hasSpikesGrouped] = td.groupElementsFlatRandomizedSingle(iR, rates, hasSpikes);

                % remove trials from each group that have no spikes
                if p.Results.removeZeroSpikeTrials
                    for iC = 1:td.nConditions
                        rateCell{iC} = rateCell{iC}(hasSpikesGrouped{iC}, :);
                    end
                end

                % compute the means
                for iC = 1:td.nConditions
                    if ~isempty(rateCell{iC})
                        [psthMat(iC, :, iR), semMat(iC, :, iR), nTrialsMat(iC, :, iR), stdMat(iC, :, iR)] = ...
                            nanMeanSemMinCount(rateCell{iC}, 1, minTrials, minTrialFraction, ...
                            'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultiplier', poissonCountMultiplier);
                    end
                end
            end
            prog.finish();
        end

        function quantilesByGroup = getSpikeRateFilteredGroupMeansRandomizedQuantiles(td, unitName, varargin)
            p = inputParser();
            p.addParameter('quantiles', [0.025 0.975], @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            meanByGroup = td.getSpikeRateFilteredGroupMeansRandomized(unitName, p.Unmatched);

            quantilesByGroup = quantile(meanByGroup, p.Results.quantiles, ndims(meanByGroup));
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                getSpikeRateFilteredGroupMeansRandomizedEachAlign(td, varargin)
            [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = ...
                td.getSpikeRateFilteredGroupMeansRandomized('eachAlign', true, varargin{:});
        end

        function [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = getPSTH(td, varargin)
            [psthMat, tvec, semMat, stdMat, nTrialsMat, whichAlign] = getSpikeRateFilteredGroupMeans(td, varargin{:});
        end

%         function [snr, range, noise] = getSpikeRateFilteredSNR(td, varargin)
%             % range / max(sem)
%             [psthMat, ~, semMat] = getPSTH(td, varargin{:});
%             noise = nanmax(semMat(:));
%             range = nanmax(psthMat(:)) - nanmin(psthMat(:));
%             snr = range / noise;
%         end

        function timesCellofCells = getSpikeTimesGrouped(td, unitName, varargin)
            timesCell = td.getSpikeTimes(unitName, varargin{:});
            timesCellofCells = td.groupElementsFlat(timesCell);
        end
        
        function timesCU = getSpikeTimesGroupedByConditionUnit(td, unitName, varargin)
            timesCell = td.getSpikeTimes(unitName, varargin{:});
            nU = size(timesCell, 2);
            timesCell = td.groupElementsFlat(timesCell);
            nC = numel(timesCell);
            
            timesCU = cell(nC, nU);
            for iC = 1:nC
                for iU = 1:nU
                    timesCU{iC, iU} = timesCell{iC}(:, iU);
                end
            end
        end

        function timesCellofCells = getSpikeTimesGroupedRandomized(td, unitName, varargin)
            timesCell = td.getSpikeTimes(unitName, varargin{:});
            timesCellofCells = td.groupElementsFlatRandomized(timesCell);
        end

        function countsCell = getSpikeCountsGrouped(td, unitName)
            counts = td.getSpikeCounts(unitName);
            countsCell = td.groupElementsFlat(counts);
        end

        function countsCell = getSpikeCountsGroupedRandomized(td, unitName)
            counts = td.getSpikeCounts(unitName);
            countsCell = td.groupElementsFlatRandomized(counts);
        end

        function [rates, durations, containsBlanked, poissonCountMultipliers] = getSpikeMeanRateAllAlign(td, unitName, varargin)
            p = inputParser();
            p.addParameter('invalidIfBlanked', false, @islogical); % if true, any trial that is partially blanked will be NaN, if false, the blanked region will be ignored and will not contribute to the time window used as the denominator for the rate calculation
            p.addParameter('combine', false, @islogical);
            p.parse(varargin{:});

            [counts, durations] = zerosvec(td.nTrials);
            containsBlanked = falsevec(td.nTrials);

            for iA = 1:td.nAlign
                counts = counts + td.useAlign(iA).getSpikeCounts(unitName, 'combine', p.Results.combine);
                [durations_this, containsBlanked_this] = td.useAlign(iA).getValidDurationsForSpikeChannel(unitName, 'combine', p.Results.combine);
                durations = durations + durations_this;
                containsBlanked = containsBlanked | containsBlanked_this;
            end

            if p.Results.invalidIfBlanked
                durations(containsBlanked) = NaN;
            end
            rates = counts ./ durations * td.timeUnitsPerSecond;
            poissonCountMultipliers = 1 ./ durations * td.timeUnitsPerSecond;
        end

        function [rateCell, durationCell, containsBlankedCell, poissonCountMultipliersCell] = getSpikeMeanRateGrouped(td, unitName, varargin)
            [rates, durations, containsBlanked, poissonCountMultipliers] = td.getSpikeMeanRate(unitName, varargin{:});
            [rateCell, durationCell, containsBlankedCell, poissonCountMultipliersCell] = td.groupElementsFlat(rates, durations, containsBlanked, poissonCountMultipliers);
        end

        function [rateCell, durationCell, containsBlankedCell, poissonCountMultipliersCell] = getSpikeMeanRateGroupedAllAlign(td, unitName, varargin)
            [rates, durations, containsBlanked, poissonCountMultipliersCell] = td.getSpikeMeanRateAllAlign(unitName, varargin{:});
            [rateCell, durationCell, containsBlankedCell, poissonCountMultipliersCell] = td.groupElementsFlat(rates, durations, containsBlanked, poissonCountMultipliersCell);
        end

        function rateCell = getSpikeMeanRateGroupedRandomized(td, unitName, varargin)
            rates = td.getSpikeMeanRate(unitName, varargin{:});
            rateCell = td.groupElementsFlatRandomized(rates);
        end

        function [meanByGroup, semByGroup, stdByGroup, nByGroup] = getSpikeMeanRateGroupMeans(td, unitName, varargin)
             p = inputParser();
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('assumePoissonStatistics', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [rateCell, ~, ~, poissonCountMultipliersCell] = td.getSpikeMeanRateGrouped(unitName, p.Unmatched);
            
            [meanByGroup, semByGroup, nByGroup, stdByGroup] = cellfun(@(rates, poissonCountMultipliers) TrialDataUtilities.Data.nanMeanSemMinCount(rates, 1, p.Results.minTrials, p.Results.minTrialFraction, ...
                'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', poissonCountMultipliers), ...
                rateCell, poissonCountMultipliersCell, 'UniformOutput', false);
            meanByGroup = cat(1, meanByGroup{:});
            semByGroup = cat(1, semByGroup{:});
            nByGroup = cat(1, nByGroup{:});
            stdByGroup = cat(1, stdByGroup{:});
        end

        function [meanByGroup, semByGroup, stdByGroup, nByGroup] = getSpikeMeanRateGroupMeansAllAlign(td, unitName, varargin)
            p = inputParser();
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('assumePoissonStatistics', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [rateCell, ~, ~, poissonCountMultipliersCell] = td.getSpikeMeanRateGroupedAllAlign(unitName, p.Unmatched);
            
            [meanByGroup, semByGroup, nByGroup, stdByGroup] = cellfun(@(rates, poissonCountMultipliers) TrialDataUtilities.Data.nanMeanSemMinCount(rates, 1, p.Results.minTrials, p.Results.minTrialFraction, ...
                'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', poissonCountMultipliers), ...
                rateCell, poissonCountMultipliersCell, 'UniformOutput', true);

%             meanByGroup = cellfun(@nanmean, rateCell);
%             semByGroup = cellfun(@nansem, rateCell);
%             stdByGroup = cellfun(@nanstd, rateCell);
%             nByGroup = cellfun(@(x) nnz(~isnan(x)), rateCell);
        end

        function [meanByGroup, semByGroup, stdByGroup, nByGroup] = getSpikeMeanRateGroupMeansRandomized(td, unitName, varargin)
            p = inputParser();
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('assumePoissonStatistics', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [rateCell, ~, ~,  poissonCountMultipliersCell] = td.getSpikeMeanRateGroupedRandomized(unitName, varargin{:});
            
            [meanByGroup, semByGroup, nByGroup, stdByGroup] = cellfun(@(rates, poissonCountMultipliers) TrialDataUtilities.Data.nanMeanSemMinCount(rates, 1, p.Results.minTrials, p.Results.minTrialFraction, ...
                'assumePoissonStatistics', p.Results.assumePoissonStatistics, 'poissonCountMultipliers', poissonCountMultipliers), ...
                rateCell, poissonCountMultipliersCell, 'UniformOutput', true);

%             meanByGroup = cellfun(@nanmean, rateCell);
%             semByGroup = cellfun(@nansem, rateCell);
%             stdByGroup = cellfun(@nanstd, rateCell);
%             nByGroup = cellfun(@numel, rateCell);
        end

        function quantilesByGroup = getSpikeMeanRateGroupMeansRandomizedQuantiles(td, unitName, varargin)
            p = inputParser();
            p.addParameter('quantiles', [0.025 0.975], @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            meanByGroup = td.getSpikeMeanRateGroupMeansRandomized(unitName, p.Unmatched);

            quantilesByGroup = quantile(meanByGroup, p.Results.quantiles, ndims(meanByGroup));
        end

        function plotTuningCurve(td, unitName, varargin)
            p = inputParser;
            p.addParameter('errorType', 'sem', @isstringlike);
            p.addParameter('LineWidth', 1, @isscalar);
            p.KeepUnmatched = true;
            p.CaseSensitive = false;
            p.parse(varargin{:});

            tdFlat = td.reshapeAxes(1); % keep only first axis to get appearances

            axh = td.getRequestedPlotAxis(p.Unmatched);
            td = td.flattenAxesExceptFirst();

            [meanByGroup, semByGroup, stdByGroup, nByGroup] = td.getSpikeMeanRateGroupMeans(unitName); %#ok<ASGLU>
            if strcmp(p.Results.errorType, 'sem')
                errorY = semByGroup;
            elseif strcmp(p.Results.errorType, 'std')
                errorY = stdByGroup;
            else
                error('Unknown errorType %s', p.Results.errorType);
            end

            nFirst = size(td.conditions, 1);
            nCurves = td.nConditions / nFirst;
            xVals = 1:nFirst;

            for iC = 1:nCurves
                app = tdFlat.conditionAppearances(iC);
                h = TrialDataUtilities.Plotting.errorline(xVals, meanByGroup(:, iC), errorY(:, iC), ...
                    'Color', app(iC).Color, 'axh', axh, 'LineWidth', p.Results.LineWidth);
                if td.nAxes > 1
                    TrialDataUtilities.Plotting.showInLegend(h, td.axisValueListsAsStringsShort{2}{iC});
                end
            end

            set(axh, 'XTick', xVals, 'XTickLabel', td.axisValueListsAsStringsShort{1});
            xlabel(axh, td.axisNames{1});
            ylabel(axh, 'Firing Rate (spikes/sec)');
            title(sprintf('Tuning Curve Unit %s', td.datasetName, unitName));

            au = AutoAxis.replace(axh);
            au.update();
        end

        function [timesCell, markCounts] = getMarkAlignedSpikeTimes(td, unitName, markIdx, window)
            % [timesCell] = getMarkAlignedSpikeTimes(td, unitName, markIdx, window)
            % timesCell is nTrials x nMarkOccurMax, only the first
            % markCounts(iTrial) are valid
            timesCell = td.getSpikeTimesUnaligned(unitName);
            timesCell = td.alignInfoActive.getMarkAlignedTimesCell(timesCell, markIdx, window);
            markCounts = td.alignInfoActive.markCountsValid;
        end

        function [timesCellByGroup, markCountsByGroup] = getMarkAlignedSpikeTimesGrouped(td, unitName, markIdx, window)
            % [timesCell] = getMarkAlignedSpikeTimes(td, unitName, markIdx, window)
            timesCell = td.getMarkAlignedSpikeTimes(unitName, markIdx, window);
            markCounts = td.alignInfoActive.markCountsValid;
            [timesCellByGroup, markCountsByGroup] = td.groupElementsFlat(timesCell, markCounts);
        end

        function [spikeCounts, markCounts] = getMarkAlignedSpikeCounts(td, unitName, markIdx, window)
            [timesCell, markCounts] = td.getMarkAlignedSpikeTimes(unitName, markIdx, window);
            spikeCounts = cellfun(@numel, timesCell);
            % NaN out spike counts for marks that didn't occur on a given trial
            for iTrial = 1:td.nTrials
                if isnan(markCounts(iTrial))
                    spikeCounts(iTrial, :) = NaN;
                else
                    spikeCounts(iTrial, markCounts(iTrial)+1:end) = NaN;
                end
            end
        end

        function [spikeCountsByGroup, markCountsByGroup] = getMarkAlignedSpikeCountsGrouped(td, unitName, markIdx, window)
            [spikeCounts, markCounts] = td.getMarkAlignedSpikeCounts(unitName, markIdx, window);
            [spikeCountsByGroup, markCountsByGroup] = td.groupElementsFlat(spikeCounts, markCounts);
        end

        function td = trimSpikeChannelToCurrentAlign(td, unitNames)
            % Timepoints that lie outside of TrialStart and TrialStop will
            % never be accessible via getTimes since they will be filtered
            % out by the AlignInfo

            td.warnIfNoArgOut(nargout);
            % default is TrialStart and TrialEnd, so just pass it along
            [startTimes, stopTimes] = td.getTimeStartStopEachTrial();
            td = td.trimSpikeChannel(unitNames, startTimes, stopTimes);
        end
    end

    methods % spike modification
        function td = setSpikeChannelWithinAlignWindow(td, ch, newTimes, varargin)
            p = inputParser();
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsInMemoryScale', false, @islogical); % if true, treat the data in values as memory class and scaling, so that it can be stored in .data as is
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            cd = td.getChannelDescriptor(ch);
            assert(isa(cd, 'SpikeChannelDescriptor') || isa(cd, 'SpikeArrayChannelDescriptor'));

            hasWaves = ~isempty(p.Results.waveforms) && cd.hasWaveforms;
            if hasWaves
               [wavesRaw, ~, timesRaw] = td.getRawSpikeWaveforms(ch, 'applyScaling', ~p.Results.waveformsInMemoryScale); % have the scaling match the input

               [waves, times] = td.replaceDataWithinAlignWindow(wavesRaw, timesRaw, ...
                   p.Results.waveforms, newTimes, 'includePadding', p.Results.includePadding);
            else
                timesRaw = td.getRawSpikeTimes(ch);

                times = td.replaceDataWithinAlignWindow(timesRaw, timesRaw, ...
                   newTimes, newTimes, 'includePadding', p.Results.includePadding);
                waves = [];
            end

            td = td.setSpikeChannel(ch, times, 'waveforms', waves, ...
                'waveformsInMemoryScale', p.Results.waveformsInMemoryScale, 'isAligned', false);
        end

        function td = maskSpikeChannelSpikes(td, ch, mask, varargin)
            p = inputParser();
            p.addParameter('keepRemovedSpikesAs', '', @isstringlike);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);
            keepAs = p.Results.keepRemovedSpikesAs;
            assert(ismatrix(mask) && size(mask, 1) == td.nTrials);

            cd = td.channelDescriptorsByName.(ch);
            assert(isa(cd, 'SpikeChannelDescriptor') || isa(cd, 'SpikeArrayChannelDescriptor'));

            timesOrig = td.getSpikeTimes(ch);
            times = timesOrig;
            notEmpty = ~cellfun(@isempty, times);
            times(notEmpty) = cellfun(@(t, m) t(m), times(notEmpty), mask(notEmpty), 'UniformOutput', false);

            hasWaves = cd.hasWaveforms;
            if hasWaves
                % maintain in memory scale to save time
                [wavesOrig, waveTvec] = td.getSpikeWaveforms(ch, 'applyScaling', false);
                waves = wavesOrig;
                waves(notEmpty) = cellfun(@(w, m) w(m, :), waves(notEmpty), mask(notEmpty), 'UniformOutput', false);
            else
                waves = {};
                waveTvec = [];
            end

            if ~isempty(p.Results.keepRemovedSpikesAs)
                timesRem = timesOrig;
                timesRem(notEmpty) = cellfun(@(t, m) t(~m), timesRem(notEmpty), mask(notEmpty), 'UniformOutput', false);
                if hasWaves
                    wavesRem = wavesOrig;
                    wavesRem(notEmpty) = cellfun(@(w, m) w(~m, :), wavesRem(notEmpty), mask(notEmpty), 'UniformOutput', false);
                else
                    wavesRem = {};
                end

                notEmpty = ~cellfun(@isempty, timesRem);
                if any(notEmpty)
                    if td.hasSpikeChannel(keepAs)
                        td = td.appendSpikesToSpikeChannel(keepAs, timesRem, 'waveforms', wavesRem, 'waveformsInMemoryScale', true);
                    elseif isa(cd, 'SpikeChannelDescriptor')
                        td = td.addSpikeChannel(keepAs, timesRem, 'waveforms', wavesRem, 'waveformsTime', waveTvec, 'waveformsInMemoryScale', true);
                    else
                        td = td.addSpikeArrayChannel(keepAs, timesRem, 'waveforms', wavesRem, 'waveformsTime', waveTvec, 'waveformsInMemoryScale', true);
                    end
                end
            end

            td = td.setSpikeChannelWithinAlignWindow(ch, times, 'waveforms', waves, ...
                'waveformsInMemoryScale', true);
        end

        function td = appendSpikesToSpikeChannel(td, ch, newTimes, varargin)
            p = inputParser();
            p.addParameter('waveforms', [], @iscell);
            p.addParameter('waveformsInMemoryScale', false, @islogical);
            p.parse(varargin{:});

            td.warnIfNoArgOut(nargout);

            td.assertHasSpikeChannel(ch);
            cd = td.channelDescriptorsByName.(ch);

            times = td.getSpikeTimes(ch);
            times = cellfun(@(t1, t2) sort(cat(1, t1, t2)), times, newTimes, 'UniformOutput', false);

            hasWaves = cd.hasWaveforms;
            if hasWaves
                waves = td.getSpikeWaveforms(ch);
                assert(~isempty(p.Results.waveforms), 'Must provide waveforms to append when existing unit has waveforms');
                waves = cellfun(@(w1, w2) cat(1, w1, w2), waves, p.Results.waveforms, 'UniformOutput', false);
            else
                waves = {};
            end

            td = td.setSpikeChannel(ch, times, 'isAligned', true, 'waveforms', waves, ...
                'waveformsInMemoryScale', p.Results.waveformsInMemoryScale);
        end
    end

    methods % spike waveforms
        function [wavesCell, waveTvec, timesCell] = getSpikeWaveforms(td, unitName, varargin)
            [wavesCell, waveTvec, timesCell] = getSpikeWaveforms@TrialData(td, unitName, varargin{:});
            [timesCell, maskCell] = td.alignInfoActive.getAlignedTimesCell(timesCell, true, 'singleTimepointTolerance', 0);
            wavesCell = cellfun(@(waves, mask) waves(mask, :), wavesCell, maskCell, 'UniformOutput', false);
            %timesCell = cellfun(@(times, mask) times(mask), timesCell, maskCell, 'UniformOutput', false);
        end

        function [wavesCell, waveTvec, timesCell] = getSpikeWaveformsGrouped(td, unitName, varargin)
            [wavesCell, waveTvec, timesCell] = td.getSpikeWaveforms(unitName, varargin{:});
            [wavesCell, timesCell] = td.groupElementsFlat(wavesCell, timesCell);
        end

        function plotSpikeWaveforms(td, unitName, varargin)
            p = inputParser();
            p.addParameter('maxToPlot', 500, @isscalar);
            p.addParameter('alpha', 0.2, @isscalar);
            p.addParameter('color', [], @(x) true);
            p.addParameter('colormap', [], @(x) isa(x, 'function_handle') || ismatrix(x) || isempty(x));
            p.addParameter('colorByCondition', false, @islogical); % if false, color by unit, if true color by condition
            p.addParameter('showThreshold', false, @islogical);
            p.addParameter('showMean', false, @islogical);
            p.addParameter('clickable', false, @islogical); % add interactive clicking to identify waveforms
            p.addParameter('fast', false, @islogical);
            p.parse(varargin{:});
            clickable = p.Results.clickable;

            if ~iscell(unitName)
                unitName = {unitName};
            end
            nUnits = numel(unitName);

            if ~p.Results.colorByCondition

                if ~isempty(p.Results.color)
                    colormap = AppearanceSpec.convertColor(p.Results.color);
                elseif isa(p.Results.colormap, 'function_handle')
                    colormap = p.Results.colormap(nUnits);
                elseif ~isempty(p.Results.colormap) && ismatrix(p.Results.colormap)
                    colormap = repmat(p.Results.colormap, ceil(nUnits / size(p.Results.colormap, 1)), 1);
                    %             elseif nUnits == 1
                    %                 colormap = [0 0 0];
                else
                    %                 colormap = distinguishable_colors(nUnits, {'w', 'k'});
                    % cbrewer set 1 with gray removed
                    if nUnits < 8
                        colormap = [0.894 0.102 0.11;0.216 0.494 0.722;0.302 0.686 0.29;0.596 0.306 0.639;1 0.498 0;1 1 0.2;0.651 0.337 0.157;0.969 0.506 0.749];
                    else
                        colormap = distinguishable_colors(nUnits, {'w'});
                    end
                end

                colormap = TrialDataUtilities.Plotting.expandWrapColormap(colormap, nUnits);

                nPlotGroups = nUnits;
                [wavesByTrialUnit, waveTvec, waveTimeCell] = td.getSpikeWaveforms(unitName, 'combine', false);
                groupNames = unitName;

                wavesCell = cellvec(nUnits);
                for iU = 1:nUnits
                    wavesCell{iU} = cat(1, wavesByTrialUnit{:, iU});
                end
                waveTvec = repmat({waveTvec}, nUnits, 1);

            else
                % plot colored by condition
                colormap = td.conditionColors;
                nPlotGroups = td.nConditions;

                [wavesGrouped, waveTvec, waveTimeCell] = td.getSpikeWaveformsGrouped(unitName, 'combine', true);
                wavesCell = cellfun(@(x) cat(1, x{:}), wavesGrouped, 'UniformOutput', false);
                waveTvec = repmat({waveTvec}, nPlotGroups, 1);
                groupNames = td.conditionNamesShort;
            end

            cdList = td.getChannelDescriptor(unitName);
            waveUnits = cdList(1).waveformsUnits;

            hMean = TrialDataUtilities.Plotting.allocateGraphicsHandleVector(nPlotGroups);

            for iU = 1:nPlotGroups
                wavesMat = wavesCell{iU};
                % wavesmat is nSpikes x nSamples;

                if clickable
                    % keep track of the details on each waveform
                    [timeInTrial, whichTrial] = TensorUtils.catWhich(1, waveTimeCell{:, iU});
                end

                if isempty(wavesMat)
                    continue;
                end

                maxToPlot = p.Results.maxToPlot;
                if maxToPlot < size(wavesMat, 1)
                    s = RandStream('mt19937ar','Seed',1);
                    idx = randsample(s, size(wavesMat, 1), maxToPlot);
                    wavesMat = wavesMat(idx, :);
                    if clickable
                        timeInTrial = timeInTrial(idx);
                        whichTrial = whichTrial(idx);
                    end
                end

                % should now be handled in TrialData getRawSpikeWaveforms
%                 if size(wavesMat, 2) < numel(waveTvec)
%                     waveTvec = waveTvec(1:size(wavesMat, 2));
%                 end

%                 h = TrialDataUtilities.Plotting.patchline(waveTvec, wavesMat', ...
%                     'Parent', axh, 'EdgeColor', colormap(iU, :), 'EdgeAlpha', p.Results.alpha);

                h = plot(waveTvec{iU}, wavesMat', 'Color', colormap(iU, :));
                TrialDataUtilities.Plotting.setLineOpacity(h, p.Results.alpha);
                TrialDataUtilities.Plotting.showFirstInLegend(h, groupNames{iU});

                if clickable
                    % generate description for each waveform
                    waveDesc = cellvec(size(wavesMat, 1));
                    if p.Results.colorByCondition
                        for iW = 1:size(wavesMat, 1)
                            waveDesc{iW} = sprintf('Unit %s\nTrial %d\nTime %s', unitName{iU}, whichTrial(iW), ...
                                td.alignInfoActive.buildStringForOffsetFromZero(timeInTrial(iW)));
                        end
                    else
                        for iW = 1:size(wavesMat, 1)
                            waveDesc{iW} = sprintf('Condition %s\nTrial %d\nTime %s', groupNames{iU}, whichTrial(iW), ...
                                td.alignInfoActive.buildStringForOffsetFromZero(timeInTrial(iW)));
                        end
                    end
                    TrialDataUtilities.Plotting.makeClickableShowDescription(h, waveDesc);
                end

                hold on;

                if p.Results.showMean
                    hMean(iU) = plot(waveTvec{iU}, mean(wavesMat, 1, 'omitnan'), '-', ...
                        'Color', colormap(iU, :), 'LineWidth', 2);
                end
            end
            if p.Results.showMean
                uistack(hMean, 'top');
            end

            if p.Results.showThreshold
                hold on
                threshEst = TrialDataUtilities.SpikeData.estimateThresholdFromSpikeWaveforms(td, unitName);
                TrialDataUtilities.Plotting.horizontalLine(threshEst, 'Color', 'r');
            end

            ht = title(sprintf('Spike Waveforms for %s', TrialDataUtilities.String.strjoin(unitName, ', ')));
            set(ht, 'Interpreter', 'none');

            axis tight;

            if ~p.Results.fast
                AutoAxis.replaceScaleBars(gca, td.timeUnitName, waveUnits);
            else
                box off;
            end
            hold off;
        end

        function plotSpikeWaveformsWithOtherUnits(td, unitName, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroUnit', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            allUnits = td.listSpikeChannelsOnSameArrayElectrodeAs(unitName, p.Results);
            otherUnits = setdiff(allUnits, unitName);
            unitNames = cat(1, otherUnits, unitName);
            cmap = [0.2*ones(numel(otherUnits), 3); 1, 0, 0];

            td.plotSpikeWaveforms(unitNames, 'colormap', cmap, p.Unmatched);
        end

        function plotSpikeWaveformsOnSameArrayElectrodeAs(td, chName, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroUnit', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unitNames = td.listSpikeChannelsOnSameArrayElectrodeAs(chName, p.Results);
            td.plotSpikeWaveforms(unitNames, p.Unmatched);
        end

        function plotSpikeWaveformsOnArrayElectrode(td, array, electrode, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroUnit', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            unitNames = td.listSpikeChannelsOnArrayElectrode(array, electrode, p.Results);
            td.plotSpikeWaveforms(unitNames, p.Unmatched);
        end

        function units = listSpikeChannelsMatchingRegex(td, units)
            if ischar(units)
                units = {units};
                wasChar = true;
            else
                wasChar = false;
            end

            allUnits = td.listSpikeChannels();
            matches = cellvec(numel(units));
            for iU = 1:numel(units)
                if ischar(units{iU})
                    units{iU} = units(iU);
                end
                matchInner = cellvec(numel(units{iU}));
                for iJ = 1:numel(units{iU})
                    matchInner{iJ} = ~cellfun(@isempty, regexp(allUnits, units{iU}{iJ}, 'match'));
                end
                 matches{iU} = allUnits(any(cat(2, matchInner{:}), 2));
            end

            if wasChar
                units = cat(1, matches{:});
            end
        end

        function units = listSpikeChannelsMatchingWildcard(td, units)
            if ischar(units)
                units = {units};
                wasChar = true;
            else
                wasChar = false;
            end

            allUnits = td.listSpikeChannels();
            matches = cellvec(numel(units));
            for iU = 1:numel(units)
                if ischar(units{iU})
                    units{iU} = units(iU);
                end
                matchInner = cellvec(numel(units{iU}));
                for iJ = 1:numel(units{iU})
                    regex = regexptranslate('wildcard', units{iU}{iJ});
                    matchInner{iJ} = ~cellfun(@isempty, regexp(allUnits, regex, 'match'));
                end
                 matches{iU} = allUnits(any(cat(2, matchInner{:}), 2));
            end

            if wasChar
                units = cat(1, matches{:});
            end
        end

        function [waveTvec, indStartByUnit, indStopByUnit] = getSpikeWaveformCommonTime(td, units)
            % for spike channels in cellstr units, find a common time
            % vector by aligning each units waveform tvec at 0

            if ischar(units)
                % easy single unit case
                waveTvec = makecol(td.channelDescriptorsByName.(units).waveformsTime);
                indStartByUnit = 1;
                indStopByUnit = numel(waveTvec);

            else
                [nPre, nPost, delta] = nanvec(numel(units));
                for iU = 1:numel(units)
                    tvec = td.channelDescriptorsByName.(units{iU}).waveformsTime;
                    [~, idxMin] = min(abs(tvec));
                    nPre(iU) = idxMin - 1;
                    nPost(iU) = numel(tvec) - idxMin;
                    delta(iU) = median(diff(tvec));
                end

                assert(max(delta) - min(delta) < 0.01 * median(delta), 'Waveform time vectors differ in their time deltas');


                waveTvec = (-max(nPre) : max(nPost)) * median(delta);
                indStartByUnit = max(nPre) - nPre + 1;
                indStopByUnit = numel(waveTvec) - (max(nPost) - nPost);
            end
        end

        function [wavesMat, waveTvec, timeWithinTrial, trialIdx, whichUnit] = getSpikeWaveformMatrix(td, units, varargin)
            % [wavesMat, waveTvec, timeWithinTrial, trialIdx, whichUnit] = getSpikeWaveformMatrix(td, units, varargin)
            % returns a matrix containing all waveforms in units. Units
            % can be a string or cellstr, where each string is a spike unit name
            % or, if paramValue 'regexp' is true, a regexp search over
            % unitNames
            p = inputParser();
            p.addParameter('regexp', false, @islogical); % if false, standard wildcard * search is done via regexptranslate

            % subselecting waveforms
            p.addParameter('maxWaves', Inf, @isscalar); % per unit
            p.addParameter('seed', 0, @isscalar);

            % alignment / interpolation. If align is true, interpolation
            % will be performed
            p.addParameter('align', false, @islogical);
            p.addParameter('alignMethod', 'upwardDownward', @ischar);
            p.addParameter('interpTroughThresh', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('interp', false, @islogical);
            p.addParameter('interpFactor', 10, @isscalar);

            % cleaning
            p.addParameter('clean', false, @islogical);
            p.addParameter('meanSDThresh', 0.4, @isscalar);
            p.addParameter('maxSDThresh', 2, @isscalar);
            p.parse(varargin{:});

            if ischar(units), units = {units}; end

            % do regexp matching
            if p.Results.regexp
                units = td.listSpikeChannelsMatchingRegex(units);
            else
                units = td.listSpikeChannels(units);
            end

            % flatten search results
            units = cat(1, units{:});

            [waveTvec, indStart, indStop] = td.getSpikeWaveformCommonTime(units);

            s = RandStream('mt19937ar','Seed',p.Results.seed);

            % get waveforms for each unit, do alignment / interpolation and
            % cleaning within each unit separately
            [waveByU, trialByU, timeByU] = deal(cellvec(numel(units)));
            for iU = 1:numel(units)
                [w, ~, t] = td.getSpikeWaveforms(units{iU});
                [waveByU{iU}, trialByU{iU}] = TensorUtils.catWhich(1, w{:});
                timeByU{iU} = cat(1, t{:});

                if size(waveByU{iU}, 1) > p.Results.maxWaves
                    mask = randsample(s, size(waveByU{iU}, 1), p.Results.maxWaves);
                    waveByU{iU} = waveByU{iU}(mask, :);
                    trialByU{iU} = trialByU{iU}(mask);
                    timeByU{iU} = timeByU{iU}(mask);
                end

                if p.Results.align
                    % will also do interpolation
                    [waveByU{iU}, waveTvec] = TrialDataUtilities.MKsort.alignSpline(waveByU{iU}, waveTvec, ...
                        'alignMethod', p.Results.alignMethod, 'interpFactor', p.Results.interpFactor, ...
                        'interpTroughThresh', p.Results.interpTroughThresh);
                elseif p.Results.interp
                    % just do interpolation
                    [waveByU{iU}, waveTvec] = interpWaves(waveByU{iU}, waveTvec, p.Results.interpFactor);
                end

                if p.Results.clean
                    [waveByU{iU}, maskKeep] = TrialDataUtilities.MKsort.cleanWaveforms(waveByU{iU}, ...
                        'meanSDThresh', p.Results.meanSDThresh, 'maxSDThresh', p.Results.maxSDThresh);
                    if ~isempty(maskKeep) && ~any(maskKeep)
                        warning('No waveforms kept by cleaning procedure');
                    end
                    trialByU{iU} = trialByU{iU}(maskKeep);
                    timeByU{iU} = timeByU{iU}(maskKeep);
                end
            end

            % concatenate waveforms together
            nWavesByUnit = cellfun(@(mat) size(mat, 1), waveByU);
            wavesMat = nan(sum(nWavesByUnit), numel(waveTvec));
            whichUnit = nan(sum(nWavesByUnit), 1);
            offset = 0;
            for iU = 1:numel(waveByU)
                wavesMat(offset + (1:nWavesByUnit(iU)), indStart(iU):indStop(iU)) = waveByU{iU};
                whichUnit(offset + (1:nWavesByUnit(iU))) = iU;
                offset = offset + nWavesByUnit(iU);
            end
            trialIdx = cat(1, trialByU{:});
            timeWithinTrial = cat(1, timeByU{:});

            function [iwaves, itvec] = interpWaves(waves, tvec, factor)
                itvec = linspace(min(tvec), max(tvec), factor*numel(tvec))';
                iwaves = zeros(size(waves, 1), numel(itvec));
                for iW = 1:size(waves, 1)
                    iwaves(iW, :) = interp1(tvec, waves(iW,:), itvec, 'spline');
                end
            end
        end

        function [waveMean, waveTvec, waveStd] = getSpikeWaveformMean(td, units, varargin)
            [wavesMat, waveTvec, ~, ~, whichUnit] = getSpikeWaveformMatrix(td, units, varargin{:});

            nU = numel(unique(whichUnit));

            [waveMean, waveStd] = deal(nan(nU, size(wavesMat, 2)));
            for iU = 1:nU
                waveMean(iU, :) = mean(wavesMat(whichUnit == iU, :), 'omitnan');
                waveStd(iU, :) = std(wavesMat(whichUnit == iU, :), 0, 'omitnan');
            end
        end
    end

    % Plotting Spike data
    methods
        function retInfo = plotPSTH(td, unitNames, varargin)
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('quick', false, @islogical);
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial fraction to average
            p.addParameter('assumePoissonStatistics', false, @islogical);
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.addParameter('errorType', '', @(s) ismember(s, {'sem', 'std', ''}));
            p.addParameter('showSem', true, @islogical); % equivalent to 'errorType', 'sem'
            p.addParameter('showRandomizedQuantiles', [], @(x) isempty(x) || isvector(x));

            p.addParameter('subtractConditionDescriptor', [], @(x) isempty(x) || isa(x, 'ConditionDescriptor'));

            p.addParameter('removeZeroSpikeTrials', false, @islogical);
            p.addParameter('useNativeScaling', false, @islogical);
            %p.addParameter('axisMarginLeft', 2.5, @isscalar);
            p.addParameter('axh', gca, @ishandle);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            retInfo = struct();
            alignIdx = p.Results.alignIdx;
            nAlignUsed = numel(alignIdx);

            sf = p.Results.spikeFilter;

            axh = TrialDataConditionAlign.getRequestedPlotAxis('axh', p.Results.axh);

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat, timeMask] = deal(cell(nAlignUsed, 1));
            for iAlign = 1:nAlignUsed
                [meanMat{iAlign}, tvecCell{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(alignIdx(iAlign)).getSpikeRateFilteredGroupMeans(unitNames, ...
                    'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                    'assumePoissonStatistics', p.Results.assumePoissonStatistics, ...
                    'spikeFilter', sf, ...
                    'combine', true, ...
                    'removeZeroSpikeTrials', p.Results.removeZeroSpikeTrials, ...
                    'subtractConditionDescriptor', p.Results.subtractConditionDescriptor, ...
                    'useNativeScaling', p.Results.useNativeScaling);

                if numel(tvecCell{iAlign}) > 1
                    [tvecTemp, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}, timeMask{iAlign}] = ...
                        TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}', meanMat{iAlign}, ...
                        semMat{iAlign}, stdMat{iAlign});
                    tvecCell{iAlign} = tvecTemp';
                end
            end

            if isempty(p.Results.errorType)
                if p.Results.showSem
                    errorMat = semMat;
                else
                    errorMat = [];
                end
            else
                switch p.Results.errorType
                    case ''
                        errorMat = [];
                    case 'sem'
                        errorMat = semMat;
                    case 'std'
                        errorMat = stdMat;
                end
            end

            if ~isempty(p.Results.showRandomizedQuantiles)
                quantileData = cell(nAlignUsed, 1);
                for iAlign = 1:nAlignUsed
                    quantileData{iAlign} = td.useAlign(alignIdx(iAlign)).getSpikeRateFilteredGroupMeansRandomizedQuantiles(unitNames, ...
                        'quantiles', p.Results.showRandomizedQuantiles, ...
                        'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                        'spikeFilter', sf, ...
                        'removeZeroSpikeTrials', p.Results.removeZeroSpikeTrials);

                    % apply same time slicing
                    quantileData{iAlign} = quantileData{iAlign}(:, timeMask{iAlign}, :);
                end
            else
                quantileData = [];
            end

            maskEmpty = cellfun(@isempty, tvecCell);
            if any(maskEmpty)
                error('No valid time window found for alignment. Perhaps there are not enough trials with spikes in all of the conditions? Try lowering minTrials or setting removeZeroSpikeTrials to false.');
            end

            retInfo = td.plotProvidedAnalogDataGroupMeans(1, 'time', tvecCell, ...
                'data', meanMat, 'dataError', errorMat, 'quantileData', quantileData, 'axh', axh, ...
                'axisInfoX', 'time', 'axisInfoY', struct('name', 'Firing Rate', 'units', 'spikes/sec'), ...
                'alignIdx', alignIdx, ...
                'quick', p.Results.quick, 'binAlignmentMode', sf.binAlignmentMode, 'binWidth', sf.timeDelta, 'retInfo', retInfo, p.Unmatched);

            if isnumeric(unitNames)
                unitNames = td.lookupSpikeChannelByIndex(unitNames);
            end
            TrialDataUtilities.Plotting.setTitleIfBlank(axh, '%s : Unit %s', td.datasetName, TrialDataUtilities.String.strjoin(unitNames, ', '));
            axis(axh, 'tight');

            ylabel('spikes/sec');
            if p.Results.quick
                xlabel('time (ms)');
            else
                au = AutoAxis(axh);
%                 au.addAutoAxisY();
                %au.axisMarginLeft = p.Results.axisMarginLeft;
                au.update();
            end

            hold(axh, 'off');
        end

        function plotRaster(td, unitNames, varargin)

            % black + cbrewer Qual Set1
            def_cmap = [0 0 0; 0.894 0.102 0.11;0.216 0.494 0.722;0.302 0.686 0.29;0.596 0.306 0.639;1 0.498 0;1 1 0.2;0.651 0.337 0.157;0.969 0.506 0.749;0.6 0.6 0.6];

            p = inputParser();
            p.addParameter('conditionIdx', 1:td.nConditions, @isvector);
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('spikeColorMap', def_cmap, @(x) true);
%             p.addParameter('spikeColor', 'k', @(x) true);
            p.addParameter('colorSpikesLikeCondition', false, @islogical);
            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.addParameter('timeAxisTickBridgeExtendToLimits', true, @islogical);
            p.addParameter('tickHeight', 1, @isscalar);
            p.addParameter('tickWidth', 1, @isscalar);
            p.addParameter('tickAlpha', 1, @isscalar);

            p.addParameter('markAlpha', 0.5, @isscalar);
            p.addParameter('markTickWidth', 2, @isscalar);
            p.addParameter('showRanges', true, @islogical);
            p.addParameter('includePadding', false, @islogical); % mostly used for debugging PSTH computation

            p.addParameter('intervalAlpha', 0.5, @isscalar);
            p.addParameter('intervalMinWidth', NaN, @isscalar); % if specified, draws intervals at least this wide to ensure visibility
            p.addParameter('gapBetweenConditions', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('plotConditionDividers', false, @islogical);

            p.addParameter('removeZeroSpikeTrials', false, @islogical);

            % for plotting interval and marks above each condition
            p.addParameter('annotateAboveEachCondition', false, @islogical);
            p.addParameter('annotationHeight', 2, @isscalar);
            p.addParameter('annotateUsingFirstTrialEachCondition', true, @islogical);

            % make room for labels using AutoAxis
            p.addParameter('axisMarginLeft', [], @isscalar);

            p.addParameter('useShortLabels', true, @islogical);

            % shade start:stop intervals in gray to show valid time interval
            p.addParameter('shadeStartStopInterval', false, @islogical);
            p.addParameter('shadeOutsideStartStopInterval', false, @islogical);
            p.addParameter('shadeBlankingRegions', false, @islogical);
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('intervalShowOnData', false, @islogical);

            % if true, draw spike waveforms instead of ticks
            p.addParameter('drawSpikeWaveforms', false, @islogical);
            p.addParameter('spikeWaveformScaleHeight', 1, @isscalar);
            p.addParameter('spikeWaveformScaleTime', 1, @isscalar);

            p.addParameter('interUnitOffset', 0, @isscalar); % add a small offset between multiply plotted units to detect simultaneous spikes

            p.addParameter('quick', false, @islogical);

            p.addParameter('axh', [], @(x) true);

            p.addParameter('combine', false, @islogical); % combine units as -one, if false, plot spikes in different colorstds
            
            p.addParameter('precomputedData', struct(), @isstruct); % special input that allows passing precomputed outputs of getSpikeTimesGrouped to save time when cycling through many units
            
%             p.KeepUnmatched = true;
            p.parse(varargin{:});
            precomputedData = p.Results.precomputedData;
            includePadding = p.Results.includePadding;

            if td.nTrialsValid == 0
                error('No valid trials found');
            end

            % for multiple unit simultaneous plotting, can be useful for
            % looking at spike sorting issues
            if isnumeric(unitNames)
                unitNames = td.lookupSpikeChannelByIndex(unitNames);
            end
            unitNames = string(unitNames);

            if p.Results.combine
                nUnits = 1;
                % wrap in another cell so that loop works
            else
                nUnits = td.computeSpikeUnitCountFromName(unitNames);
            end

            axh = td.getRequestedPlotAxis('axh', p.Results.axh);

            conditionIdx = p.Results.conditionIdx;
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            alignIdx = p.Results.alignIdx;
            nAlignUsed = numel(alignIdx);

            times = cell(nConditionsUsed, nUnits);
            if p.Results.drawSpikeWaveforms
                waves = cell(nConditionsUsed, nUnits);
            end

            % compute x-axis offsets for each align
            timePointsCell = cell(nAlignUsed, 1);

            if nUnits == 0
                % special mode with no units just for visualizing marks and rasters
                [startData, stopData] = deal(cell(nAlignUsed, nConditionsUsed, 1));
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);

                    if ~includePadding
                        [start, stop] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrial();
                    else
                        [start, stop] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrialWithPadding();
                    end
                    mask = start <= stop;
                    timePointsCell{iAlign} = [start(mask); stop(mask)];

                    % and store start/stop by trial in each align/condition
                    % cell
                    for iCond = 1:nConditionsUsed
                        idxCond = conditionIdx(iCond);
                        startData{iAlign, iCond} = start(td.listByCondition{idxCond});
                        stopData{iAlign, iCond} = stop(td.listByCondition{idxCond});
                    end
                end
            else
                % regular mode where we grab spikes
                % get grouped spike times by alignment
                if isfield(precomputedData, 'spikeTimesGroupedByConditionUnit')
                    times = precomputedData.spikeTimesGroupedByConditionUnit;
                    times = times(conditionIdx, :);
                else
                    times = td.getSpikeTimesGroupedByConditionUnit(unitNames, 'combineAligns', true, 'alignIdx', alignIdx, 'combine', p.Results.combine, 'includePadding', includePadding); % nCond {nTrials x nUnits}
                    times = times(conditionIdx, :);
                end

                if p.Results.drawSpikeWaveforms
                    error('This needs to be fixed for new combineAlign access pattern');
%                     [wavesC, wavesTvec] = td.useAlign(idxAlign).getSpikeWaveformsGrouped(unitNames, 'combine', p.Results.combine);
%                     idxCond = conditionIdx(iC);
%                     for iU = 1:nUnits
%                         waves{iC, iU} = wavesC{idxCond}(:, iU);
%                     end
                end

                [startData, stopData] = deal(cell(nAlignUsed, nConditionsUsed));
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);

                    % figure out time validity window for this alignment
                    % TODO might want to update this for the selected
                    % conditions only
                    if ~includePadding
                        [start, stop] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrial();
                    else
                        [start, stop] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrialWithPadding();
                    end
                    mask = start <= stop;
                    timePointsCell{iAlign} = [start(mask); stop(mask)];

                    % and store start/stop by trial in each align/condition
                    % cell
                    for iCond = 1:nConditionsUsed
                        idxCond = conditionIdx(iCond);
                        startData{iAlign, iCond} = start(td.listByCondition{idxCond});
                        stopData{iAlign, iCond} = stop(td.listByCondition{idxCond});
                    end
                end
            end
            [tOffsetByAlign, tLimitsByAlign] = td.getAlignPlottingTimeOffsets([], 'includePadding', p.Results.includePadding);

            % optionally filter out trials that have zero spikes
            listByCondition = td.listByCondition(conditionIdx);

            if p.Results.removeZeroSpikeTrials && nUnits > 0
                % C x U
                hasSpikesData = cellfun(@(thisC) ~cellfun(@isempty, thisC), times, 'UniformOutput', false);
                % C x 1
                hasSpikesAnyUnit = arrayfun(@(iC) TensorUtils.anyMultiDim(cat(2, hasSpikesData{iC, :}), 2), (1:nConditionsUsed)', 'UniformOutput', false);

                listByConditionMask = hasSpikesAnyUnit;
            else
                listByConditionMask = cellfun(@(x) true(size(x)), listByCondition, 'UniformOutput', false);
            end

            trialCounts = cellfun(@nnz, listByConditionMask);
            nTrialsTotal = sum(trialCounts);

            % default gap depends on presence of annotations between conditions
            % it is set such that 10% of the space is split among the gaps
            % equally
            if isempty(p.Results.gapBetweenConditions)
                gap = nTrialsTotal * 0.1 / nConditionsUsed;
                if p.Results.annotateAboveEachCondition
                    gap = gap + p.Results.annotationHeight + 2; % 1 above, 1 below (plus annotation height)
                    dividerOffset = p.Results.annotationHeight + 1.5;
                else
                    dividerOffset = gap / 2;
                end

            else
                gap = p.Results.gapBetweenConditions;
                if p.Results.annotateAboveEachCondition
                    dividerOffset = p.Results.annotationHeight + 1 + gap/2;
                else
                    dividerOffset = gap/2;
                end
            end

            % compute y-axis offsets for each condition
            % note that trials proceed from the top downwards, and the
            % offsets here refer to the top of the first trial's spikes

            yOffsetByCondition = zeros(nConditionsUsed, 1);
            yLimsByCondition = nan(2, nConditionsUsed);
            currentOffset = 0;
%             lastTrialCount = 0;
            for iC = nConditionsUsed:-1:1
                % top of this block of trials occurs here
                yOffsetByCondition(iC) = currentOffset + trialCounts(iC);
                % min max (for last, first trials)
                yLimsByCondition(:, iC) = [currentOffset; currentOffset + trialCounts(iC)];
                if trialCounts(iC) > 0
                    currentOffset = currentOffset + trialCounts(iC) + gap;
%                     lastTrialCount = trialCounts(iC);
                end
            end
            yCentersByCondition = mean(yLimsByCondition, 1)';

            yDividersByCondition = yOffsetByCondition(2:end) + dividerOffset;

            % if we're drawing waveforms, normalize all waveforms to the [0 1] range]
            if p.Results.drawSpikeWaveforms
                [maxW, minW] = deal(nan(nConditionsUsed, nUnits));
                for iU = 1:nUnits
                    for iC = 1:nConditionsUsed
                        waves = cat(1, waves{iC, iU}{listByConditionMask{iC}});
                        if isempty(waves)
                            continue;
                        end
                        maxW(iC, iU) = max(waves(:), [], 'omitnan');
                        minW(iC, iU) = min(waves(:), [], 'omitnan');
                    end
                end

                maxW = max(maxW(:), [], 'omitnan');
                minW = min(minW(:), [], 'omitnan');

                waves = cellfun(@(wc) cellfun(@(w) (w - minW) / (maxW - minW), wc, 'UniformOutput', false), waves, 'UniformOutput', false);
            end

            % draw tick rasters in a grid pattern (conditions vertically,
            % aligns horizontally)

            offsetByUnit = (0:nUnits-1) * p.Results.interUnitOffset;
            for iU = 1:nUnits
                for iC = 1:nConditionsUsed
                    app = td.conditionAppearances(conditionIdx(iC));
                    if p.Results.colorSpikesLikeCondition
                        color = app.Color;
                    else
                        idx = mod(iU-1, size(p.Results.spikeColorMap, 1))+1;
                        color = p.Results.spikeColorMap(idx, :);
                    end

                    if p.Results.drawSpikeWaveforms && ~p.Results.quick
                        % draw waveforms in lieu of ticks
                            TrialDataUtilities.Plotting.drawTickRaster(times{iC, iU}(listByConditionMask{iC}), ...
                            'xOffset', offsetByUnit(iU), 'yOffset', yOffsetByCondition(iC), ...
                            'color', color, 'alpha', p.Results.tickAlpha, 'axh', axh, ...
                            'waveCell', waves{iC, iU}(listByConditionMask{iC}), 'waveformTimeRelative', wavesTvec, ...
                            'normalizeWaveforms', false, ... % already normalized to [0 1]
                            'waveScaleHeight', p.Results.spikeWaveformScaleHeight, 'waveScaleTime', p.Results.spikeWaveformScaleTime);
                    else
                        % draw vertical ticks
                        TrialDataUtilities.Plotting.drawTickRaster(times{iC, iU}(listByConditionMask{iC}), ...
                            'xOffset', offsetByUnit(iU), 'yOffset', yOffsetByCondition(iC), ...
                            'color', color, 'alpha', p.Results.tickAlpha, 'lineWidth', p.Results.tickWidth, ...
                            'tickHeight', p.Results.tickHeight, 'axh', axh);
                    end
                    hold(axh, 'on');
                end
            end

            % if we're shading blanking regions, we create a temporary interval to hold them
            if p.Results.shadeBlankingRegions
                blankRegions = td.getSpikeBlankingRegions(unitNames, 'combine', true);
                app = AppearanceSpec('Color', [0.7 0.7 0.7]);
                td = td.intervalManual(blankRegions, 'Blanking Regions', ...
                    'appear', app, 'showOnData', true, 'showOnAxis', false);
            end
            
            if p.Results.plotConditionDividers
                tLims = [min(tLimitsByAlign(:)); max(tLimitsByAlign(:))];
                for iC = 1:nConditionsUsed
                    app = td.conditionAppearances(conditionIdx(iC));
                    line(axh, tLims, yLimsByCondition([1 1], iC), 'Color', app.Color);
                    line(axh, tLims, yLimsByCondition([2 2], iC), 'Color', app.Color);
                end
            end

            if ~p.Results.quick
                % draw marks and intervals on each raster
                if p.Results.annotateAboveEachCondition
                    if p.Results.annotateUsingFirstTrialEachCondition
                        % draw first trial's marks and intervals above the
                        % spikes for that condition. First trial is used when
                        % only the gist of the marks/intervals is required and
                        % is better seen from an individual trial rather than
                        % the average
                        for iAlign = 1:nAlignUsed
                            idxAlign = alignIdx(iAlign);
                            for iCond = 1:nConditionsUsed
                                if isempty(listByCondition{iCond}(listByConditionMask{iCond}))
                                    continue;
                                end
                                % use first trial to draw marks and intervals
                                % on
                                firstTrialIdx = find(listByConditionMask{iCond}, 1);
                                td.alignInfoSet{idxAlign}.drawOnRasterByTrial(...
                                    'showMarks', p.Results.markShowOnData, ...
                                    'showIntervals', p.Results.intervalShowOnData, ...
                                    'startByTrial', startData{iAlign, iCond}(firstTrialIdx), ...
                                    'stopByTrial', stopData{iAlign, iCond}(firstTrialIdx), ...
                                    'trialIdx', listByCondition{iCond}(firstTrialIdx), ...
                                    'showInLegend', iCond == 1, 'tOffsetZero', tOffsetByAlign(iAlign), ...
                                    'yOffsetTop', yOffsetByCondition(iCond) + p.Results.annotationHeight + 1, ...
                                    'tickHeight', p.Results.annotationHeight, ...
                                    'intervalMinWidth', p.Results.intervalMinWidth, ...
                                    'axh', axh, 'intervalAlpha', p.Results.intervalAlpha, ...
                                    'markAlpha', p.Results.markAlpha, 'markTickWidth', p.Results.markTickWidth);
                            end
                        end
                    else
                        error('Not yet implemented');
                    end
                else
                    % draw marks and intervals on each trial
                    for iAlign = 1:nAlignUsed
                        idxAlign = alignIdx(iAlign);
                        for iCond = 1:nConditionsUsed
                            if isempty(td.listByCondition{conditionIdx(iCond)}), continue; end
                            td.alignInfoSet{idxAlign}.drawOnRasterByTrial(...
                                'showMarks', p.Results.markShowOnData, ...
                                'showIntervals', p.Results.intervalShowOnData, ...
                                'startByTrial', startData{iAlign, iCond}(listByConditionMask{iCond}), ...
                                'stopByTrial', stopData{iAlign, iCond}(listByConditionMask{iCond}), ...
                                'trialIdx', listByCondition{iCond}(listByConditionMask{iCond}), ...
                                'showInLegend', iCond == 1, 'tOffsetZero', tOffsetByAlign(iAlign), ...
                                'yOffsetTop', yOffsetByCondition(iCond), ...
                                'tickHeight', p.Results.tickHeight, ...
                                'intervalMinWidth', p.Results.intervalMinWidth, ...
                                'axh', axh, 'intervalAlpha', p.Results.intervalAlpha, ...
                                'shadeStartStopInterval', p.Results.shadeStartStopInterval, ...
                                'shadeOutsideStartStopInterval', p.Results.shadeOutsideStartStopInterval, ...
                                'fullTimeLimits', tLimitsByAlign(iAlign, :), ...
                                'markAlpha', p.Results.markAlpha, 'markTickWidth', p.Results.markTickWidth);
                        end
                    end
                end
            end

            yLims = [min(yLimsByCondition(:)), max(yLimsByCondition(:))];
            if p.Results.annotateAboveEachCondition
                yLims(2) = yLims(2) + p.Results.annotationHeight + 1;
            end

            % setup y axis condition labels
            if td.isGrouped
                colors = cat(1, td.conditionAppearances(conditionIdx).Color);
                if p.Results.useShortLabels
                    conditionNames = td.conditionNamesShort(conditionIdx); % todo change this to multiline?
                else
                    conditionNames = td.conditionNamesMultiline(conditionIdx);
                end

                if ~p.Results.quick
                    % only include conditions with at least 1 trial
                    mask = trialCounts(conditionIdx) > 0;
                    au = AutoAxis(axh);
                    au.addLabeledSpan('y', 'span', yLimsByCondition(:, mask), 'label', ...
                        conditionNames(mask), 'color', colors(mask, :));

                    mask = trialCounts(conditionIdx(2:end)) > 0;
                    set(axh, 'YTick', flipud(yDividersByCondition(mask)));
                else
                    mask = trialCounts(conditionIdx) > 0;
                    
%                     yTickLabels = flipud(conditionNames(mask));
                    yTickLabelsColored = strings(numel(conditionNames), 1);
                    for iC = 1:numel(conditionNames)
                        yTickLabelsColored(iC) = sprintf("\\color[rgb]{%.3f,%.3f,%.3f}%s", colors(iC, :), conditionNames{iC});
                    end
   
                    set(axh, 'YTick', flipud(yCentersByCondition(mask)), 'YTickLabel', flipud(yTickLabelsColored(mask)));
                end
            else
                set(axh, 'YTick', []);
            end

            % setup time axis markers
            if p.Results.quick
                set(axh, 'TickDir', 'out');
                td.alignSummarySet{1}.setupTimeAutoAxis('which', 'x', ...
                    'style', 'quick', 'axh', axh);
            % elseif p.Results.annotateAboveEachCondition
            %     if iAlign == nAlignUsed
            %         % just use horizontal scale bar
            %         au = AutoAxis(axh);
            %         au.xUnits = td.timeUnitName;
            %         au.addAutoScaleBarX();
            %         au.update();
            %         au.installCallbacks();
            %     end
            else
                % use marks or tickBridges via AlignSummary
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    td.alignSummarySet{idxAlign}.setupTimeAutoAxis('axh', axh, 'which', 'x', ...
                        'style', p.Results.timeAxisStyle, 'tOffsetZero', tOffsetByAlign(iAlign), ...
                        'showRanges', p.Results.showRanges, 'tickBridgeExtendToLimits', p.Results.timeAxisTickBridgeExtendToLimits);
                end
            end

            if nUnits > 0
                if isstring(unitNames)
                    unitNameStr = TrialDataUtilities.String.strjoin(unitNames, ',');
                else
                    unitNameStr = unitName;
                end
%                 TrialDataUtilities.Plotting.setTitleIfBlank(axh, '%s %s', td.datasetName, unitNameStr);
                title(axh, sprintf('%s %s', td.datasetName, unitNameStr), 'Interpreter', 'none');
            end

            tLims = [min(tLimitsByAlign(:)), max(tLimitsByAlign(:))];
            if tLims(2) > tLims(1)
                set(axh, 'XLim', tLims);
            end
            if yLims(2) > yLims(1)
                set(axh, 'YLim', yLims);
            end

            axisMarginLeft = p.Results.axisMarginLeft;
            if isempty(axisMarginLeft)
                scale = getFigureSizeScale();
                if td.isGrouped
                    axisMarginLeft = 1.5 * scale;
                else
                    axisMarginLeft = 0.5 * scale;
                end
            end
            if ~p.Results.quick
                au = AutoAxis(axh);
                au.axisMarginLeft = axisMarginLeft; % make room for left hand side labels
                axis(axh, 'off');
                au.update();
            else
                box(axh, 'off');
            end
            hold(axh, 'off');
        end
        
        function plotMarkRaster(td, varargin)
             td.plotRaster([], 'markShowOnData', true, 'intervalShowOnData', true, 'shadeOutsideStartStopInterval', true, varargin{:});
        end
    end

    % Plotting Analog each trial
    methods
        function aa = setupTimeAxis(td, varargin)
            p = inputParser();
            p.addParameter('trialIdx', [], @(x) isempty(x) || isvector(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            trialIdx = p.Results.trialIdx;
            if isempty(trialIdx)
                alignSummaryActive = td.alignSummaryActive;
            else
                alignSummarySet = td.buildAlignSummarySetWithTrials(trialIdx);
                alignSummaryActive = alignSummarySet{td.alignInfoActiveIdx};
            end

            aa = alignSummaryActive.setupTimeAutoAxis('style', 'tickBridge', p.Unmatched);
        end

        function [offsets, lims] = getAlignPlottingTimeOffsets(td, tvecCell, varargin)
            % when plotting multiple alignments of data simultaneously,
            % timeseries are plotted side by side along an axis, separated
            % by gaps given by .interAlignGaps. This returns the list of
            % time offsets where each alignment's zero event would be
            % plotted. The time vectors in tvec Cell are used to determine
            % the "width" of each alignments and to find the appropriate
            % location to start the next alignment.
            %
            % tvecCell is a nAlign x 1 cell of:
            %     time vectors, for common time across trials, or
            %     nTrials x 1 cell of time vectors, for different times
            %     across trials.
            %
            % Parameters:
            %     alignIdx : vector indicating which alignments to include

            p = inputParser();
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('includePadding', false, @islogical);
            p.addParameter('alignTimeOffsets', [], @(x) isempty(x) || isvector(x));
            p.parse(varargin{:});
            
            % select alignment indices
            alignIdx = 1:td.nAlign;
            alignIdx = alignIdx(p.Results.alignIdx);
            nAlign = numel(alignIdx);

            includePadding = p.Results.includePadding;

            % compute start/stop of each alignment
            [mins, maxs] = deal(nanvec(nAlign));
            if isempty(tvecCell)
                % use min / max per align (full valid window)
                for iAlign = 1:nAlign
                    if ~includePadding
                        [mins_this, maxs_this] = td.alignInfoSet{alignIdx(iAlign)}.getStartStopRelativeToZeroByTrial();
                    else
                        [mins_this, maxs_this] = td.alignInfoSet{alignIdx(iAlign)}.getStartStopRelativeToZeroByTrialWithPadding();
                    end
                    mins(iAlign) = min(mins_this, [], 'omitnan');
                    maxs(iAlign) = max(maxs_this, [], 'omitnan');
                end
            else
                if ~iscell(tvecCell)
                    tvecCell = {tvecCell};
                end
            
                for iAlign = 1:nAlign
                    if iscell(tvecCell{iAlign}) && ~isempty(tvecCell{iAlign})
                        mins(iAlign) = min(cellfun(@nanminNanEmpty, tvecCell{iAlign}), [], 'omitnan');
                        maxs(iAlign) = max(cellfun(@nanmaxNanEmpty, tvecCell{iAlign}), [], 'omitnan');
                    elseif isempty(tvecCell{iAlign})
                        error('Time cell for alignment %d is empty', iAlign);
                    else
                        mins(iAlign) = min(tvecCell{iAlign}, [], 'omitnan');
                        maxs(iAlign) = max(tvecCell{iAlign}, [], 'omitnan');
                    end
                end
            end

            if isempty(td.manualAlignTimeOffsets) && isempty(p.Results.alignTimeOffsets)
                offsets = nan(nAlign, 1);
                offsets(1) = 0;
                currentOffset = 0;

                if nAlign == 1
                    % no gaps needed
                    offsets = 0;
                    lims = [mins(1), maxs(1)];
                    return;
                end

                % determine the inter alignment gaps
                if isempty(td.interAlignGaps)
                    % no inter align gap specified, determine automatically as
                    % 2% of total span
                    T = sum(maxs - mins);
                    gaps = repmat(0.02 * T, td.nAlign-1, 1);
                else
                    gaps = td.interAlignGaps;
                end

                for iAlign = 2:nAlign
                    currentOffset = currentOffset + maxs(iAlign-1) + ...
                        gaps(iAlign-1) - mins(iAlign);
                    offsets(iAlign) = currentOffset;
                end
            else
                if ~isempty(p.Results.alignTimeOffsets)
                    offsets = p.Results.alignTimeOffsets;
                    if numel(offsets) == nAlign - 1
                        offsets = [0; offsets];
                    end
                    offsets = offsets(alignIdx);
                else
                    offsets = td.manualAlignTimeOffsets(alignIdx);
                end
            end

            lims = [mins + offsets, maxs + offsets];

            function r = nanmaxNanEmpty(v1)
                if isempty(v1)
                    r = NaN;
                else
                    r = max(v1, [], 'omitnan');
                end
            end
        end
        
        function [tvec, whichAlign] = catTimeVecOverAlignWithSeparator(td, tvecCell, varargin)
            p = inputParser();
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.parse(varargin{:});
            alignIdx = p.Results.alignIdx;
            assert(numel(tvecCell) == numel(alignIdx));
            
            nPerAlign = cellfun(@numel, tvecCell);
            nTotal = sum(nPerAlign) + numel(alignIdx) - 1;
            
            tOffsets = td.getAlignPlottingTimeOffsets(tvecCell, 'alignIdx', p.Results.alignIdx);
            
            tvec = nan(nTotal, 1);
            whichAlign = nan(nTotal, 1);
            offset = 1; 
            for iAlign = 1:numel(tvecCell)
                this = tvecCell{iAlign};
                insert_at = offset + (0:numel(this)-1);
                tvec(insert_at) = this + tOffsets(iAlign);
                whichAlign(insert_at) = alignIdx(iAlign);
                offset = offset + numel(this) + 1;
            end 
        end
        
        function data = catDataOverAlignWithSeparator(td, time_dim, dataCell) %#ok<INUSL>
            sz = size(dataCell{1});
            
            separator_sz = sz;
            separator_sz(time_dim) = 1;
            separator = nan(separator_sz, 'like', dataCell{1});
            
            catArgs = cell(numel(dataCell)*2 - 1, 1);
            catArgs(1:2:end) = dataCell;
            [catArgs{2:2:end}] = deal(separator);
            
            data = cat(time_dim, catArgs{:});
        end

        function plotProvidedAnalogDataGroupedEachTrial(td, D, varargin)
            % common utility function for drawing analog data grouped by
            % condition, used by the plotAnalog functions below.
            % Can handle multiple alignments. If a single alignment is
            % specified, defaults to the active alignment. Override by
            % specifying parameter 'alignIdx'. If nAlign alignments are
            % specified, uses all of them.
            %
            % Assumes that data are grouped according to the current
            % grouping of td, that is, the configuration of the data
            % matches the current state of td.
            %
            % time is either:
            %     common time vector
            %     cell with size (nConditions or 1) x nAlign cell containing either
            %       nTrials x 1 cell of time vectors or
            %       time vector for all trials in that condition / alignment.
            %
            % for D==1: 1-D data vs. time as matrix
            %   data is nConditions x nAlign cell containing nTrials x 1
            %     cell of T x D matrices, or is nTrials x T x D data matrix
            % for D==2 or D==3: 2-D or 3-D data:
            %   data are nConditions x nAlign cell containing nTrials x 1
            %     cell of T x D matrices or containing nTrials x T x D data matrices

            p = inputParser();
            p.addParameter('time', [], @(x) isvector(x) || iscell(x)); % for D == 1,2,3 (for marking)
            p.addParameter('data', {}, @iscell); % for D == 1,2,3
            p.addParameter('yOffsetBetweenTrials', 0, @isscalar);
            p.addParameter('yOffsetBetweenConditions', [], @isscalar);

            p.addParameter('axisInfoX', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoZ', [], @(x) isempty(x) || isstringlike(x) || isa(x, 'ChannelDescriptor'));

            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @isstringlike);
            p.addParameter('axisStyleY', 'tickBridge', @isstringlike);

            p.addParameter('conditionIdx', 1:td.nConditions, @isnumeric);
            p.addParameter('alignIdx', [], @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('lineWidth', NaN, @isscalar);
            p.addParameter('showRangesOnAxis', true, @islogical); % show ranges for marks below axis

            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('markShowOnAxis', true, @islogical);
            p.addParameter('markShowInLegend', false, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 8, @isscalar);
            %%
            p.addParameter('markOutline', true, @islogical);
            p.addParameter('markOutlineAlpha', 1, @isscalar);


            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalShowOnAxis', true, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);

            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('useTranslucentMark3d', false, @islogical);

            p.addParameter('sortOrderMode', 'byTrial', @isstringlike); % byTrial, byCondition

            p.addParameter('quick', false, @islogical);
            p.addParameter('clickable', false, @islogical); % make interactive and clickable
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            quick = p.Results.quick;

            if td.nTrialsValid == 0
                error('No valid trials to plot');
            end

            axh = td.getRequestedPlotAxis(p.Unmatched);

            conditionIdx = p.Results.conditionIdx;
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            app = td.conditionAppearances(conditionIdx);

            if ~ismember('scaleBars', p.UsingDefaults) && p.Results.scaleBars
                axisStyleX = 'scaleBar';
                axisStyleY = 'scaleBar';
            else
                axisStyleX = p.Results.axisStyleX;
                axisStyleY = p.Results.axisStyleY;
            end

            % splay analog traces out vertically for D==1
            yOffsetTrial = p.Results.yOffsetBetweenTrials;
            if isempty(p.Results.yOffsetBetweenConditions)
                yOffsetCondition = 5 * yOffsetTrial;
            else
                yOffsetCondition = p.Results.yOffsetBetweenConditions;
            end

            plotOptions = [{'Clipping', 'off'} p.Results.plotOptions{:}];

            if ~isnan(p.Results.lineWidth)
                % this will overwrite condition specific line width if specified
                plotOptions = [{'LineWidth', p.Results.lineWidth}, plotOptions{:}];
            end
            
            % based on dimensionality, check sizes of provided
            if D == 1
                data = p.Results.data;
                time = p.Results.time;
                if isempty(data) || isempty(time)
                    error('Must provide data, time for D==1');
                end

                assert(isvector(time) || iscell(time), 'Time cell must be vector or cell');
                assert(iscell(data) , 'Data cell must be be cell');

                % check size of data cell
                if size(data, 1) == nConditionsUsed
                    % okay as is
                elseif size(data, 1) == td.nConditions
                    % needs to be masekd
                    data = data(conditionIdx, :);
                end
                nAlignUsed = size(data, 2);

                assert(nAlignUsed == 1 || iscell(time), 'Time must be vector if nAlign==size(data, 2) > 1');

                % check length of time cell, should be nConditionsUsed x
                % nAlign at the end
                if iscell(time)
                    if size(time, 1) == 1
                        % clone across conditions
                        time = repmat(time, nConditionsUsed, 1);
                    elseif size(time, 1) == nConditionsUsed
                        % okay as is
                    elseif size(time, 1) == td.nConditions
                        % needs to be masekd
                        time = time(conditionIdx, :);
                    end

                    assert(size(time, 2) == nAlignUsed, 'size(time, 2) must be nAlign');
                else
                    % convert to cell
                    time = repmat({time}, nConditionsUsed, nAlignUsed);
                end

             elseif D == 2 || D == 3
                data = p.Results.data;
                time = p.Results.time;
                if isempty(data) || isempty(time)
                    error('Must provide data, time for D==2, 3');
                end

                nAlignUsed = size(data, 2);

                % check length of data cell
                if size(data,1) == nConditionsUsed
                    % okay as is
                elseif size(data,1) == td.nConditions
                    % needs to be masekd
                    data = data(conditionIdx, :);
                end

                % check length of time cell, should be nConditionsUsed x
                % nAlign at the end
                assert(nAlignUsed == 1 || iscell(time), 'Time must be vector if nAlign==size(data, 2) > 1');

                if iscell(time)
                    if size(time, 1) == 1
                        % clone across conditions
                        time = repmat(time, nConditionsUsed, 1);
                    elseif size(time, 1) == nConditionsUsed
                        % okay as is
                    elseif size(time, 1) == td.nConditions
                        % needs to be masekd
                        time = time(conditionIdx, :);
                    end

                    assert(size(time, 2) == nAlignUsed, 'size(time, 2) must be nAlign');
                else
                    % convert to cell
                    time = repmat({time}, nConditionsUsed, nAlignUsed);
                end

                assert(size(data, 2) == nAlignUsed, 'size(dataY, 2) must match nAlign==size(dataX, 2)');
            else
                error('D must be 1,2,3');
            end

            % figure out which alignments are used
            if isempty(p.Results.alignIdx)
                if nAlignUsed == 1
                    alignIdx = td.alignInfoActiveIdx;
                else
                    alignIdx = 1:td.nAlign;
                end
            else
                alignIdx = 1:td.nAlign;
                alignIdx = alignIdx(p.Results.alignIdx);
            end

            % compute time offsets between successive alignments
            if D == 1
                timeByAlign = cell(numel(alignIdx), 1);
                for iA = 1:numel(alignIdx)
                    timeByAlign{iA} = cat(1, time{:, iA});
                end
                timeOffsetByAlign = td.getAlignPlottingTimeOffsets(timeByAlign, 'alignIdx', alignIdx);
            else
                timeOffsetByAlign = zeros(nAlignUsed, 1);
            end

            % store handles as we go
            hData = cell(nConditionsUsed, nAlignUsed);

            yOffsetsByCondition = cell(nConditionsUsed, 1);

            % plot data traces
            for iAlign = 1:nAlignUsed
                yOffsetCurrent = 0;

                % plot one condition at a time
                for iCond = 1:nConditionsUsed
%                     zlevel = iCond / (2*nConditionsUsed);

                    % determine whether we can plot all trials simultaneously
                    % (matrix format with shared time vector) or one at a time
                    % (cell element for each trial). Grab the data for this
                    % condition.
                    plotSimultaneously = false;
                    if D == 1
                        timeC = time{iCond, iAlign};
                        dataC = data{iCond, iAlign};
                        if iscell(dataC)
                            nTrialsC = numel(dataC);
                        else
                            nTrialsC = size(dataC, 1);
                            plotSimultaneously = true;
                        end

                        if nTrialsC > 0
                            yOffsets = makecol(yOffsetCurrent + (0:(nTrialsC-1)) * yOffsetTrial);
                            yOffsetCurrent = max(yOffsets) + yOffsetCondition;
                            yOffsetsByCondition{iCond} = yOffsets;
                        end

                    elseif D == 2 || D == 3
                        dataC = data{iCond, iAlign};

                        if iscell(dataC)
                            nTrialsC = numel(dataC);
                        else
                            nTrialsC = size(dataC, 1);
                            plotSimultaneously = true;
                        end
                    end

                    if plotSimultaneously
                        % plot all trials from this condition simultaneously
                        % timeC, dataC will be nTrials x T x D
                        hold(axh, 'on');
                        if D == 1
                            tOffset = timeOffsetByAlign(iAlign);
                            timeC = timeC + tOffset;
                            dataC = bsxfun(@plus, dataC, yOffsets);
%                             if p.Results.alpha < 1
%                                hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(timeC', dataC', ...
%                                    ones(size(timeC)) * zlevel, ...
%                                    'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
%                                    'LineWidth', app(iCond).LineWidth, plotOptions{:});
%                             else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, timeC', dataC', '-', ...
                                    plotArgs{:}, plotOptions{:});
                                if p.Results.alpha < 1
                                    TrialDataUtilities.Plotting.setLineOpacity(hData{iCond, iAlign}, p.Results.alpha);
                                end
%                             end

                        elseif D == 2
%                             if p.Results.alpha < 1
%                                hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dataC(:, :, 1)', dataC(:, :, 2)', ...
%                                    ones(size(dataC(:, :, 1)')) * zlevel, ...
%                                    'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
%                                    'LineWidth', app(iCond).LineWidth, plotOptions{:});
%                             else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, dataC(:, :, 1)', dataC(:, :, 2)', '-', ...
                                    plotArgs{:}, plotOptions{:});
                                if p.Results.alpha < 1
                                    TrialDataUtilities.Plotting.setLineOpacity(hData{iCond, iAlign}, p.Results.alpha);
                                end
%                             end

                        elseif D == 3
%                             if p.Results.alpha < 1
%                                 hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dataC(:, :, 1)', dataC(:, :, 2)', dataC(:, :, 3)', ...
%                                    'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
%                                    'LineWidth', app(iCond).LineWidth, 'Parent', axh, plotOptions{:});
%                             else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot3(axh, dataC(:, :, 1)', dataC(:, :, 2)', dataC(:, :, 3)', '-', ...
                                    plotArgs{:}, plotOptions{:});
                                if p.Results.alpha < 1
                                    TrialDataUtilities.Plotting.setLineOpacity(hData{iCond, iAlign}, p.Results.alpha);
                                end
%                             end
                        end
                    else
                        % plot each trial from this condition individually
                        hData{iCond, iAlign} = gobjects(nTrialsC, 1);
                        for iTrial = 1:nTrialsC
                            if D == 1
                                tvec = timeC{iTrial};
                                dvec = dataC{iTrial};

                                tOffset = timeOffsetByAlign(iAlign);
                                yOffset = yOffsets(iTrial);

                                if isstringlike(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                                    xvec = dvec;
                                    yvec = tvec;
                                else
                                    yvec = dvec;
                                    xvec = tvec;
                                end

                                if ~isempty(tvec) && ~isempty(dvec)
%                                     if p.Results.alpha < 1
%                                        hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(xvec + tOffset, yvec + yOffset, ...
%                                            ones(size(yvec)) * zlevel, ...
%                                            'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
%                                            'LineWidth', app(iCond).LineWidth, plotOptions{:});
%                                     else
                                        plotArgs = app(iCond).getPlotArgs();
                                        hData{iCond, iAlign}(iTrial) = plot(axh, xvec + tOffset, yvec + yOffset, '-', ...
                                            plotArgs{:}, plotOptions{:});
                                        if p.Results.alpha < 1
                                            TrialDataUtilities.Plotting.setLineOpacity(hData{iCond, iAlign}(iTrial), p.Results.alpha);
                                        end
%                                     end
                                end

                            elseif D==2
                                assert(size(dataC{iTrial}, 2) == 2, 'Expecting data matrix per trial with two columns corresponding to dim 1 and 2');
                                dxvec = dataC{iTrial}(:, 1);
                                dyvec = dataC{iTrial}(:, 2);

                                if ~isempty(dxvec) && ~isempty(dyvec)
%                                     if p.Results.alpha < 1
%                                        
%                                        hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(dxvec, dyvec, ...
%                                            ones(size(dxvec)) * zlevel, ...
%                                            'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
%                                            'LineWidth', app(iCond).LineWidth, plotOptions{:});
%                                     else
                                    plotArgs = app(iCond).getPlotArgs();
                                    hData{iCond, iAlign}(iTrial) = plot(axh, dxvec, dyvec, '-', ...
                                        plotArgs{:}, plotOptions{:});
                                    hData{iCond, iAlign}(iTrial).Color(4) = p.Results.alpha;
%                                     if p.Results.alpha < 1
%                                         TrialDataUtilities.Plotting.setLineOpacity(hData{iCond, iAlign}(iTrial), p.Results.alpha);
%                                     end
                                end

                            elseif D==3
                                dxvec = dataC{iTrial}(:, 1);
                                dyvec = dataC{iTrial}(:, 2);
                                dzvec = dataC{iTrial}(:, 3);

                                if ~isempty(dxvec) && ~isempty(dyvec) && ~isempty(dzvec)
% %                                     if p.Results.alpha < 1
% %                                        hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline3(dxvec, dyvec, dzvec, ...
% %                                            'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
% %                                            'LineWidth', app(iCond).LineWidth, plotOptions{:});
% %                                     else
                                    plotArgs = app(iCond).getPlotArgs();
                                    hData{iCond, iAlign}(iTrial) = plot3(axh, dxvec, dyvec, dzvec, '-', ...
                                        plotArgs{:}, plotOptions{:});
                                    hData{iCond, iAlign}(iTrial).Color(4) = p.Results.alpha;
%                                     end
                                end
                            end

                            if iTrial == 1
                                hold(axh, 'on');
                            end
                        end
                    end

                    % update name for inclusion in legend
                    TrialDataUtilities.Plotting.showFirstInLegend(hData{iCond, iAlign}, td.conditionNamesShort{iCond});
                end
            end

            switch p.Results.sortOrderMode
                case 'byCondition'
                    % resort the handles by condition
                    hDataVec = cat(1, hData{:});
                    TrialDataUtilities.Plotting.graphicsSortInPlaceAs(hDataVec);

                case 'byTrial'
                    % resort the handles by trial number
                    hDataVec = cat(1, hData{:});
                    trialIdxVec = cat(1, td.listByCondition{conditionIdx});

                    [~, sortIdx] = sort(trialIdxVec);
                    hDataVec = hDataVec(sortIdx);
                    TrialDataUtilities.Plotting.graphicsSortInPlaceAs(hDataVec);

                otherwise
                    error('Unknown sortOrderMode %s. Valid options are byCondition and byTrial.', p.Results.sortOrderMode);
            end

            if ~quick
                % plot marks and intervals on top of data
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    for iCond = 1:nConditionsUsed
                        if isempty(td.listByCondition{conditionIdx(iCond)}), continue; end
                        timeC = time{iCond, iAlign};
                        dataC = data{iCond, iAlign};

                        if D==1
                            if iscell(dataC)
                                for iTrial = 1:numel(dataC)
                                    dataC{iTrial} = dataC{iTrial} + yOffsetsByCondition{iCond}(iTrial);
                                end
                            else
                                dataC = bsxfun(@plus, dataC, yOffsetsByCondition{iCond});
                            end
                        end

                        if D == 3
                            markIntervalClipping = 'off';
                        else
                            markIntervalClipping = 'on';
                        end

                        td.alignInfoSet{idxAlign}.drawOnDataByTrial('time', timeC, 'data', dataC, ...
                            'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                            'showInLegend', p.Results.markShowInLegend && iCond == 1, 'tOffsetZero', timeOffsetByAlign(iAlign), ...
                            'axh', axh, 'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize, ...
                            'clipping', markIntervalClipping, ...
                            'intervalAlpha', p.Results.intervalAlpha, ...
                            'showMarks', p.Results.markShowOnData, ...
                            'markOutline', p.Results.markOutline, ...
                            'markOutlineAlpha', p.Results.markOutlineAlpha, ...
                            'showIntervals', p.Results.intervalShowOnData);
                    end

                    % setup time axis for this align
                    if D == 1
                        td.alignSummarySet{idxAlign}.setupTimeAutoAxis('which', 'x', 'style', p.Results.timeAxisStyle, ...
                            'tOffsetZero', timeOffsetByAlign(iAlign), 'showMarks', p.Results.markShowOnAxis, 'showRanges', p.Results.showRangesOnAxis);
                    end
                end
            end

            % setup non-time axes
            if quick
                axisStyleX = 'label';
                axisStyleY = 'label';
%                 axisStyleZ = 'label';
            end
            if D == 1
                if isstringlike(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                    % x is data, y is time
                    if ~isempty(p.Results.axisInfoX)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                            'which', 'x', 'style', axisStyleX);
                    end
                else
                    % y is data, x is time
                    if ~isempty(p.Results.axisInfoY)
%                         if yOffsetTrial ~= 0 || yOffsetCondition ~= 0
%                             scaleBar = true;
%                         else
%                             scaleBar = false;
%                         end
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                            'which', 'y', 'style', axisStyleY);
                    end
                end

            elseif D == 2
                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, 'which', 'x', ...
                        'style', axisStyleX);
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, 'which', 'y', ...
                        'style', axisStyleY);
                end

            elseif D == 3
                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, 'which', 'x', 'useAutoAxis', false);
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, 'which', 'y', 'useAutoAxis', false);
                end
                if ~isempty(p.Results.axisInfoZ)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoZ, 'which', 'z', 'useAutoAxis', false);
                end

                if p.Results.useThreeVector
                    ThreeVector(axh);
                    axis(axh, 'off');
                    axis(axh, 'vis3d');
                end
            end

            % make plot interactive if requested
            if p.Results.clickable
                descGrouped = td.getTrialDescriptionsGrouped('multiline', true);
                descGrouped = descGrouped(conditionIdx);

                for iAlign = 1:nAlignUsed
                    for iCond = 1:nConditionsUsed
                        hThis = hData{iCond, iAlign};
                        descThis = descGrouped{iCond};
                        assert(numel(hThis) == numel(descThis), 'Internal issue with condition descriptions - count does not match');

                        for iH = 1:numel(hThis)
                            set(hThis(iH), 'Description', descThis{iH});
                        end
                    end
                end

                hSetFull = cat(1, hData{:});
                TrialDataUtilities.Plotting.makeClickableShowDescription(hSetFull);
            end

            box(axh, 'off');
            axis(axh, 'tight');

            if D < 3 && ~quick
                au = AutoAxis(axh);
                if strcmp(axisStyleY, 'scaleBar')
                    au.axisMarginLeft = 0.8;
                end
    %             if strcmp(axisStyleX, 'scaleBar')
    %                 au.axisMarginBottom = 1;
    %             end
                au.update();
                au.installCallbacks();
            end
            hold(axh, 'off');

            set(axh, 'SortMethod', 'childorder'); % fix rendering issues
        end

        function plotAnalogEachTrial(td, name, varargin)
            td.ungroupMaintainValid.plotAnalogGroupedEachTrial(name, varargin{:});
        end

        function plotAnalogGroupedEachTrial(td, name, varargin)
            p = inputParser;
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.addParameter('ensureUniformSampling', false, @islogical);
            p.addParameter('timeDelta', []);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);

            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % grab raw data (for marking) and grouped data (for plotting)
            if td.hasAnalogChannel(name) || td.hasAnalogChannelGroup(name)
                [dataByGroup, timeByGroup] = td.getAnalogGroupedEachAlign(name, p.Results);
                cd = td.getChannelDescriptor(name); % in case analogGroup(5)

            elseif td.hasSpikeChannel(name)
                [dataByGroup, timeByGroup] = td.getSpikeRateFilteredGrouped(name);
                cd = td.channelDescriptorsByName.(name);
            elseif td.hasChannel(name)
                error('Unknown channel type');
            else
                error('Channel %s not found', name);
            end

            td.plotProvidedAnalogDataGroupedEachTrial(1, 'time', timeByGroup, ...
                'data', dataByGroup, 'axisInfoY', cd, ...
                p.Unmatched);
        end

        function plotAnalogGroupedEachTrial2D(td, name1, name2, varargin)
            p = inputParser;
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.addParameter('timeDelta', []);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            names = [string(name1), string(name2)];
            cds = td.getChannelDescriptor(names);
            [dataCell, timeCell] = td.getAnalogMultiCommonTimeGrouped(names, p.Results);
            td.plotProvidedAnalogDataGroupedEachTrial(2, ...
                'time', timeCell(:), 'data', dataCell(:), ...
                'axisInfoX', cds(1), ...
                'axisInfoY', cds(2), p.Unmatched);
        end

        function plotAnalogGroupedEachTrial2DvsTime(td, name1, name2, varargin)
            p = inputParser;
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [dataCell, time] = td.getAnalogMultiCommonTimeGrouped({name1, name2}, p.Results);
            td.plotProvidedAnalogDataGroupedEachTrial(2, ...
                'time', time(:), 'data', dataCell(:), ...
                'axisInfoZ', 'time', ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), p.Unmatched);
        end

        function plotAnalogGroupedEachTrial3D(td, name1, name2, name3, varargin)
            p = inputParser;
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [dataCell, timeCell] = td.getAnalogMultiCommonTimeGrouped({name1, name2, name3}, p.Results);
            td.plotProvidedAnalogDataGroupedEachTrial(3, ...
                'time', timeCell(:), 'data', dataCell(:), ...
                'axisInfoX', td.getChannelDescriptor(name1), ...
                'axisInfoY', td.getChannelDescriptor(name2), ...
                'axisInfoZ', td.getChannelDescriptor(name3), p.Unmatched);
        end
    end

    % Plotting Analog means by group
    methods
        function addConditionLabels(td, varargin)
            td.conditionInfo.addColoredLabels(varargin{:});
        end

        function retInfo = plotProvidedAnalogDataGroupMeans(td, D, varargin)
            p = inputParser();
            p.addParameter('time', [], @(x) true);
            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.addParameter('alignIdx', 1:td.nAlign, @(x) true);
            p.addParameter('retInfo', struct(), @isstruct);
            p.addParameter('alignTimeOffsets', [], @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            alignTimeOffsets = td.getAlignPlottingTimeOffsets(p.Results.time, 'alignIdx', p.Results.alignIdx, 'alignTimeOffsets', p.Results.alignTimeOffsets);
            res = rmfield(p.Results, 'timeAxisStyle');
            res.xAxisStyle = p.Results.timeAxisStyle;
            
            retInfo = TrialDataConditionAlign.plotConditionAlignedAnalogDataGroupMeans(D, p.Results, p.Unmatched, ...
                'conditionDescriptor', td.conditionInfo, ...
                'alignSummarySet', td.alignSummarySet, ...
                'alignInfoActiveIdx', td.alignInfoActiveIdx, ...
                'alignTimeOffsets', alignTimeOffsets);


            retInfo.alignTimeOffsets = alignTimeOffsets;
        end

        function plotAnalogGroupMeans(td, name, varargin)
            % plot the mean and sem for an analog channel vs. time within
            % each condition

            p = inputParser();
            p.addParameter('alignIdx', 1:td.nAlign, @(x) true);
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp
            p.addParameter('interpolateMethod', 'linear', @isstringlike); % see interp1 for details
            p.addParameter('assumeUniformSampling', false, @islogical);

            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.addParameter('errorType', '', @(s) ismember(s, {'sem', 'std', ''}));
            p.addParameter('showSem', true, @islogical);
            p.addParameter('showRandomizedQuantiles', [], @(x) isempty(x) || isvector(x));
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionDescriptor', [], @(x) isempty(x) || isa(x, 'ConditionDescriptor'));

            p.addParameter('label', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            alignIdx = p.Results.alignIdx;
            nAlignUsed = numel(alignIdx);

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat, timeMask] = deal(cell(nAlignUsed, 1));
            for iAlign = 1:nAlignUsed
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(alignIdx(iAlign)).getAnalogGroupMeans(name, 'minTrials', p.Results.minTrials, ...
                    'minTrialFraction', p.Results.minTrialFraction, ...
                    'subtractTrialBaseline', p.Results.subtractTrialBaseline, ...
                    'subtractTrialBaselineAt', p.Results.subtractTrialBaselineAt, ...
                    'subtractConditionBaselineAt', p.Results.subtractConditionBaselineAt, ...
                    'subtractConditionDescriptor', p.Results.subtractConditionDescriptor, ...
                    'timeDelta', p.Results.timeDelta, ...
                    'timeReference', p.Results.timeReference, ...
                    'binAlignmentMode', p.Results.binAlignmentMode, ...
                    'resampleMethod', p.Results.resampleMethod, ...
                    'assumeUniformSampling', p.Results.assumeUniformSampling);
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}, timeMask{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, ...
                    semMat{iAlign}, stdMat{iAlign});
            end

            switch p.Results.errorType
                case ''
                    errorMat = [];
                case 'sem'
                    errorMat = semMat;
                case 'std'
                    errorMat = stdMat;
            end

            if p.Results.showSem && isempty(errorMat)
                errorMat = semMat;
            end


            if ~isempty(p.Results.showRandomizedQuantiles)
                quantileData = cell(nAlignUsed, 1);
                for iAlign = 1:nAlignUsed
                    quantileData{iAlign} = td.useAlign(alignIdx(iAlign)).getAnalogGroupMeansRandomizedQuantiles(name, ...
                        'quantiles', p.Results.showRandomizedQuantiles, ...
                         'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                        'subtractTrialBaseline', p.Results.subtractTrialBaseline, ...
                        'subtractTrialBaselineAt', p.Results.subtractTrialBaselineAt, ...
                        'subtractConditionBaselineAt', p.Results.subtractConditionBaselineAt);

                    % apply same time slicing
                    quantileData{iAlign} = quantileData{iAlign}(:, timeMask{iAlign}, :);
                end
            else
                quantileData = [];
            end

            cd = td.getChannelDescriptor(name);
            td.plotProvidedAnalogDataGroupMeans(1, 'time', tvecCell, ...
                'data', meanMat, 'dataError', errorMat, p.Unmatched, ...
                'quantileData', quantileData, ...
                'axisInfoX', 'time', 'axisInfoY', cd, 'labelY', p.Results.label, ...
                'alignIdx', alignIdx, ...
                p.Unmatched);
        end

        function plotAnalogGroupMeans2D(td, name1, name2, varargin)
            % plot the mean and sem for an analog channel vs. time within
            % each condition

            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getAnalogMultiGroupMeans({name1, name2}, ...
                    'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                    'subtractTrialBaseline', p.Results.subtractTrialBaseline, ...
                    'subtractTrialBaselineAt', p.Results.subtractTrialBaselineAt, ...
                    'subtractConditionBaselineAt', p.Results.subtractConditionBaselineAt);
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, ...
                    semMat{iAlign}, stdMat{iAlign});
            end
            td.plotProvidedAnalogDataGroupMeans(2, 'time', tvecCell, ...
                'data', meanMat, 'dataError', semMat, p.Unmatched, ...
                'axisInfoX', td.getChannelDescriptor(name1), ...
                'axisInfoY', td.getChannelDescriptor(name2), p.Unmatched);
        end

        function plotAnalogGroupMeans3D(td, name1, name2, name3, varargin)
            % plot the mean and sem for an analog channel vs. time within
            % each condition

            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('minTrialFraction', 0, @isscalar); % minimum fraction of trials required for average
            p.addParameter('subtractTrialBaseline', [], @(x) true);
            p.addParameter('subtractTrialBaselineAt', '', @isstringlike);
            p.addParameter('subtractConditionBaselineAt', '', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getAnalogMultiGroupMeans({name1, name2, name3}, ...
                    'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                    'subtractTrialBaseline', p.Results.subtractTrialBaseline, ...
                    'subtractTrialBaselineAt', p.Results.subtractTrialBaselineAt, ...
                    'subtractConditionBaselineAt', p.Results.subtractConditionBaselineAt);
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, ...
                    semMat{iAlign}, stdMat{iAlign});
            end

            td.plotProvidedAnalogDataGroupMeans(3, 'time', tvecCell, ...
                'data', meanMat, p.Unmatched, ...
                'axisInfoX', td.getChannelDescriptor(name1), ...
                'axisInfoY', td.getChannelDescriptor(name2), ...
                'axisInfoZ', td.getChannelDescriptor(name3));
        end

        function plotAnalogGroupMeansByConditionAxes(td, name, varargin)
            % plot the mean and sem for an analog channel vs. time within
            % each condition according to the arrangment of td.conditions
            % The first dimension will be plotted along each row, the
            % second along the columns, and subsequent dimensions within
            % the same axes on top of the others
            import TrialDataUtilities.Plotting.hideInLegend;
            import TrialDataUtilities.Plotting.showInLegend;
            [pan, args] = td.getRequestedPlotPanel(varargin{:});

            import TrialDataUtilities.Plotting.errorshade;
            p = inputParser();
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('legend', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(args);

            showLegend = p.Results.legend;

            % determine number of subplots needed
            szPanel = td.conditionInfo.conditionsSize(1:2);
            pan.pack(szPanel(1), szPanel(2));
            pan.units = 'cm';
            pan.margin = 0;
            pan.de.margin = 0;
            cInds = td.conditionInfo.conditionsAsLinearInds;
            app = td.conditionAppearances;

            % loop over alignments and gather time vectors
            [meanMat, semMat, tvecCell] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}] = td.useAlign(iAlign).getAnalogGroupMeans(name);
            end

            % get global time limits
            % and determine the x-offsets at which to plot the data
            [tOffsets, xlims] = td.getAlignPlottingTimeOffsets(tvecCell);

            for row = 1:szPanel(1)
                for col = 1:szPanel(2)
                    axh = pan(row, col).select();
                    au = AutoAxis();
                    au.axisInset = [2 0.4 0.4 0.4];

                    cIndsThisPanel = squeeze(cInds(row, col, :));

                    % track legend handles so that legend auto appears nicely
                    legHandleByCond = nan(numel(cIndsThisPanel), 1);

                    for iCond = 1:numel(cIndsThisPanel)
                        c = cIndsThisPanel(iCond);
                        for iAlign = 1:td.nAlign
                            h = errorshade(tvecCell{iAlign} + tOffsets(iAlign), ...
                                meanMat{iAlign}(c, :), ...
                                semMat{iAlign}(c, :), app(c).Color, 'axh', axh);
                            hold(axh, 'on');

                            % handle legend entries
                            if ~isnan(h)
                                hideInLegend(h);
                                set(h, 'DisplayName', td.conditionNames{c});
                                legHandleByCond(iCond) = h;
                            end

                            % handle draw on data for this condition
                        end

                    end

                    tdTheseCond = td.selectConditions(cIndsThisPanel);
                    for iAlign = 1:td.nAlign
                        tdTheseCond.alignSummarySet{iAlign}.setupTimeAutoAxis('tOffsetZero', tOffsets(iAlign), 'style', 'tickBridge');
                    end
                    au.axisInset(2) = 2;

                    showInLegend(removenan(legHandleByCond));
                    box(axh, 'off');
                    axis(axh, 'tight');
                    set(axh, 'XLim', xlims);
                    axis(axh, 'off');

                    % setup y axis
                    if col == 1 || true
                        TrialDataUtilities.Plotting.setupAxisForChannel(td.channelDescriptorsByName.(name));
                    end

                    if showLegend
                        legend(axh, 'show');
                        legend(axh, 'boxoff');
                    end
                end
            end
        end
    end

    % Plotting single trials
    methods
        function plotSingleTrialAnalogChannels(td, varargin)
            p = inputParser();
            p.addOptional('channels', td.listAnalogChannels(), @(x) isstringlike(x) || iscellstr(x) || isstring(x));
            p.addParameter('validTrialIdx', [], @isvector); % selection into valid trials
            p.addParameter('trialIdx', [], @isvector); % selection into all trials
            p.addParameter('normalize', true, @islogical);
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('quick', false, @islogical);
            p.addParameter('commonTime', false, @islogical);
            p.addParameter('channelLabels', {}, @iscellstr);
            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if isempty(p.Results.validTrialIdx)
                if isempty(p.Results.trialIdx)
                    % auto choose first valid trial
                    idx = find(td.valid, 1);
                else
                    idx = p.Results.trialIdx;
                end
            else
                idx = mi2ui(p.Results.validTrialIdx, td.valid);
            end

            chList = string(p.Results.channels);

            cdCell = td.getChannelDescriptorMulti(chList);
            if ~isempty(p.Results.channelLabels)
                channelLabels = p.Results.channelLabels;
            else
                channelLabels = chList;
            end
            dataUnits = arrayfun(@(cd) cd.unitsPrimary, cdCell, 'UniformOutput', false);

            td = td.selectTrials(idx);
            if p.Results.commonTime
                [data, time] = td.getAnalogMultiCommonTime(chList, 'timeDelta', p.Results.timeDelta);
                data = data{1};
                time = time{1};
            else
                [data, time] = td.getAnalogMulti(chList, 'timeDelta', p.Results.timeDelta);
                time = time';
                data = data';
            end

            TrialDataUtilities.Plotting.plotStackedTraces(time, data, 'labels', channelLabels, ...
                'showLabels', true, 'dataUnits', dataUnits, 'normalize', p.Results.normalize, 'quick', p.Results.quick, ...
                p.Unmatched);
            xlabel('');
            if ~p.Results.quick
                td.alignSummaryActive.setupTimeAutoAxis('style', p.Results.timeAxisStyle, 'labelFirstMarkOnly', true);

                ax = AutoAxis();
                ax.axisMarginLeft = 5;
                ax.update();
            end
        end

        function plotSingleTrialAnalogChannelGroup(td, groupName, varargin)
            td.plotSingleTrialAnalogChannels(td.listAnalogChannelsInGroup(groupName), varargin{:});
        end

        function plotSingleTrialRaster(td, varargin)
            p = inputParser();
            p.addParameter('trialIdx', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x))); % selection into all trials
            p.addParameter('validTrialIdx', [], @isvector); % selection into valid trials
            p.addParameter('arrayNames', [], @isstringlike);
            p.addParameter('unitNames', [], @(x) isstringlike(x));
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.addParameter('tickHeight', 1, @isscalar);
%             % if true, draw spike waveforms instead of ticks
            p.addParameter('drawSpikeWaveforms', false, @islogical);
            p.addParameter('spikeWaveformScaleHeight', 1, @isscalar);
            p.addParameter('spikeWaveformScaleTime', 1, @isscalar);

            p.addParameter('quick', false, @islogical);
            p.addParameter('axh', gca, @ishandle);
            p.parse(varargin{:});

            if td.nTrialsValid == 0
                error('No valid trials found');
            end

            if isempty(p.Results.validTrialIdx)
                if isempty(p.Results.trialIdx)
                    % auto choose first valid trial
                    idx = find(td.valid, 1);
                else
                    idx = p.Results.trialIdx;
                end
            else
                idx = mi2ui(p.Results.validTrialIdx, td.valid);
            end

            td = td.selectTrials(idx);

            % for multiple unit simultaneous plotting, can be useful for
            % looking at spike sorting issues
            unitNames = string(p.Results.unitNames);
            if isempty(unitNames)  
                unitNames = td.listSpikeChannels();
            end
            nUnits = numel(unitNames);

            axh = td.getRequestedPlotAxis('axh', p.Results.axh);

            alignIdx = p.Results.alignIdx;
            nAlignUsed = numel(alignIdx);

            %timesByAlign = cell(nAlignUsed, nUnits);
            if p.Results.drawSpikeWaveforms
                wavesByAlign = cell(nAlignUsed, nUnits);
            end

            % compute x-axis offsets for each align
            timePointsCell = cell(nAlignUsed, 1);

            [start, stop] = nanvec(nAlignUsed);
            for iAlign = 1:nAlignUsed
                idxAlign = alignIdx(iAlign);
                % figure out time validity window for this alignment
                % TODO might want to update this for the selected
                % conditions only
                [start(iAlign), stop(iAlign)] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrial();
                timePointsCell{iAlign} = [start(iAlign); stop(iAlign)];
            end

            % get grouped spike times by alignment
            timesByAlign = td.useAlign(idxAlign).getSpikeTimesEachAlign(unitNames);

%             if p.Results.drawSpikeWaveforms
%                 [waves, wavesTvec] = td.useAlign(idxAlign).getSpikeWaveformsEachAlign(unitName);
%             end
            tOffsetByAlign = zerosvec(nAlignUsed);
%             tOffsetByAlign = td.getAlignPlottingTimeOffsets(timePointsCell);

            % if we're drawing waveforms, normalize all waveforms to the [0 1] range]
%             if p.Results.drawSpikeWaveforms
%                 [maxW, minW] = deal(nan(nAlignUsed, nConditionsUsed, nUnits));
%                 for iU = 1:nUnits
%                     for iA = 1:nAlignUsed
%                             waves = cat(1, wavesByAlign{iA, iU});
%                             if isempty(waves)
%                                 continue;
%                             end
%                             maxW(iA, iC, iU) = nanmax(waves(:));
%                             minW(iA, iC, iU) = nanmin(waves(:));
%                     end
%                 end
%
%                 maxW = nanmax(maxW(:));
%                 minW = nanmin(minW(:));
%
%                 wavesByAlign = cellfun(@(wc) cellfun(@(w) (w - minW) / (maxW - minW), wc, 'UniformOutput', false), wavesByAlign, 'UniformOutput', false);
%             end

            % draw tick rasters in a grid pattern (conditions vertically,
            % aligns horizontally)
            for iAlign = 1:nAlignUsed
                color = 'k';
                timesThis = squeeze(timesByAlign(1, :, iAlign));
                if p.Results.drawSpikeWaveforms && ~p.Results.quick
                    wavesThis = squeeze(wavesByAlign(1, :, iAlign));
                    % draw waveforms in lieu of ticks
                        TrialDataUtilities.Plotting.drawTickRaster(timesThis, ...
                        'xOffset', tOffsetByAlign(iAlign), ...
                        'color', color, ...
                        'waveCell', wavesThis, 'waveformTimeRelative', wavesTvec, ...
                        'normalizeWaveforms', false, ... % alrerady normalized to [0 1]
                        'waveScaleHeight', p.Results.spikeWaveformScaleHeight, 'waveScaleTime', p.Results.spikeWaveformScaleTime);
                else
                    % draw vertical ticks
                    TrialDataUtilities.Plotting.drawTickRaster(timesThis, ...
                        'xOffset', tOffsetByAlign(iAlign), ...
                        'color', color, ...
                        'tickHeight', p.Results.tickHeight);
                end
                hold(axh, 'on');
            end


            if ~p.Results.quick
                au = AutoAxis(axh);
            else
%                 set(gca, 'YTick', []);
            end

            % setup time axis markers
            if p.Results.quick
                td.alignSummarySet{1}.setupTimeAutoAxis('which', 'x', ...
                    'style', 'quick');
            else
                % use marks or tickBridges via AlignSummary
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    td.alignSummarySet{idxAlign}.setupTimeAutoAxis('which', 'x', ...
                        'style', p.Results.timeAxisStyle, 'tOffsetZero', tOffsetByAlign(iAlign));
                end
            end

%             if iscell(unitNames)
%                 unitNameStr = TrialDataUtilities.String.strjoin(unitNames, ',');
%             else
%                 unitNameStr = unitName;
%             end
%             TrialDataUtilities.Plotting.setTitleIfBlank(axh, '%s Unit %s', td.datasetName, unitNameStr);

            axis(axh, 'tight');
            if ~p.Results.quick
%                 au = AutoAxis(axh);
%                 au.axisMarginLeft = p.Results.axisMarginLeft; % make room for left hand side labels
%                 axis(axh, 'off');
                au.update();
            else
                box(axh, 'off');
            end
            hold(axh, 'off');
       end
    end

    methods(Static) % Utility drawing methods
        function retInfo = plotConditionAlignedAnalogDataGroupMeans(D, varargin)
            % common utility function for drawing analog data averaged by
            % condition, used by the plotAnalogGroupMeans functions below.
            % Can handle multiple alignments. If a single alignment is
            % specified, defaults to the active alignment. Override by
            % specifying parameter 'alignIdx'. If nAlign alignments are
            % specified, uses all of them.
            %
            % time is either:
            %     vector, for nAlign == 1
            %     nAlign x 1 cell of ime vectors for each align
            %
            %  data, dataError is either
            %     nConditions x T x D matrix
            %     nAlign x 1 cell of nConditions x T x D matrices
            %
            %  quantilesData is either
            %      nConditions x T x nRandomized matrix
            %      nAlign x 1 cell of nConditions x T x nRandomized matrices

            p = inputParser();
            p.addParameter('time', [], @(x) isvector(x) || iscell(x)); % for D == 1,2,3 (for marking)
            p.addParameter('data', {}, @(x) isnumeric(x) || iscell(x)); % for D == 1,2,3
            p.addParameter('dataError', {}, @(x) isnumeric(x) || iscell(x));

            p.addParameter('conditionDescriptor', [], @(x) isa(x, 'ConditionDescriptor'));
            p.addParameter('alignSummarySet', [], @iscell);
            p.addParameter('alignInfoActiveIdx', 1, @isscalar);
            p.addParameter('alignTimeOffsets', [], @(x) isempty(x) || isvector(x));

            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('zOffset', 0, @isscalar);

            p.addParameter('axisInfoX', [], @(x) isempty(x) || isstringlike(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || isstringlike(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoZ', [], @(x) isempty(x) || isstringlike(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));

            p.addParameter('labelX', '', @isstringlike);
            p.addParameter('labelY', '', @isstringlike);
            p.addParameter('labelZ', '', @isstringlike);

            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @isstringlike);
            p.addParameter('axisStyleY', 'tickBridge', @isstringlike);

            p.addParameter('style', 'line', @isstringlike); % line, stairs
            p.addParameter('binAlignmentMode', BinAlignmentMode.Acausal, @(x) isa(x, 'BinAlignmentMode')); % used for stairs style to draw the staircase appropriately
            p.addParameter('binWidth', 0, @isscalar); % used for stairs to draw the staircase appropriately

            p.addParameter('lineWidth', get(0, 'DefaultLineLineWidth'), @isscalar);
            p.addParameter('conditionIdx', [], @isnumeric);
            p.addParameter('alignIdx', [], @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));

            p.addParameter('alpha', 1, @isscalar); % alpha for main traces
            p.addParameter('errorAlpha', 0.5, @isscalar); % alpha for surrounding error fills

            p.addParameter('quantileData', [], @(x) isnumeric(x) || iscell(x));

            p.addParameter('showMarkLabels', true, @islogical);
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('markShowOnAxis', true, @islogical);
            p.addParameter('markShowInLegend', true, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 4, @isscalar);
            p.addParameter('markOutline', true, @islogical);
            p.addParameter('markOutlineColor', 'w', @(x) true);

            p.addParameter('showIntervalLabels', true, @islogical);
            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalShowOnAxis', true, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);

            p.addParameter('showRangesOnData', false, @islogical); % show ranges for marks on traces
            p.addParameter('showRangesOnAxis', true, @islogical); % show ranges for marks below axis

            p.addParameter('timeAxisStyle', 'tickBridge', @isstringlike);
            p.addParameter('timeScaleBarWidth', NaN, @isscalar);

            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('useTranslucentMark3d', true, @islogical);

            p.addParameter('quick', false, @islogical);
            p.addParameter('clickable', false, @islogical);
            p.addParameter('Clipping', false, @islogical);

            p.addParameter('axh', [], @(x) true); % pass thru to getRequestedPlotAxis
            p.addParameter('deferUpdate', false, @islogical);

            p.addParameter('retInfo', struct(), @isstruct);

            p.KeepUnmatched = false;
            p.parse(varargin{:});

            retInfo = p.Results.retInfo;

            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            zOffset = p.Results.zOffset;

            clipping = p.Results.Clipping;
            plotOptions = [{'Clipping', clipping} p.Results.plotOptions{:}];

            % plot the mean and sem for an analog channel vs. time within
            % each condition
            import TrialDataUtilities.Plotting.errorshade;

            axh = TrialDataConditionAlign.getRequestedPlotAxis('axh', p.Results.axh);

            if ~ismember('scaleBars', p.UsingDefaults) && p.Results.scaleBars
                axisStyleX = 'scaleBar';
                axisStyleY = 'scaleBar';
            else
                axisStyleX = p.Results.axisStyleX;
                axisStyleY = p.Results.axisStyleY;
            end

            cd = p.Results.conditionDescriptor;

            conditionIdx = p.Results.conditionIdx;
            if isempty(conditionIdx)
                conditionIdx = 1:cd.nConditions;
            elseif islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            app = cd.appearances(conditionIdx);

            data = p.Results.data;
            dataError = p.Results.dataError;
            quantileData = p.Results.quantileData;
            time = p.Results.time;
            if isempty(data) || isempty(time)
                error('Must provide data and time');
            end

            plotErrorY = ~isempty(dataError);

            if iscell(time)
                assert(isvector(time) && isvector(data) && (isempty(dataError) || isvector(dataError)), ...
                    'Arguments must be nAlign x 1 cell or C x T x D matrices');
                nAlignUsed = numel(time);
            else
                assert(~iscell(time) && ~iscell(data) && (isempty(dataError) || ~iscell(dataError)) && (isempty(quantileData) || ~iscell(quantileData)), ...
                    'Arguments must be nAlign x 1 cell or C x T x D matrices');
                nAlignUsed = 1;

                time = {time};
                data = {data};
                if ~isempty(dataError)
                    dataError = {dataError};
                end
                if ~isempty(quantileData)
                    quantileData = {quantileData};
                end
            end

            % check size of time/data cell contents and mask conditions if
            % needed
            for iA = 1:nAlignUsed
                if size(time{iA}, 1) == nConditionsUsed
                    % okay as is
                elseif size(time{iA}, 1) == cd.nConditions
                    % needs to be masekd
                    time{iA} = time{iA}(conditionIdx, :, :);
                end
                if size(data{iA}, 1) == nConditionsUsed
                    % okay as is
                elseif size(data{iA}, 1) == cd.nConditions
                    % needs to be masekd
                    data{iA} = data{iA}(conditionIdx, :, :);
                end
                if ~isempty(dataError)
                    if size(dataError{iA}, 1) == nConditionsUsed
                        % okay as is
                    elseif size(dataError{iA}, 1) == cd.nConditions
                        % needs to be masekd
                        dataError{iA} = dataError{iA}(conditionIdx, :, :);
                    end
                end
                if ~isempty(quantileData)
                    if size(quantileData{iA}, 1) == nConditionsUsed
                        % okay as is
                    elseif size(quantileData{iA}, 1) == cd.nConditions
                        % needs to be masked
                        quantileData{iA} = quantileData{iA}(conditionIdx, :, :);
                    end
                end
            end

            % figure out which alignments are used
            alignSummarySet = p.Results.alignSummarySet;
            if isempty(p.Results.alignIdx)
                if nAlignUsed == 1
                    alignIdx = p.Results.alignInfoActiveIdx;
                else
                    alignIdx = 1:numel(alignSummarySet);
                end
            else
                alignIdx = 1:numel(alignSummarySet);
                alignIdx = alignIdx(p.Results.alignIdx);
            end
            assert(numel(data) == nAlignUsed, 'Number of alignments must match numel(data)');
            if ~isempty(dataError)
                assert(numel(dataError) == nAlignUsed, 'Number of alignments must match numel(dataError)');
            end
            if ~isempty(quantileData)
                assert(numel(quantileData) == nAlignUsed, 'Number of alignments must match numel(quantileData)');
            end
            assert(numel(alignIdx) == nAlignUsed, 'numel(time) must match numel(alignIdx)');

            % compute time offsets between successive alignments
            if D == 1
                timeOffsetByAlign = p.Results.alignTimeOffsets;
            else
                timeOffsetByAlign = zeros(nAlignUsed, 1);
            end

            % store handles as we go
            hData = gobjects(nConditionsUsed, nAlignUsed);

            % plot data traces
            hold(axh, 'on');
            for iAlign = 1:nAlignUsed
                for iCond = 1:cd.nConditions
                    tvec = time{iAlign};
                    % these will be T X D
                    dmat = squeeze(data{iAlign}(iCond, :, :));
                    if ~isempty(dataError)
                        errmat = squeeze(dataError{iAlign}(iCond, :, :));
                    end

                    if all(isnan(dmat(:))), continue; end

                    if D == 1
                        tOffset = timeOffsetByAlign(iAlign);
                        plotArgs = app(iCond).getPlotArgs();

                        if ~isempty(quantileData)
                            qmat = squeeze(quantileData{iAlign}(iCond, :, :)); % should be T by nQuantiles
                            % plot quantiles first
                            hQuant = plot(axh, tvec + tOffset + xOffset, qmat + yOffset, '-', ...
                                'LineWidth', p.Results.lineWidth / 2, 'Parent', axh, ...
                                plotArgs{:}, plotOptions{:});
                            TrialDataUtilities.Plotting.hideInLegend(hQuant);
                        end

                        if strcmp(p.Results.style, 'line')
                            if numel(tvec) > 1
                                if plotErrorY
                                    hShade = TrialDataUtilities.Plotting.errorshade(tvec + tOffset + xOffset, dmat + yOffset, ...
                                        errmat, app(iCond).Color, 'axh', axh, 'Clipping', clipping, ...
                                        'alpha', p.Results.errorAlpha, 'z', zOffset, 'showLine', false); % we'll plot the mean line ourselves
                                    TrialDataUtilities.Plotting.hideInLegend(hShade);
                                end
                                hData(iCond, iAlign) = plot(axh, tvec + tOffset + xOffset, dmat + yOffset, '-', ...
                                    'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                    plotArgs{:}, plotOptions{:});
                            else
                                if plotErrorY
                                    [hData(iCond, iAlign), hShade] = TrialDataUtilities.Plotting.errorline(tvec + tOffset + xOffset, dmat + yOffset, errmat, ...
                                        'Color', app(iCond).Color, 'axh', axh, 'LineWidth', p.Results.lineWidth, 'MarkerSize', 8);
                                    TrialDataUtilities.Plotting.hideInLegend(hShade);
                                else
                                    hData(iCond, iAlign) = plot(tvec + tOffset + xOffset, dmat + yOffset, 'o', 'Parent', axh, ...
                                        'MarkerFaceColor', app(iCond).Color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
                                end
                                TrialDataUtilities.Plotting.hideInLegend(hShade);
                            end
                            if p.Results.alpha < 1
                                TrialDataUtilities.Plotting.setLineOpacity(hData(iCond, iAlign), p.Results.alpha);
                            end

                        elseif strcmp(p.Results.style, 'stairs')
                            % offset the plot so as to resemble the binning
                            % mode used
                            xBinOffset = p.Results.binAlignmentMode.getBinStartOffsetForBinWidth(p.Results.binWidth);

                            if plotErrorY
                                 [hData(iCond, iAlign), hShade] = TrialDataUtilities.Plotting.stairsError(...
                                     tvec + tOffset + xOffset + xBinOffset, dmat + yOffset, errmat, ...
                                     'axh', axh, 'errorAlpha', p.Results.errorAlpha, 'color', app(iCond).Color, 'alpha', p.Results.alpha, ...
                                     'errorStyle', 'fill', 'errorColor', app(iCond).Color, 'lineWidth', p.Results.lineWidth);
                                TrialDataUtilities.Plotting.hideInLegend(hShade);
                            else
                                hData(iCond, iAlign) = TrialDataUtilities.Plotting.stairs(...
                                 tvec + tOffset + xOffset + xBinOffset, dmat + yOffset, ...
                                 'axh', axh, 'color', app(iCond).Color, 'alpha', p.Results.alpha, ...
                                 'lineWidth', p.Results.lineWidth);
                            end
                        end

                    elseif D == 2
                        if p.Results.alpha < 1
                           hData(iCond, iAlign) = TrialDataUtilities.Plotting.patchline(dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, ...
                               'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                               'LineWidth', p.Results.lineWidth, 'Parent', axh, plotOptions{:});
                        else
                            plotArgs = app(iCond).getPlotArgs();
                            hData(iCond, iAlign) = plot(axh, dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, '-', ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                plotArgs{:}, plotOptions{:});
                        end

                    elseif D == 3
                        if p.Results.alpha < 1
                            hData(iCond, iAlign) = TrialDataUtilities.Plotting.patchline(dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, dmat(:, 3) + zOffset, ...
                               'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                               'LineWidth', p.Results.lineWidth,  'Parent', axh, plotOptions{:});
                        else
                            plotArgs = app(iCond).getPlotArgs();
                            hData(iCond, iAlign) = plot3(axh, dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, dmat(:, 3) + zOffset, '-', ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                plotArgs{:}, plotOptions{:});
                        end
                    end

                    % update name for inclusion in legend
                    if iAlign == 1
                        TrialDataUtilities.Plotting.showFirstInLegend(hData(iCond, iAlign), cd.namesShort{iCond});
                    else
                        TrialDataUtilities.Plotting.hideInLegend(hData(iCond, iAlign));
                    end
                end

                box(axh, 'off');
                axis(axh, 'tight');
            end

            % draw marks and intervals on data
            if ~p.Results.quick && ~isempty(alignSummarySet)
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    alignSummarySet{idxAlign}.drawOnDataByCondition(time{iAlign}, ...
                        permute(data{iAlign}, [2 3 1]), ...  % data needs to be T x D x C x N, currently C x T x D
                        'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
                        'xOffset', xOffset, 'yOffset', yOffset, 'zOffset', zOffset, ...
                        'tOffsetZero', timeOffsetByAlign(iAlign), ...
                        'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize, ...
                        'markOutline', p.Results.markOutline, 'markOutlineColor', p.Results.markOutlineColor, ...
                        'intervalAlpha', p.Results.intervalAlpha, ...
                        'showRanges', p.Results.showRangesOnData, ...
                        'showInLegend', p.Results.markShowInLegend', ...
                        'tMin', min(time{iAlign}), 'tMax', max(time{iAlign}), ...
                        'style', p.Results.style);
                end
            end

            % setup time axes for each alignment
            if ~p.Results.quick && ~isempty(alignSummarySet)
                if D == 1
                    for iAlign = 1:nAlignUsed
                        idxAlign = alignIdx(iAlign);

                        % setup x axis
                        alignSummarySet{idxAlign}.setupTimeAutoAxis('axh', axh, 'tOffsetZero', timeOffsetByAlign(iAlign) + xOffset, ...
                            'tMin', min(time{iAlign}), 'tMax', max(time{iAlign}), 'timeScaleBarWidth', p.Results.timeScaleBarWidth, ...
                            'style', p.Results.timeAxisStyle,  'showIntervals', p.Results.intervalShowOnAxis, ...
                            'showMarks', p.Results.markShowOnAxis, 'showRanges', p.Results.showRangesOnAxis, ...
                            'showMarkLabels', p.Results.showMarkLabels, 'showIntervalLabels', p.Results.showIntervalLabels);
                    end
                end
            end

            % setup non-time axes
            if ~p.Results.quick
                if D == 1
                    if isstringlike(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                        % x is data, y is time
                        if ~isempty(p.Results.axisInfoX)
                            TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                                'which', 'x', 'style', axisStyleX, 'label', p.Results.labelX);
                        end
                    else
                        % y is data, x is time
                        if ~isempty(p.Results.axisInfoY)
                            TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                                'which', 'y', 'style', axisStyleY, 'label', p.Results.labelY);
                        end
                    end

                elseif D == 2
                    if ~isempty(p.Results.axisInfoX)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                            'which', 'x', 'style', axisStyleX, 'label', p.Results.labelX);
                    end
                    if ~isempty(p.Results.axisInfoY)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                            'which', 'y', 'style', axisStyleY, 'label', p.Results.labelY);
                    end

                elseif D == 3
                    if ~isempty(p.Results.axisInfoX)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                            'which', 'x', 'useAutoAxis', false, 'label', p.Results.labelX);
                    end
                    if ~isempty(p.Results.axisInfoY)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                            'which', 'y', 'useAutoAxis', false, 'label', p.Results.labelY);
                    end
                    if ~isempty(p.Results.axisInfoZ)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoZ, ...
                            'which', 'z', 'useAutoAxis', false, 'label', p.Results.labelZ);
                    end

                    if p.Results.useThreeVector
                        ThreeVector(axh);
                        axis(axh, 'off');
                        axis(axh, 'vis3d');
                    end
                end
            end

            if p.Results.clickable
                for iAlign = 1:nAlignUsed
                    TrialDataUtilities.Plotting.makeClickableShowDescription(hData(:, iAlign), cd.names(:));
                end
            end

            box(axh, 'off');
            axis(axh, 'tight');

            if ~p.Results.quick && ~p.Results.deferUpdate
                au = AutoAxis(axh);
                au.update();
            end

            hold(axh, 'off');
            set(axh, 'SortMethod', 'childorder');
        end
    end

    % Axis randomization settings
    methods
        % wipe out any randomization
        function td = withoutRandomization(td, varargin)
            td.warnIfNoArgOut(nargout);
            % restore original condition info, but don't flush randomized
            td.conditionInfo = td.conditionInfo.noRandomization(varargin{:});
            td = td.postUpdateConditionInfo(false);
        end

        function assertHasRandomizationSpecified(td)
            assert(td.hasRandomizationSpecified, 'TrialData has no randomization active. Try .setRandomizedResampleTrialsWithinConditions()');
        end

        function td = withRandomized(td, idxRandom)
            % make this TD instance act as a specific draw of the
            % randomization, so as to visualize specific reshuffles,
            % resamplings, etc.
            % this affects data coming from the trial data instance itself,
            % all methods below only affect the random data returned by
            % methods ending in groupedRandomized
            td.warnIfNoArgOut(nargout);
            td.assertHasRandomizationSpecified();
            if nargin < 2
                idxRandom = 1;
            end
            % use conditionInfoRandomized's seed + (idxRandom-1)
            % so as to match the results of random sampling
            td.conditionInfo = td.conditionInfoRandomized.setRandomSeed(td.conditionInfoRandomized.randomSeed + idxRandom-1);
            td = td.postUpdateConditionInfo(false);
        end

        function td = setRandomSeed(td, seed)
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.setRandomSeed(seed);
        end

        function td = newRandomSeed(td, varargin)
            % no args
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.newRandomSeed(varargin{:});
        end

        function td = setRandomizedResampleTrialsWithinConditions(td, varargin)
            % no args
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.resampleTrialsWithinConditions(varargin{:});
        end

        function td = setRandomizedAxisNoRandomization(td, varargin)
            % args: idxOrAttr
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.axisNoRandomization(varargin{:});
        end

        function td = setRandomizedAxisShuffle(td, varargin)
            % args (idxOrAttr, replace)
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.axisShuffle(varargin{:});
        end

        function td = setRandomizedAxisResampleFromSpecifiedValueListIndices(td, varargin)
            % args: (axisIdxOrAttr, resampleFromIndices, replace)
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.axisResampleFromSpecifiedValueListIndices(varargin{:});
        end

        function td = axisResampleFromSpecifiedValues(td, varargin)
            % args: (axisIdxOrAttr, resampleFromValueStructMatch, replace)
            td.warnIfNoArgOut(nargout);
            td.conditionInfoRandomized = td.conditionInfoRandomized.axisResampleFromSpecifiedValues(varargin{:});
        end

        function tf = get.hasRandomizationActive(td)
            tf = td.conditionInfo.hasRandomization;
        end

        function tf = get.hasRandomizationSpecified(td)
            tf = td.conditionInfoRandomized.hasRandomization;
        end

        function v = get.randomizationDescription(td)
            v = td.conditionInfoRandomized.randomizationDescription;
        end

        function td = setNumRandomized(td, v)
            td.warnIfNoArgOut(nargout);
            if v ~= td.nRandomized
                td.randomizedListsByCondition = [];
            end
            td.nRandomized = v;
        end

        function lists = buildRandomizedListsByCondition(td)
            if ~td.hasRandomizationSpecified
                lists = [];
            else
                lists = td.conditionInfoRandomized.generateMultipleRandomizedListByCondition(td.nRandomized, 'showProgress', true);
                lists = reshape(lists, [td.conditionsSizeNoExpand td.nRandomized]); % make (conditionsSize) x nRandomized
            end
        end

        % given a cellvec or nmeric vector, group its elements according to
        % a specific randomized list
        function varargout = groupElementsRandomizedSingle(td, idxRandom, varargin)
            assert(isscalar(idxRandom));
            td.assertHasRandomizationSpecified();
            lists = TensorUtils.selectAlongDimension(td.randomizedListsByCondition, td.conditionInfoRandomized.nAxes+1, idxRandom, true); % squeeze result
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) TensorUtils.selectAlongDimension(data, 1, idx, false), ...
                    lists, 'UniformOutput', false);
            end
        end

        % given a cellvec or nmeric vector, group its elements according to
        % a specific randomized list
        function varargout = groupElementsRandomized(td, varargin)
            td.assertHasRandomizationSpecified();
            lists = td.randomizedListsByCondition;
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) TensorUtils.selectAlongDimension(data, 1, idx, false), ...
                    lists, 'UniformOutput', false);
            end
        end
    end
end
