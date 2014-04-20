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
    end
    
    % On Demand Cache handle
    properties(Access=protected)
        odc % TrialDataConditionAlignOnDemandCache instance
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
    end

    % Properties which read through to ConditionInfo
    properties(Dependent, SetAccess=protected)
        nConditions
        listByCondition
        conditionIdx
        conditionSubs
        conditions
        conditionsAsStrings
        conditionNames
        conditionAppearances
        
        axisValueLists
        axisValueListsAsStrings
        
        conditionAppearanceFn
    end
    
    % Simple dependent properties
    properties(Dependent, SetAccess=protected)
        nAlign

        alignInfoActive % alignInfo currently active

        alignSummaryActive % alignSummary currently active
    end

    % Properties which read through to AlignInfo
    properties(Dependent, SetAccess=protected)
        minTimeDelta
    end
    
    % Initializing and building
    methods
        function td = TrialDataConditionAlign(varargin)
            td = td@TrialData(varargin{:});
            td.odc = TrialDataConditionAlignOnDemandCache();
            td = td.initializeConditionInfo();
            td = td.initializeAlignInfo();
            td = td.updateValid();
        end
    end
    
    % Simple dependent getters
    methods
        function v = get.nAlign(td)
            v = numel(td.alignInfoSet);
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
    end
    
    % build* methods for properties stored in odc
    methods
        function buildEventData(td)
            % build a cached version of event times into a matrix for easy alignment
            % sets .eventData and .eventCounts
            
            evStruct = td.getRawEventFlatStruct();
            evList = fieldnames(evStruct);
            nEvents = numel(evList);
            for iE = 1:nEvents
                ev = evList{iE};
                times = evStruct.(ev);
                if iscell(times)
                    % event may happen zero, one, or multiple times
                    % convert to nTrials x nOccur matrix
                    counts = cellfun(@numel, times); 
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

        function buildAlignSummarySet(td)
            alignSummarySet = cell(td.nAlign, 1);
            for i = 1:td.nAlign
                alignSummarySet{i} = AlignSummary.buildFromConditionAlignInfo(td.conditionInfo, td.alignInfoSet{i});
                alignSummarySet{i}.timeUnitName = td.timeUnitName;
            end

            c = td.odc;
            c.alignSummarySet = alignSummarySet;
        end
    end
    
    % General utilites
    methods
        % print a short description
        function disp(td)
            td.printDescriptionShort();
            
            td.conditionInfo.printDescription();
            for iA = 1:td.nAlign
                td.alignInfoSet{iA}.printDescription();
            end
            
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end

        % synchronize valid between AlignInfo and ConditionINfo
        function td = updateValid(td)
            td.warnIfNoArgOut(nargout);

            if isempty(td.manualValid)
                td.manualValid = truevec(td.nTrials);
            end

            valid = td.buildValid();

            td.conditionInfo = td.conditionInfo.setInvalid(~valid);
            for iA = 1:td.nAlign
                td.alignInfoSet{iA} = td.alignInfoSet{iA}.setInvalid(~valid);
            end
        end
        
        function valid = buildValid(td)
            % valid is the intersection of manualValid, conditionInfo valid,
            % and all AlignInfo valid
            cvalid = td.conditionInfo.computedValid;
            avalid = truevec(td.nTrials);
            for iA = 1:td.nAlign
                avalid = avalid & td.alignInfoSet{iA}.computedValid;
            end

            if ~isempty(td.manualValid)
                valid = td.manualValid & cvalid & avalid;
            else
                valid = cvalid & avalid;
            end
        end
        
        function cause = buildInvalidCause(td)
            cause = cell(td.nTrials, 1);
            explained = false(td.nTrials, 1);
            
            cause(~td.manualValid & ~explained) = {'marked invalid manually'};
            
            % invalid by condition info
            cmask = ~td.conditionInfo.computedValid & ~explained;
            cause(cmask) = cellfun(@(s) ['ConditionInfo: ' s], ...
                td.conditionInfo.invalidCause(cmask), 'UniformOutput', false);
            
            % invalid by each align info
            for iA = 1:td.nAlign
                amask = ~td.alignInfoSet{iA}.computedValid & ~explained;
                cause(amask) = cellfun(@(s) ['AlignInfo ', num2str(iA), ': ', s], ...
                    td.alignInfoSet{iA}.invalidCause(amask), 'UniformOutput', false);
            end
            
            cause(td.valid) = {''};
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
            names = wrapCell(names);
            
            % check whether any of the alignInfo events and error if so
            for iA = 1:td.nAlign
                alignEvents = td.alignInfoSet{iA}.getEventList();
                mask = ismember(names, alignEvents);
                if any(mask)
                    error('TrialData alignment depends on event %s', ...
                        strjoin(alignEvents(mask)));
                end
            end
            
            % check whether any of the events are in condition info
            conditionParams = td.conditionInfo.attributeNames;
            mask = ismember(names, conditionParams);
            if any(mask)
                error('TrialData conditioning depends on params %s', ...
                    strjoin(conditionParams(mask)));
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
    end

    % ConditionInfo control
    methods
        % parameters that are either scalar or strings
        function names = listConditionInfoCompatibleParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor') && ...
                (cd.isScalarByField{1} || cd.isStringByField{1}), ...
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
        
        % given data with dimension 1 with size nTrials, group by condition
        % and map out{i} = fn(group{i})
        function out = mapByGroup(td, fn, varargin)
            dataByGroup = td.groupElements(varargin{:});
            out = cellfun(fn, dataByGroup, 'UniformOutput', false);
        end

        
        function td = setConditionDescriptor(td, cd)
            td.warnIfNoArgOut(nargout);
            
            % manually accept a condition descriptor instance
            assert(isequal(class(cd), 'ConditionDescriptor'), 'Must be a ConditionDescriptor instance');
            
            % grab the param data to feed to the condition descriptor
            paramData = td.getRawChannelDataAsStruct(cd.attributeRequestAs);
            
            % build condition info from condition descriptor
            td.conditionInfo = ConditionInfo.fromConditionDescriptor(cd, paramData);
            td = td.postUpdateConditionInfo();
        end
        
        function td = setConditionAppearanceFn(td, fn)
            % Update the appearanceFn callback of conditionDescriptor
            % without invalidating any of the other cached info
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.appearanceFn = fn;
        end
        
        function td = colorByAttributes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.colorByAttributes(varargin{:});
        end
        
        function td = colorByAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.colorByAxes(varargin{:});
        end


        function td = selectTrials(td, mask)
            % Apply a logical mask or index selection to the list of trials
            % within, appropriately notifying the condition descriptor and
            % align descriptor within
            td.warnIfNoArgOut(nargout);
            td = selectTrials@TrialData(td, mask);
            td.conditionInfo = td.conditionInfo.selectTrials(mask);
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
                
            td.alignSummarySet = [];
        end
        
        function td = addAttribute(td, names)
            % add attributes in names that aren't already in ConditionInfo
            new = setdiff(names, td.conditionInfo.attributeNames);
            
            for i = 1:numel(new)
                td.conditionInfo = td.conditionInfo.addAttribute(new{i}, ...
                    'values', td.getParamRaw(new{i}));
            end
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
        
        function td = ungroup(td) 
            % this only undoes the grouping axes, NOT the value list
            % filtering. use reset condition info for that
            td.warnIfNoArgOut(nargout);
            td = td.groupBy();
        end
        
        function td = reshapeAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.reshapeAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = flattenAxes(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.flattenAxes(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
       function td = selectConditions(td, varargin)
            % select specific conditions by linear index or mask
            % and return a single-axis condition descriptor with just those
            % conditions selected
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.selectConditions(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = selectConditionsAlongAxis(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.selectConditionsAlongAxis(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = matchSelectConditionsAlongAxis(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.matchSelectConditionsAlongAxis(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = resetConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = [];
            td = td.initializeConditionInfo();
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
            td.conditionInfo = td.conditionInfo.setAttributeValueListAuto(attr);
            td = td.postUpdateConditionInfo();
        end
        
        function td = binAttribute(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttribute(varargin{:});
            td = td.postUpdateConditionInfo();
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
        
        function td = setAttributeValueList(td, attrName, valueList)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(attrName);
            td.conditionInfo = td.conditionInfo.setAttributeValueList(attrName, valueList);
            td = td.postUpdateConditionInfo();
        end
        
        function valueList = getAttributeValueList(td, attrName)
            td = td.addAttrbute(attrName);
            valueList = td.conditionInfo.getAttributeValueList(attrName);
        end

        function td = setAxisValueList(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAxisValueList(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function valueList = getAxisValueList(td, varargin)
            valueList = td.conditionInfo.getAxisValueList(varargin{:});
        end
        
        function td = postUpdateConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
            td.alignSummarySet = [];
        end
        
        % filter trials that are valid based on ConditionInfo
        function td = filterValidTrialsConditionInfo(td, varargin)
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
            td = td.filterValidTrialsConditionInfo();
        end

        % filter trials based on matching attribute values, 
        % pass thru to ConditionInfo
        function td = filteredByAttributeStruct(td, varargin)
            % call on ConditionInfo, and filter based on its valid trials
            td.warnIfNoArgOut(nargout);
            % no need to copy, filteredBy already copies the conditionInfo
            td.conditionInfo = td.conditionInfo.filteredByAttributeStruct(varargin{:});
            td = td.filterValidTrialsConditionInfo();
        end

        function tdCell = getTrialDataByCondition(td)
            tdCell = arrayfun(@(condition) td.filteredByAttributeStruct(condition), ...
                td.conditions, 'UniformOutput', false);
        end

        function v = get.conditions(td)
            v = td.conditionInfo.conditions;
        end
        
        function v = get.conditionsAsStrings(td)
            v = td.conditionInfo.conditionsAsStrings;
        end

        function n = get.nConditions(td)
            n = td.conditionInfo.nConditions;
        end
        
        function v = get.axisValueLists(td)
            v = td.conditionInfo.axisValueLists;
        end
        
        function v = get.axisValueListsAsStrings(td)
            v = td.conditionInfo.axisValueListsAsStrings;
        end

        function v = get.conditionNames(td)
            v = td.conditionInfo.names;
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
        
        function v = get.listByCondition(td)
            v = td.conditionInfo.listByCondition;
        end

        function v = get.conditionIdx(td)
            v = td.conditionInfo.conditionIdx;
        end

        function v = get.conditionSubs(td)
            v = td.conditionInfo.conditionSubs;
        end
        
        function minTimeDelta = get.minTimeDelta(td)
            minTimeDelta = td.alignDescriptor.minTimeDelta;
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
            td = td.postUpdateAlignInfo();
        end
        
        function td = postUpdateAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
            % cause align summary to be recomputed
            td.alignSummarySet = [];
        end

        function td = align(td, varargin)
            td.warnIfNoArgOut(nargout);

            adSet = cell(numel(varargin), 1);
            for i = 1:numel(varargin)
                ad = varargin{i};

                if iscell(ad)
                    error('Please provide alignDescriptors as successive arguments');
                end
                if ischar(ad)
                    adSet{i} = AlignInfo(ad);
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
        
        function td = addAlign(td, varargin)
            % whereas align replaces the existing align descriptor set,
            % this appends an align descriptor
            td.warnIfNoArgOut(nargout);

            adSet = cat(1, td.alignInfoSet, makecol(varargin));
            td = td.align(adSet{:});
            
            % and for convenience switch to make this new alignment active
            td = td.useAlign(td.nAlign);
        end
        
        function td = unalign(td)
            td.warnIfNoArgOut(nargout);
            td = td.align('TrialStart:TrialEnd');
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
        
        % the following methods pass-thru to alignInfo:
        
        function td = pad(td, window)
            % add a padding window to the AlignInfo
            % may change which trials are valid
            % usage: pad([pre post]) or pad(pre, post)
            % pre > 0 means add padding before the start (typical case)
            td.warnIfNoArgOut(nargout);
            
            updated = false;
            for iA = 1:td.nAlign
                if td.alignInfoSet{iA}.padPre ~= window(1) || td.alignInfoSet{iA}.padPost ~= window(2)
                    td.alignInfoSet{iA} = td.alignInfoSet{iA}.pad(window);
                    updated = true;
                end
            end

            if updated
                % to synchronize any changes in alignment validity
                td = td.postUpdateAlignInfo();
            end
        end
        
        function td = padForSpikeFilter(td, sf)
            % Pad trial data alignment for spike filter
            td.warnIfNoArgOut(nargout);
            td = td.pad([sf.preWindow sf.postWindow]);
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
        
        function td = mark(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.mark(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = interval(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfoActive = td.alignInfoActive.interval(varargin{:});
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
            assert(isscalar(idx) && isnumeric(idx) && idx >= 1 && idx <= td.nAlign);
            td.alignInfoActiveIdx = idx;
        end

        function durations = getValidDurations(td)
            durations = td.td.alignInfoActive.getValidDurationByTrial();
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
            [tMinByTrial, tMaxByTrial] = td.alignInfoActive.getStartStopRelativeToZeroByTrial();

            tMinByTrial(~td.valid) = NaN;
            tMaxByTrial(~td.valid) = NaN;
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

        function offsets = getTimeOffsetsFromZeroEachTrial(td)
            offsets = td.alignInfoActive.getZeroByTrial();
            offsets(~td.valid) = NaN;
        end

        function offsets = getTimeOffsetsFromZeroEachTrialEachAlign(td)
            offsets = nan(td.nTrials, td.nAlign);
            for iA = 1:pset.nAlign
                offsets(:, iA) = td.alignInfoSet{iA}.getZeroByTrial();
            end
            offsets(~td.valid, :) = NaN;
        end
    end

    % Analog channel access
    methods
        % return aligned analog channel
        function [data, time] = getAnalog(td, name)
            [data, time] = getAnalog@TrialData(td, name);
            [data, time] = td.alignInfoActive.getAlignedTimeseries(data, time, false);
        end
        
        function tvec = getAnalogCommonTimeVector(nameCell, varargin)
            p = inputParser();
            p.addParamValue('timeDelta', [], @isscalar); 
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
            timeCell = cell(td.nTrials, C);
            for c = 1:C
                [~, timeCell(:, c)] = td.getAnalog(nameCell{c});
            end
            
            tvec = TrialDataUtilies.Data.inferCommonTimeVectorForTimeseriesData(timeCell, ...
                'timeDelta', timeDelta);
        end
        
        function [mat, tvec] = getAnalogAsMatrix(td, name, varargin)
            % return aligned analog channel, resampled and interpolated to
            % a uniformly spaced time vector around t=0 such that the
            % result can be embedded in a nTrials x nTime matrix. time will
            % be chosen to encapsulate the min / max timestamps across all
            % trials. Missing samples will be returned as NaN
            %
            % if name is a cellstr of multiple channels, mat will be nTrials x nTime x nChannels
            % tensor of values
            
            p = inputParser;
            p.addParamValue('timeDelta', [], @isscalar);
            p.addParamValue('tvec', [], @isvector);
            p.parse(varargin{:});

            if ischar(name)
                name = {name};
            end
            
            if isempty(p.Results.tvec)
                % infer timeDelta to use
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = td.alignInfoActive.minTimeDelta;
                    if isempty(timeDelta)
                        timeDelta = td.getAnalogTimeDelta(name);
                        warning('timeDelta auto-computed from analog timestamps. Specify manually or call .round for consistent results');
                    end
                end
            end
            
            % build nTrials x nChannels cell of data/time vectors
            C = numel(name);
            [dataCell, timeCell] = deal(cell(td.nTrials, C));
            for c = 1:C
                [dataCell(:, c), timeCell(:, c)] = td.getAnalog(name{c});
            end

            % interpolate to common time vector
            [mat, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell, ...
                'timeDelta', timeDelta, 'timeReference', 0, 'tvec', p.Results.tvec);
        end

        function [data, tvec] = getMultiAnalogAsMatrix(td, name, varargin)
            p = inputParser;
            p.addParamValue('timeDelta', [], @isscalar);
            p.addParamValue('tvec', [], @isvector);
            p.parse(varargin{:});

            if ischar(name)
                name = {name};
            end
            
            % build nTrials x nChannels cell of data/time vectors
            C = numel(name);
            [dataCell, timeCell] = deal(cell(td.nTrials, C));
            for c = 1:C
                [dataCell(:, c), timeCell(:, c)] = td.getAnalog(name{c});
            end
            
            if isempty(p.Results.tvec)
                % infer timeDelta to use
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = td.alignInfoActive.minTimeDelta;
                    if isempty(timeDelta)
                        timeDelta = td.getAnalogTimeDelta(name);
                        warning('timeDelta auto-computed from analog timestamps. Specify manually or call .round for consistent results');
                    end
                end
                
                % find common time vector
                tvec = TrialDataUtilities.Data.inferCommonTimeVectorForTimeseriesData(timeCell, ...
                    'timeDelta', timeDelta);
            else
                tvec = p.Results.tvec;
            end
            
            % interpolate to common time vector
            data = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell, 'tvec', tvec);
        end
        
        function [dCell, tCell] = getAnalogGrouped(td, name)
            [dataCell, timeCell] = td.getAnalog(name);
            [dCell, tCell] = td.groupElements(dataCell, timeCell);
        end
        
        function [dCell, tvec] = getAnalogAsMatrixGrouped(td, name)
            [mat, tvec] = td.getAnalogAsMatrix(name);
            dCell = td.groupElements(mat);
        end
        
        function [dataCell, tvec] = getMultiAnalogAsMatrixGrouped(td, nameCell, varargin)
            [data, tvec] = td.getMultiAnalogAsMatrix(nameCell);
            dataCell = td.groupElements(data);
        end
        
        function [meanMat, semMat, tvec, nTrialsMat, stdMat] = getAnalogGroupMeans(td, name, varargin)
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParamValue('minTrials', 1, @isscalar); % minimum trial count to average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            [dCell, tvec] = td.getAnalogAsMatrixGrouped(name);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec)));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :), semMat(iC, :), nTrialsMat(iC, :), stdMat(iC, :)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials);
                end
            end 
        end
        
        function [dataCell, timeCell] = getAnalogSampleGrouped(td, name, varargin)
            [dataVec, timeVec] = td.getAnalogSample(name);
            [dataCell, timeCell] = td.groupElements(dataVec, timeVec);
        end
    end

    % Event channel access
    methods
        % return aligned event times
        function timesCell = getEvent(td, name)
            timesCell = getEvent@TrialData(td, name);
            timesCell = td.alignInfoActive.getAlignedTimesCell(timesCell, false);
        end

        function dCell = getEventGrouped(td, name)
            dCell = td.groupElements(td.getEvent(name));
        end

        function td = addEvent(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = addEvent@TrialData(td, varargin{:});
            
            % force .eventData and .eventCounts to be recomputed
            td.eventData = [];
            td.eventCounts = [];
            td = td.applyAlignInfoSet();
            
        end
    end

    % Param channel access
    methods 
        function dCell = getParamGrouped(td, name)
            dCell = td.groupElements(td.getParam(name));
        end
    end

    % Spike data
    methods
        % return aligned unit spike times
        function [timesCell] = getSpikeTimes(td, unitName)
            timesCell = getSpikeTimes@TrialData(td, unitName);
            timesCell = td.alignInfoActive.getAlignedTimesCell(timesCell, true);
        end
        
        function [rateCell, timeCell] = getSpikeRateFiltered(td, unitName, varargin)
            p = inputParser;
            p.addParamValue('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            spikeCell = td.getSpikeTimes(unitName);
            timeInfo = td.alignInfoActive.timeInfo;
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond);
        end
           
        function [rates, tvec] = getSpikeRateFilteredAsMatrix(td, unitName, varargin)
            p = inputParser;
            p.addParamValue('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.addParamValue('timeDelta', 1, @isscalar);
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            timeDelta = p.Results.timeDelta;
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            spikeCell = td.getSpikeTimes(unitName);
            timeInfo = td.alignInfoActive.timeInfo;
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rates, tvec] = sf.filterSpikeTrainsWindowByTrialAsMatrix(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond, ...
                'timeDelta', timeDelta);
        end

        function [rateCell, timeCell] = getSpikeRateFilteredGrouped(td, unitName, varargin)
            [rateCell, timeCell] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElements(rateCell);
            timeCell = td.groupElements(timeCell);
        end
        
        function [rateCell, tvec] = getSpikeRateFilteredAsMatrixGrouped(td, unitName, varargin)
            [rates, tvec] = td.getSpikeRateFilteredAsMatrix(unitName, varargin{:});
            rateCell = td.groupElements(rates);
        end
        
        function [psthMatrix, tvec, semMatrix] = getSpikeRateFilteredMeanByGroup(td, unitName, varargin)
            [rateCell, tvec] = getSpikeRateFilteredAsMatrixGrouped(td, unitName, varargin{:});
            % flatten the rateCell
            rateCell = rateCell(:);
            psthMatrix = cell2mat(cellfun(@(r) nanmean(r, 1), rateCell, 'UniformOutput', false));
            semMatrix =  cell2mat(cellfun(@(r) nansem(r, 1),  rateCell, 'UniformOutput', false));
        end
        
        function [psthMatrix, tvec, semMatrix] = getPSTH(td, varargin)
            [psthMatrix, tvec, semMatrix] = getSpikeRateFilteredMeanByGroup(td, unitName, varargin{:});
        end
        
        function timesCellofCells = getSpikeTimesGrouped(td, unitName)
            timesCell = td.getSpikeTimes(unitName);
            timesCellofCells = td.groupElements(timesCell);
        end
        
        function countsCell = getSpikeCountsGrouped(td, unitName)
            counts = td.getSpikeCounts(unitName);
            countsCell = td.groupElements(counts);
        end
        
        function rateCell = getSpikeMeanRateGrouped(td, unitName)
            rates = td.getSpikeMeanRate(unitName);
            rateCell = td.groupElements(rates);
        end
    end

    % Plotting
    methods
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
            %    time vectors, for common time across trials, or
            %    nTrials x 1 cell of time vectors, for different times
            %    across trials.
            %  
            % Parameters:
            %  alignIdx : vector indicating which alignments to include
            
            p = inputParser();
            p.addParamValue('alignIdx', 1:td.nAlign, @isvector);
            p.parse(varargin{:});
            
            % select alignment indices
            alignIdx = 1:td.nAlign;
            alignIdx = alignIdx(p.Results.alignIdx);
            nAlign = numel(alignIdx);
            
            offsets = nan(nAlign, 1);
            offsets(1) = 0;
            currentOffset = 0;
            
            % compute start/stop of each alignment
            [mins, maxs] = deal(nanvec(nAlign));
            for iAlign = 1:nAlign
                if iscell(tvecCell{iAlign})
                    mins(iAlign) = nanmin(cellfun(@nanminNanEmpty, tvecCell{iAlign}));
                    maxs(iAlign) = nanmax(cellfun(@nanmaxNanEmpty, tvecCell{iAlign}));
                else
                    mins(iAlign) = nanmin(tvecCell{iAlign});
                    maxs(iAlign) = nanmax(tvecCell{iAlign});
                end
            end
            
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
                currentOffset = currentOffset + nanmax(tvecCell{iAlign-1}) + ...
                    gaps(iAlign-1) - nanmin(tvecCell{iAlign});
                offsets(iAlign) = currentOffset;
            end
            
            lims = [mins(1), maxs(end) + offsets(end)];
        end
        
        function plotProvidedAnalogDataGroupedEachTrial(td, D, varargin)
            % common utility function for drawing analog data grouped by
            % condition, used by the plotAnalog functions below.
            % Can handle multiple alignments. If a single alignment is
            % specified, defaults to the active alignment. Override by
            % specifying parameter 'alignIdx'. If nAlign alignments are
            % specified, uses all of them.
            % 
            % time is either:
            %     common time vector
            %     cell with size (nConditions or 1) x nAlign cell containing either
            %       nTrials x 1 cell of time vectors or 
            %       time vector for all trials in that condition / alignment. 
 
            % for D==1: 1-D data vs. time as matrix
            %   data is nConditions x nAlign cell containing nTrials x 1
            %     cell of T x D matrices, or is nTrials x T x D data matrix
            % for D==2 or D==3: 2-D or 3-D data:
            %   data are nConditions x nAlign cell containing nTrials x 1
            %     cell of T x D matrices or containing nTrials x T x D data matrices
            
            p = inputParser();
            p.addParamValue('time', [], @(x) isvector(x) || iscell(x)); % for D == 1,2,3 (for marking)
            p.addParamValue('data', {}, @iscell); % for D == 1,2,3
            
            p.addParamValue('axisInfoX', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            p.addParamValue('axisInfoY', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            p.addParamValue('axisInfoZ', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            
            p.addParamValue('conditionIdx', 1:td.nConditions, @isnumeric);
            p.addParamValue('alignIdx', [], @isnumeric);
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.addParamValue('alpha', 1, @isscalar);
            p.addParamValue('markAlpha', 1, @isscalar);
            p.addParamValue('markSize', 5, @isscalar);
            p.addParamValue('timeAxisStyle', 'tickBridge', @ischar);
            p.KeepUnmatched;
            p.parse(varargin{:});
            
            axh = td.getRequestedPlotAxis(p.Unmatched);
            
            conditionIdx = p.Results.conditionIdx;
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            app = td.conditionAppearances(conditionIdx);
            
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
                    alignIdx = 1:nAlign;
                end
            else
                alignIdx = 1:td.nAlign;
                alignIdx = alignIdx(p.Results.alignIdx);
            end
            
            % compute time offsets between successive alignments
            if D == 1
                timeOffsetByAlign = td.getAlignPlottingTimeOffsets(time, 'alignIdx', alignIdx);
            end

            % store handles as we go
            hData = cell(nConditionsUsed, nAlignUsed);
            
            for iAlign = 1:nAlignUsed
                % plot one condition at a time
                for iCond = 1:nConditionsUsed
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
                            if p.Results.alpha < 1
                               hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(timeC' + tOffset, dataC', ...
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                            else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, tvec, dmat, '-', ...
                                    plotArgs{:}, p.Results.plotOptions{:});
                            end

                        elseif D == 2
                            if p.Results.alpha < 1
                               hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dataC(:, :, 1)', dataC(:, :, 2)', ...
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                            else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, dataC(:, :, 1)', dataC(:, :, 2)', '-', ...
                                    plotArgs{:}, p.Results.plotOptions{:});
                            end

                        elseif D == 3
                            if p.Results.alpha < 1
                                hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline3(dataC(:, :, 1)', dataC(:, :, 2)', dataC(:, :, 3)', ... 
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, 'axh', axh, p.Results.plotOptions{:});
                            else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot3(axh, dataC(:, :, 1)', dataC(:, :, 2)', dataC(:, :, 3)', '-', ...
                                    plotArgs{:}, p.Results.plotOptions{:});
                            end
                        end
                    else
                        % plot each trial from this condition individually
                        hData{iCond, iAlign} = nan(nTrialsC, 1);
                        for iTrial = 1:nTrialsC
                            if D == 1
                                tvec = timeC{iTrial};
                                dvec = dataC{iTrial};
                                tOffset = timeOffsetByAlign(iAlign);

                                if ischar(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                                    xvec = dvec;
                                    yvec = tvec;
                                else
                                    yvec = dvec;
                                    xvec = tvec;
                                end                         

                                if ~isempty(tvec) && ~isempty(dvec)
                                    if p.Results.alpha < 1
                                       hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(xvec + tOffset, yvec, ...
                                           'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                           'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                                    else
                                        plotArgs = app(iCond).getPlotArgs();
                                        hData{iCond, iAlign}(iTrial) = plot(axh, xvec, yvec, '-', ...
                                            plotArgs{:}, p.Results.plotOptions{:});
                                    end
                                end

                            elseif D==2
                                dxvec = dataC{iTrial}(:, 1);
                                dyvec = dataC{iTrial}(:, 2);

                                if ~isempty(dxvec) && ~isempty(dyvec)
                                    if p.Results.alpha < 1
                                       hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(dxvec, dyvec, ...
                                           'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                           'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                                    else
                                        plotArgs = app(iCond).getPlotArgs();
                                        hData{iCond, iAlign}(iTrial) = plot(axh, dxvec, dyvec, '-', ...
                                            plotArgs{:}, p.Results.plotOptions{:});
                                    end
                                end
                                
                            elseif D==3
                                dxvec = dataC{iTrial}(:, 1);
                                dyvec = dataC{iTrial}(:, 2);
                                dzvec = dataC{iTrial}(:, 3);

                                if ~isempty(dxvec) && ~isempty(dyvec) && ~isempty(dzvec)
                                    if p.Results.alpha < 1
                                       hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline3(dxvec, dyvec, dzvec, ...
                                           'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                           'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                                    else
                                        plotArgs = app(iCond).getPlotArgs();
                                        hData{iCond, iAlign}(iTrial) = plot3(axh, dxvec, dyvec, dzvec, '-', ...
                                            plotArgs{:}, p.Results.plotOptions{:});
                                    end
                                end
                            end

                            if iTrial == 1
                                hold(axh, 'on'); 
                            end
                        end
                    end
                    
                    % update name for inclusion in legend
                    TrialDataUtilities.Plotting.showFirstInLegend(hData{iCond, iAlign}, td.conditionNames{iCond});
                end
            end
            
            for iAlign = 1:nAlignUsed
                for iCond = 1:nConditionsUsed
                    timeC = time{iCond, iAlign};
                    dataC = data{iCond, iAlign};
                    if(D == 1)
                        % draw marks and intervals on timeseries
                        td.alignInfoActive.drawOnDataByTrial('time', timeC, 'data', dataC, ...
                            'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                            'showInLegend', iCond == 1, ...
                            'axh', axh, 'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize);

                    elseif D == 2
                        % draw marks and intervals on timeseries
                        td.alignInfoActive.drawOnDataByTrial('time', timeC, 'data', dataC, ...
                            'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                            'showInLegend', iCond == 1, ...
                            'axh', axh, 'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize);

                    elseif D == 3
                        % draw marks and intervals on timeseries
                        td.alignInfoActive.drawOnDataByTrial('time', timeC, 'data', dataC, ...
                            'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                            'showInLegend', iCond == 1, ...
                            'axh', axh, 'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize);
                    end
                end

                if D == 1
                    if ischar(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                        % x is data, y is time
                        td.alignSummarySet{iAlign}.setupTimeAutoAxis('which', 'y', 'style', p.Results.timeAxisStyle);
                        if ~isempty(p.Results.axisInfoX)
                            TrialDataUtilities.Plotting.setupAxisForChannel(td.channelDescriptorsByName.(name), 'which', 'x');
                        end
                    else
                        % y is data, x is time
                        td.alignSummarySet{iAlign}.setupTimeAutoAxis('which', 'x', 'style', p.Results.timeAxisStyle);
                        if ~isempty(p.Results.axisInfoY)
                            TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, 'which', 'y');
                        end
                    end
                end
            end

            box(axh, 'off');
            axis(axh, 'tight');
        end
        
        function plotAnalogGroupedEachTrial(td, name, varargin) 
            % grab raw data (for marking) and grouped data (for plotting)
            [dataByGroup, timeByGroup] = td.getAnalogGrouped(name);     
            td.plotProvidedAnalogDataGroupedEachTrial(1, 'time', timeByGroup(:), ...
                'data', dataByGroup(:), 'axisInfoY', td.channelDescriptorsByName.(name), ...
                varargin{:});
        end
        
        function plotAnalogGroupedEachTrial2D(td, name1, name2, varargin) 
            [dataCell, tvec] = td.getMultiAnalogAsMatrixGrouped({name1, name2});
            td.plotProvidedAnalogDataGroupedEachTrial(2, ...
                'time', tvec, 'data', dataCell(:), ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), varargin{:});
        end
        
        function plotAnalogGroupedEachTrial3D(td, name1, name2, name3, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            dataByGroup1 = td.getAnalogGrouped(name1);  
            dataByGroup2 = td.getAnalogGrouped(name2);
            dataByGroup3 = td.getAnalogGrouped(name3);
            app = td.conditionAppearances;

            for iCond = 1:td.nConditions
                dataCell1 = dataByGroup1{iCond};
                dataCell2 = dataByGroup2{iCond};
                dataCell3 = dataByGroup3{iCond};
                for iTrial = 1:numel(dataCell1)
                    plot3(axh, dataCell1{iTrial}, dataCell2{iTrial}, dataCell3{iTrial}, ...
                        '-', 'Color', app(iCond).color, ...
                        'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                    if iTrial == 1, hold(axh, 'on'); end
                end
            end
            box(axh, 'off');
            axis(axh, 'tight');
            axis(axh, 'vis3d');
            
            xlabel(td.getAxisLabelForChannel(name1));
            ylabel(td.getAxisLabelForChannel(name2));
            zlabel(td.getAxisLabelForChannel(name3));
        end
        
        function plotAnalogGroupMeans(td, name, varargin) 
            % plot the mean and sem for an analog channel vs. time within
            % each condition
            import TrialDataUtilities.Plotting.errorshade;
            
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.addParamValue('minTrials', 1, @isscalar);
            p.addParamValue('alpha', 1, @isscalar);
            p.addParamValue('markAlpha', 1, @isscalar);
            p.addParamValue('timeAxisStyle', 'tickBridge', @ischar);
            p.addParamValue('markShowRanges', true, @islogical); % show ranges for marks 
           % p.KeepUnmatched = true;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);
            cla(axh);
            app = td.conditionAppearances;

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getAnalogGroupMeans(name, 'minTrials', p.Results.minTrials);     
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign});
            end
            
            % determine the x-offsets at which to plot the data
            tOffsets = td.getAlignPlottingTimeOffsets(tvecCell);
            
            for iAlign = 1:td.nAlign
                for iCond = 1:td.nConditions
                    h = errorshade(tvecCell{iAlign} + tOffsets(iAlign), ...
                        meanMat{iAlign}(iCond, :), ...
                        semMat{iAlign}(iCond, :), app(iCond).Color, 'axh', axh, ...
                        'alpha', p.Results.alpha);
                    hold(axh, 'on');
                    if ~isnan(h)
                        set(h, 'DisplayName', td.conditionNames{iCond});
                    end
                end
                box(axh, 'off');
                axis(axh, 'tight');
                
                td.alignSummarySet{iAlign}.drawOnTimeseriesByCondition(tvecCell{iAlign}, ...
                    permute(meanMat{iAlign}, [2 3 1 4]), ...  % data needs to be T x D x C x N, currently C x T
                    'tOffsetZero', tOffsets(iAlign), 'alpha', p.Results.alpha, ...
                    'markAlpha', p.Results.markAlpha, ...
                    'tMin', min(tvecCell{iAlign}), 'tMax', max(tvecCell{iAlign}));
                
                % setup x axis 
                td.alignSummarySet{iAlign}.setupTimeAutoAxis('tOffsetZero', tOffsets(iAlign), ...
                    'tMin', min(tvecCell{iAlign}), 'tMax', max(tvecCell{iAlign}), ...
                    'style', p.Results.timeAxisStyle);
            end

            % setup y axis
            TrialDataUtilities.Plotting.setupAxisForChannel(td.channelDescriptorsByName.(name));
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
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.addParamValue('legend', false, @islogical);
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
                        tdTheseCond.alignSummarySet{iAlign}.setupTimeAutoAxis('tOffsetZero', tOffsets(iAlign), 'style', 'marker');
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
end
