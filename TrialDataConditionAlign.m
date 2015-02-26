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
        conditionsAxisAttributesOnly
        conditionsAsStrings
        conditionNames
        conditionNamesShort
        conditionNamesMultiline
        conditionAppearances
        conditionColors
        conditionsSize
        
        nAxes
        axisValueLists
        axisValueListsAsStrings
        axisValueListsAsStringsShort
        axisNames
        axisDescriptions
        
        conditionAppearanceFn
        conditionNameFn
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
            td = td.rebuildOnDemandCache();
            td = td.initializeConditionInfo();
            td = td.initializeAlignInfo();
            td = td.updateValid();
        end
        
        function td = rebuildOnDemandCache(td)
            td.warnIfNoArgOut(nargout);
            td.odc = TrialDataConditionAlignOnDemandCache(); 
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
    methods(Access=protected)
        function buildEventData(td)
            % build a cached version of event times into a matrix for easy alignment
            % sets .eventData and .eventCounts
            
            evStruct = td.getRawEventFlatStruct();
            evList = fieldnames(evStruct);
            nEvents = numel(evList);
            eventCounts = struct();
            eventData = struct();
            
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
    end
    
    % General utilites
    methods
         % synchronize valid between AlignInfo and ConditionINfo
         % shouldn't need to call this manually, but just in case
        function td = updateValid(td)
            td.warnIfNoArgOut(nargout);
            td = updateValid@TrialData(td);

            if isempty(td.manualValid)
                td.manualValid = truevec(td.nTrials);
            end

            valid = td.buildValid();

            td.conditionInfo = td.conditionInfo.setInvalid(~valid);
            for iA = 1:td.nAlign
                td.alignInfoSet{iA} = td.alignInfoSet{iA}.setInvalid(~valid);
            end
            
            % cause align summary to be recomputed
            td.alignSummarySet = [];
        end
        
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
                mask = ismember(alignEvents, names);
                if any(mask)
                    error('TrialData alignment depends on event %s', ...
                        strjoin(alignEvents(mask)));
                end
            end
            
            % check whether any of the events are in condition info
            conditionParams = td.conditionInfo.attributeNames;
            mask = ismember(conditionParams, names);
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
        
        % given a cellvec or nmeric vector, group its elements along
        % dimension 1 but without preserving the conditions tensor shape
        function varargout = groupElementsFlat(td, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');  
                varargout{i} = cellfun(@(idx) data(idx,:), td.listByCondition(:), ...
                    'UniformOutput', false);
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
            if isempty(cd.attributeRequestAs)
                paramData = emptyStructArray(td.nTrials);
            else
                paramData = td.getRawChannelDataAsStruct(cd.attributeRequestAs);
            end
            
            % build condition info from condition descriptor
            td.conditionInfo = ConditionInfo.fromConditionDescriptor(cd, paramData);
            td = td.postUpdateConditionInfo();
        end
        
        function td = setConditionAppearanceFn(td, fn)
            % Update the appearanceFn callback of conditionDescriptor
            % without invalidating any of the other cached info
            assert(isa(fn, 'function_handle'));
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.appearanceFn = fn;
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
            
            if ischar(names)
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
                
                cd = td.channelDescriptorsByName.(name);
                if isa(cd, 'EventChannelDescriptor')
                    % the event time will be referenced from the current
                    % zero event
                    nameMod = matlab.lang.makeValidName(sprintf('%s_from_%s', name, td.alignInfoActive.zeroLabel));
                    values = td.getEventFirst(name);
                    
                elseif isa(cd, 'AnalogChannelDescriptor')
                    % take analog sample at zero time
                    nameMod = matlab.lang.makeValidName(sprintf('%s_at_%s', name, td.alignInfoActive.zeroLabel));
                    values = td.getAnalogSample(name);
                    
                elseif isa(cd, 'ParamChannelDescriptor')
                    % params get used as is
                    nameMod = name;
                    values = td.getParamRaw(name);
                else
                    error('Unable to create attribute for channel type %s', class(cd));
                end
                
                td.conditionInfo = td.conditionInfo.addAttribute(nameMod, 'values', values);
                namesModified{i} = nameMod;
            end
            
            if wasChar
                namesModified = namesModified{1};
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
        
        function td = addAxis(td, attrList, varargin)
            td.warnIfNoArgOut(nargout);
            
            % add any needed attributes to condition info
            td = td.addAttribute(attrList);
            
            td.conditionInfo = td.conditionInfo.addAxis(attrList, varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = ungroup(td) 
            % this only undoes the grouping axes, NOT the value list
            % filtering. use reset condition info for that
            td.warnIfNoArgOut(nargout);
            td = td.groupBy();
        end
        
        function td = reset(td)
            td.warnIfNoArgOut(nargout);
            td = td.resetConditionInfo();
            td = td.unalign();
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
        
        function td = flattenAxesExceptFirst(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.flattenAxesExceptFirst(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = sortWithinConditionsBy(td, attrList, varargin)
            td.warnIfNoArgOut(nargout);
            
            if ischar(attrList)
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
        
        % Axis randomization settings
        function td = setRandomSeed(td, varargin)
            % args: seed
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setRandomSeed(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = newRandomSeed(td, varargin)
            % no args
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.newRandomSeed(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = noRandomization(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.axisNoRandomization(varargin{:});
            td = td.postUpdateConditionInfo();
        end  
        
        function td = resampleTrialsWithinConditions(td, varargin)
            % no args
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.resampleTrialsWithinConditions(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = axisNoRandomization(td, varargin)
            % args: idxOrAttr
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.axisNoRandomization(varargin{:});
            td = td.postUpdateConditionInfo();
        end
                   
        function td = axisShuffle(td, varargin)
            % args (idxOrAttr, replace)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.axisShuffle(varargin{:});
            td = td.postUpdateConditionInfo();
        end

        function td = axisResampleFromSpecifiedValueListIndices(td, varargin)
            % args: (axisIdxOrAttr, resampleFromIndices, replace) 
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.axisResampleFromSpecifiedValueListIndices(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = axisResampleFromSpecifiedValues(td, varargin)
            % args: (axisIdxOrAttr, resampleFromValueStructMatch, replace) 
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.axisResampleFromSpecifiedValues(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = postUpdateConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
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
        
        function v = get.conditionsAxisAttributesOnly(td)
            v = td.conditionInfo.conditionsAxisAttributesOnly;
        end
        
        function v = get.conditionsAsStrings(td)
            v = td.conditionInfo.conditionsAsStrings;
        end

        function n = get.nConditions(td)
            n = td.conditionInfo.nConditions;
        end
        
        function sz = get.conditionsSize(td)
            sz = td.conditionInfo.conditionsSize;
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
            td = td.align('TrialStart:TrialEnd @ TimeZero');
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
        
        function td = lag(td, varargin)
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
        
        function [data, time] = getAnalogEachAlign(td, name, varargin)
            % data, time are nTrials x nAlign cells
            [data, time] = deal(cell(td.nTrials, td.nAlign));
            for iA = 1:td.nAlign
                [data(:, iA), time(:, iA)] = td.useAlign(iA).getAnalog(name);
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
            timeCell = cell(td.nTrials, C);
            for c = 1:C
                [~, timeCell(:, c)] = td.getAnalog(nameCell{c});
            end
            
            tvec = TrialDataUtilies.Data.inferCommonTimeVectorForTimeseriesData(timeCell, ...
                'timeDelta', timeDelta);
        end
        
        function [dataUnif, timeUnif, delta] = getAnalogUniformlySampled(td, name, varargin)
            [dataUnif, timeUnif, delta] = getAnalogUniformlySampled@TrialData(td, name, varargin{:});
            [dataUnif, timeUnif] = td.alignInfoActive.getAlignedTimeseries(dataUnif, timeUnif, false);
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
            % case, use getMultiAnalogAsMatrix
            
            p = inputParser;
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('tvec', [], @isvector);
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
        

        function [data, tvec] = getMultiAnalogAsMatrix(td, name, varargin)
            p = inputParser;
            p.addParameter('timeDelta', [], @isscalar);
            p.addParameter('tvec', [], @isvector);
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
        
        function [mat, tvec, alignIdx] = getMultiAnalogAsMatrixEachAlign(td, name, varargin)
            % similar to getMultiAnalogAsMatrix, except each alignment will be
            % concatenated in time
            
            matCell = cellvec(td.nAlign);
            tvecCell = cellvec(td.nAlign);
            for iA = 1:td.nAlign
                [matCell{iA}, tvecCell{iA}] = td.useAlign(iA).getMultiAnalogAsMatrix(name, varargin{:});
            end
            
            [mat, alignIdx] = TensorUtils.catWhich(2, matCell{:});
            tvec = cat(2, tvecCell{:});
        end
        
        function [dCell, tCell] = getAnalogGrouped(td, name)
            [dataCell, timeCell] = td.getAnalog(name);
            [dCell, tCell] = td.groupElements(dataCell, timeCell);
        end
        
        function [dCell, tCell] = getAnalogGroupedEachAlign(td, name, varargin)
            [dataCell, timeCell] = td.getAnalogEachAlign(name, varargin{:});
            [dCell, tCell] = deal(cell(td.nConditions, td.nAlign));
            for iA = 1:td.nAlign
                [dCell(:, iA), tCell(:, iA)] = td.useAlign(iA).groupElementsFlat(dataCell(:, iA), timeCell(:, iA));
            end
        end
        
        function [dCell, tvec] = getAnalogAsMatrixGrouped(td, name)
            [mat, tvec] = td.getAnalogAsMatrix(name);
            dCell = td.groupElements(mat);
        end
        
        function [dCell, tvec, alignIdx] = getAnalogAsMatrixGroupedEachAlign(td, name)
            [mat, tvec, alignIdx] = td.getAnalogAsMatrixEachAlign(name);
            dCell = td.groupElements(mat);
        end
        
        function [dataCell, tvec] = getMultiAnalogAsMatrixGrouped(td, nameCell, varargin)
            % dataCell will be size(td.conditions)
            % contents will be nTrials x T x nChannels
            [data, tvec] = td.getMultiAnalogAsMatrix(nameCell);
            dataCell = td.groupElements(data);
        end
        
        function [dCell, tvec, alignIdx] = getMultiAnalogAsMatrixGroupedEachAlign(td, nameCell)
            [mat, tvec, alignIdx] = td.getMultiAnalogAsMatrixEachAlign(nameCell);
            dCell = td.groupElements(mat);
        end
        
        function [meanMat, semMat, tvec, stdMat, nTrialsMat] = getAnalogGroupMeans(td, name, varargin)
            % *Mat will be nConditions x T matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
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
        
        function [meanMat, semMat, tvec, stdMat, nTrialsMat, alignIdx] = getAnalogGroupMeansEachAlign(td, name, varargin)
            % *Mat will be nConditions x T x nChannels matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            if iscell(name)
                nC = numel(name);
            else
                nC = 1;
            end
            [dCell, tvec, alignIdx] = td.getAnalogAsMatrixGroupedEachAlign(name);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nC));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), stdMat(iC, :, :)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials);
                end
            end 
        end
        
         function [meanMat, semMat, tvec, nTrialsMat, stdMat] = ...
                 getMultiAnalogGroupMeans(td, nameCell, varargin)
             % *Mat will be nConditions x T x nChannels tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            assert(iscellstr(nameCell), 'Must be cellstr of channel names');
            nChannels = numel(nameCell);
            
            % dCell is size(conditions) with nTrials x T x nChannels inside
            [dCell, tvec] = td.getMultiAnalogAsMatrixGrouped(nameCell);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nChannels));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), ...
                        stdMat(iC, :, :)] = nanMeanSemMinCount(dCell{iC}, 1, minTrials);
                end
            end 
         end
         
         function [meanMat, semMat, tvec, nTrialsMat, stdMat, alignIdx] = ...
                 getMultiAnalogGroupMeansEachAlign(td, nameCell, varargin)
             % *Mat will be nConditions x T x nChannels tensors
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            assert(iscellstr(nameCell), 'Must be cellstr of channel names');
            nChannels = numel(nameCell);
            
            % dCell is size(conditions) with nTrials x T x nChannels inside
            [dCell, tvec, alignIdx] = td.getMultiAnalogAsMatrixGroupedEachAlign(nameCell);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.nConditions, numel(tvec), nChannels));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC, :, :), semMat(iC, :, :), nTrialsMat(iC, :, :), ...
                        stdMat(iC, :, :)] = nanMeanSemMinCount(dCell{iC}, 1, minTrials);
                end
            end 
         end
        
         function [dataVec, timeVec] = getAnalogSample(td, name, varargin)
             % same as in TrialData, except issues warning if alignment has
             % more than one sample. Pulls a single time point from the t=0
             % aligned time for each trial. Guaranteed to be scalar per
             % trial, i.e. vector over trials
             if ~td.alignInfoActive.isStartStopEqual
                 warning('Getting analog sample when alignment is not for single timepoint');
             end
             [dataVec, timeVec] = getAnalogSample@TrialData(td, name, varargin{:});
        end 
        
        function [dataCell, timeCell] = getAnalogSampleGrouped(td, name, varargin)
            [dataVec, timeVec] = td.getAnalogSample(name);
            [dataCell, timeCell] = td.groupElements(dataVec, timeVec);
        end
        
        % differentiation
        function [diffData, time, diffUnits] = differentiateAnalogChannel(td, name, varargin)
            % data will be in units / second, checking .timeUnitsPerSecond
            % in order to normalize appropriately
            p = inputParser;
            p.addParameter('delta', [], @isscalar);
            p.addParameter('interpolationMethod', 'linear', @ischar);
            p.addParameter('smoothing', 7, @(x) isscalar(x) && mod(x, 2) == 1);
            p.addParameter('order', 1, @isscalar);
            p.addParameter('polynomialOrder', 2, @isscalar);
            p.parse(varargin{:});
            
            [data, time, delta] = td.getAnalogUniformlySampled(name, ...
                'delta', p.Results.delta, 'method', p.Results.interpolationMethod);
            diffData = cellvec(td.nTrials);
            
            w = -1 / (delta / td.timeUnitsPerSecond) ^ p.Results.order;
            prog = ProgressBar(td.nTrials, 'Smoothing/Differentiating %s', name);
            for iT = 1:td.nTrials
                prog.update(iT);
                if isempty(data{iT})
                    continue;
                end
                
                diffData{iT} = w * TrialDataUtilities.Data.savitzkyGolayFilt( ...
                    data{iT}, p.Results.polynomialOrder, p.Results.order, p.Results.smoothing)'; 
            end
            prog.finish();
            
            diffUnits = sprintf('%s/sec', td.getChannelUnitsPrimary(name));
        end
        
        function td = addDifferentiatedAnalogChannel(td, name, diffName, varargin)
            % drops alignment and grouping before differentiating since
            % this is typically what the user wants
            
            td.warnIfNoArgOut(nargout);
            [diffData, time, diffUnits] = td.reset().differentiateAnalogChannel(name, varargin{:});
            td = td.addAnalog(diffName, diffData, time, 'units', diffUnits, 'isAligned', false);
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
            % get parameter values grouped by condition
            % dCell will be size(conditions) cell tensor with value arrays
            % within
            dCell = td.groupElements(td.getParam(name));
        end
        
        function values = getParamUniqueGrouped(td, name)
            vCell = td.getParamGrouped(name);
            values = cellfun(@getUnique, vCell, 'UniformOutput', false);
            
            function values = getUnique(vals)
                if ~iscell(vals)
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
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            dCell = td.getParamGrouped(name);
            [meanMat, semMat, nTrialsMat, stdMat] = deal(nan(td.conditionsSize));
            for iC = 1:td.nConditions
                if ~isempty(dCell{iC})
                    [meanMat(iC), semMat(iC), nTrialsMat(iC), stdMat(iC)] = ...
                        nanMeanSemMinCount(dCell{iC}, 1, minTrials);
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
            cdX = td.channelDescriptorsByName.(nameX);
            cdY = td.channelDescriptorsByName.(nameY);
            td.plotProvidedGroupedScatterData(dataX, dataY, ...
                'axisInfoX', cdX, 'axisInfoY', cdY, varargin{:});
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
            %p.addOptional('dataZ', [], @(x) isnumeric(x) || iscell(x)); 
            
            p.addParameter('axisInfoX', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            %p.addParameter('axisInfoZ', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            
            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @ischar);
            p.addParameter('axisStyleY', 'tickBridge', @ischar);
            
            p.addParameter('conditionIdx', 1:td.nConditions, @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('alpha', 0.7, @isscalar);
            p.addParameter('edgeAlpha', [], @isscalar);
            p.addParameter('markerSize', 5, @isscalar);
            p.addParameter('useThreeVector', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            alpha = p.Results.alpha;
            if isempty(p.Results.edgeAlpha)
                edgeAlpha = alpha;
            else
                edgeAlpha = p.Results.edgeAlpha;
            end
            
            axh = td.getRequestedPlotAxis(p.Unmatched);
            hold(axh, 'on');
            
            if ~ismember('scaleBars', p.UsingDefaults) && p.Results.scaleBars
                axisStyleX = 'scaleBar';
                axisStyleY = 'scaleBar';
            else
                axisStyleX = p.Results.axisStyleX;
                axisStyleY = p.Results.axisStyleY;
            end
            
            conditionIdx = p.Results.conditionIdx;
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
            
            h = nan(nConditionsUsed, 1);
            for iC = 1:nConditionsUsed
                idxC = conditionIdx(iC);
                args = app(idxC).getMarkerPlotArgs();
                if isempty(dataX{iC}) || isempty(dataY{iC})
                    continue;
                end
                h(iC) = plot(dataX{iC}, dataY{iC}, 'o', 'MarkerSize', p.Results.markerSize, ...
                    args{:}, p.Results.plotOptions{:});
                SaveFigure.setMarkerOpacity(h(iC), alpha, edgeAlpha);
                
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
        end
        
        function plotParamVsParamMeanLinePlot(td, nameX, nameY, varargin)
            % for each bin or grouping of nameX, plot the mean of paramY
            % as a line plot with vertical error bars conveying either sem
            % or std (according to 'errorType' parameter). One line plot
            % for each group defined in td (before adding nameX to the
            % list).
            
            p = inputParser();
            p.addParameter('errorType', 'sem', @ischar);
            p.addParameter('LineWidth', 1, @isscalar);
            p.KeepUnmatched = true;
            p.CaseSensitive = false;
            p.parse(varargin{:});
            
            axh = td.getRequestedPlotAxis(p.Unmatched);
            hold(axh, 'on');
            
            % cache before messing with axis shapes
            conditionNames = td.conditionNamesShort;
            
            % add as the last axis
            td = td.addAxis(nameX);
            
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
            h = nan(nCondOther, 1);
            for iC = 1:nCondOther
                h(iC) = TrialDataUtilities.Plotting.errorline(xValues, meanY(:, iC), errorY(:, iC), ...
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

    % Spike data
    methods
        % return aligned unit spike times
        function [timesCell] = getSpikeTimes(td, unitName, includePadding)
            if nargin < 3
                includePadding = false;
            end
            timesCell = getSpikeTimes@TrialData(td, unitName);
            timesCell = td.alignInfoActive.getAlignedTimesCell(timesCell, includePadding);
        end
        
        function [rateCell, timeCell, hasSpikes] = getSpikeRateFiltered(td, unitName, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.addParameter('removeZeroSpikeTrials', true, @islogical);
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            spikeCell = td.getSpikeTimes(unitName, true);
            timeInfo = td.alignInfoActive.timeInfo;
            
            % provide an indication as to which trials have spikes
            hasSpikes = ~cellfun(@isempty, spikeCell);
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond);
        end
           
        function [rates, tvec, hasSpikes] = getSpikeRateFilteredAsMatrix(td, unitName, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            if isempty(sf)
                sf = SpikeFilter.getDefaultFilter();
            end
            timeDelta = p.Results.timeDelta;
            if isempty(timeDelta)
                timeDelta = td.alignInfoActive.minTimeDelta;
                if isempty(timeDelta)
                    timeDelta = 1;
                    warning('Using timeDelta=%d for spike rate timepoints. Specify ''timeDelta'' or call .round for consistent results', timeDelta);
                end
            end
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            % critical to include the spike times in the padded window
            spikeCell = td.getSpikeTimes(unitName, true);
            [tMinByTrial, tMaxByTrial] = td.alignInfoActive.getStartStopRelativeToZeroByTrial();
            
            % provide an indication as to which trials have spikes
            hasSpikes = ~cellfun(@isempty, spikeCell);
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            [rates, tvec] = sf.filterSpikeTrainsWindowByTrialAsMatrix(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond, ...
                'timeDelta', timeDelta);
            tvec = makecol(tvec);
        end
        
        function [rates, tvec, hasSpikes, alignVec] = getSpikeRateFilteredAsMatrixEachAlign(td, unitName, varargin)
            ratesCell = cellvec(td.nAlign);
            tvecCell = cellvec(td.nAlign);
            hasSpikesMat = nan(td.nTrials, td.nAlign);
            for iA = 1:td.nAlign
                [ratesCell{iA}, tvecCell{iA}, hasSpikesMat(:, iA)] = td.useAlign(iA).getSpikeRateFilteredAsMatrix(unitName);
            end
            
            [rates, alignVec] = TensorUtils.catWhich(2, ratesCell{:});
            tvec = cat(1, tvecCell{:});
            hasSpikes = any(hasSpikesMat, 2);
        end

        function [rateCell, timeCell, hasSpikesGrouped] = getSpikeRateFilteredGrouped(td, unitName, varargin)
            [rateCell, timeCell] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElements(rateCell);
            timeCell = td.groupElements(timeCell);
            hasSpikesGrouped = td.groupElements(hasSpikes);
        end
        
        function [rateCell, tvec, hasSpikesGrouped] = getSpikeRateFilteredAsMatrixGrouped(td, unitName, varargin)
            [rates, tvec, hasSpikes] = td.getSpikeRateFilteredAsMatrix(unitName, varargin{:});
            rateCell = td.groupElements(rates);
            hasSpikesGrouped = td.groupElements(hasSpikes);
        end
        
        function [rateCell, tvec, hasSpikesGrouped, alignVec] = getSpikeRateFilteredAsMatrixGroupedEachAlign(td, unitName, varargin)
            [rates, tvec, hasSpikes, alignVec] = td.getSpikeRateFilteredAsMatrixEachAlign(unitName, varargin{:});
            rateCell = td.groupElements(rates);
            hasSpikesGrouped = td.groupElements(hasSpikes);
        end
        
        function [psthMat, tvec, semMat, stdMat, nTrialsMat] = ...
                getSpikeRateFilteredMeanByGroup(td, unitName, varargin)
            % *Mat will be nConditions x T matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('removeZeroSpikeTrials', true, @islogical);
            p.parse(varargin{:});
            
            if isempty(p.Results.minTrials)
                minTrials = 1;
            else
                minTrials = p.Results.minTrials;
            end
            
            if isempty(p.Results.minTrialFraction)
                minTrialFraction = 1;
            else
                minTrialFraction = p.Results.minTrialFraction;
            end
            
            [rateCell, tvec, hasSpikesGrouped] = td.getSpikeRateFilteredAsMatrixGrouped(unitName, ...
                'timeDelta', p.Results.timeDelta, 'spikeFilter', p.Results.spikeFilter);
            
            % remove trials from each group that have no spikes
            if p.Results.removeZeroSpikeTrials
                for iC = 1:td.nConditions
                    rateCell{iC} = rateCell{iC}(hasSpikesGrouped{iC}, :);
                end
            end
            
            % comute the means
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec)));
            for iC = 1:td.nConditions
                if ~isempty(rateCell{iC})
                    [psthMat(iC, :), semMat(iC, :), nTrialsMat(iC, :), stdMat(iC, :)] = ...
                        nanMeanSemMinCount(rateCell{iC}, 1, minTrials, minTrialFraction);
                end
            end
        end
        
        function [psthMat, tvec, semMat, stdMat, nTrialsMat] = ...
                getSpikeRateFilteredMeanByGroupEachAlign(td, unitName, varargin)
            % *Mat will be nConditions x T matrices
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar); % minimum trial count to average
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('removeZeroSpikeTrials', true, @islogical);
            p.parse(varargin{:});
            minTrials = p.Results.minTrials;
            
            [rateCell, tvec, hasSpikesGrouped] = td.getSpikeRateFilteredAsMatrixGroupedEachAlign(unitName, ...
                'timeDelta', p.Results.timeDelta, 'spikeFilter', p.Results.spikeFilter);
            
            % remove trials from each group that have no spikes
            if p.Results.removeZeroSpikeTrials
                for iC = 1:td.nConditions
                    rateCell{iC} = rateCell{iC}(hasSpikesGrouped{iC}, :);
                end
            end
            
            % comute the means
            [psthMat, semMat, stdMat, nTrialsMat] = deal(nan(td.nConditions, numel(tvec)));
            for iC = 1:td.nConditions
                if ~isempty(rateCell{iC})
                    [psthMat(iC, :), semMat(iC, :), nTrialsMat(iC, :), stdMat(iC, :)] = ...
                        nanMeanSemMinCount(rateCell{iC}, 1, minTrials);
                end
            end
        end
            
        function [psthMat, tvec, semMat, stdMat, nTrialsMat] = getPSTH(td, varargin)
            [psthMat, tvec, semMat, stdMat, nTrialsMat] = getSpikeRateFilteredMeanByGroup(td, varargin{:});
        end
        
        function timesCellofCells = getSpikeTimesGrouped(td, unitName, includePadding)
            if nargin < 3
                includePadding = false;
            end
            timesCell = td.getSpikeTimes(unitName, includePadding);
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
        
        function [meanByGroup, semByGroup, stdByGroup, nByGroup] = getSpikeMeanRateGroupMeans(td, unitName) 
            rateCell = td.getSpikeMeanRateGrouped(unitName);
            
            meanByGroup = cellfun(@nanmean, rateCell);
            semByGroup = cellfun(@nansem, rateCell);
            stdByGroup = cellfun(@nanstd, rateCell);
            nByGroup = cellfun(@numel, rateCell);
        end
        
        function plotTuningCurve(td, unitName, varargin)
            p = inputParser;
            p.addParameter('errorType', 'sem', @ischar);
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
            [timesCellByGroup, markCountsByGroup] = td.groupElements(timesCell, markCounts);
        end
        
        function [spikeCounts, markCounts] = getMarkAlignedSpikeCounts(td, unitName, markIdx, window)
            [timesCell, markCounts] = td.getMarkAlignedSpikeTimes(unitName, markIdx, window);
            spikeCounts = cellfun(@numel, timesCell);
            for iTrial = 1:td.nTrials
                spikeCounts(iTrial, markCounts(iTrial)+1:end) = NaN;
            end
        end
        
        function [spikeCountsByGroup, markCountsByGroup] = getMarkAlignedSpikeCountsGrouped(td, unitName, markIdx, window)
            [spikeCounts, markCounts] = td.getMarkAlignedSpikeCounts(unitName, markIdx, window);
            [spikeCountsByGroup, markCountsByGroup] = td.groupElements(spikeCounts, markCounts);
        end
    end
    
    methods % spike waveforms
        function [wavesCell, waveTvec, timesCell] = getSpikeWaveforms(td, unitName, varargin)
            [wavesCell, waveTvec, timesCell] = getSpikeWaveforms@TrialData(td, unitName);
            [~, maskCell] = td.alignInfoActive.getAlignedTimesCell(timesCell, true);
            wavesCell = cellfun(@(waves, mask) waves(mask, :), wavesCell, maskCell, 'UniformOutput', false);
            timesCell = cellfun(@(times, mask) times(mask), timesCell, maskCell, 'UniformOutput', false);
        end
        
        function h = plotSpikeWaveforms(td, unitName, varargin)
            p = inputParser();
            p.addParameter('maxToPlot', 200, @isscalar);
            p.addParameter('alpha', 0.4, @isscalar);
            p.addParameter('color', 'k', @(x) ischar(x) || isvector(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [axh, unmatched] = td.getRequestedPlotAxis(p.Unmatched);
            
            [wavesCell, waveTvec] = td.getSpikeWaveforms(unitName, unmatched);
            wavesMat = cat(1, wavesCell{:});
            
            maxSamples = p.Results.maxToPlot;
            if maxSamples < size(wavesMat, 1)
                s = RandStream('mt19937ar','Seed',1);
                idx = randsample(s, size(wavesMat, 1), maxSamples);
                wavesMat = wavesMat(idx, :);
            end
            
            h = nanvec(size(wavesMat, 1));
            for i = 1:size(wavesMat, 1)
                h(i) = TrialDataUtilities.Plotting.patchline(waveTvec, wavesMat(i, :), ...
                    'Parent', axh, 'EdgeColor', p.Results.color, 'EdgeAlpha', p.Results.alpha);
            end
            
            ht = title(sprintf('Spike Waveforms for %s', unitName));
            set(ht, 'Interpreter', 'none');
            
            axis tight;
            AutoAxis.replaceScaleBars(gca, td.timeUnitName, td.channelDescriptorsByName.(unitName).waveformsUnits);
            
        end
        
        function [wavesMat, timeWithinTrial, trialIdx, whichUnit] = getSpikeWaveformMatrix(td, units, varargin)
            % returns a matrix containing all waveforms in units. Units
            % can be a string or cellstr, where each string is a spike unit name
            % or, if paramValue 'regexp' is true, a regexp search over
            % unitNames
            p = inputParser();
            p.addParameter('regexp', false, @islogical);
            p.parse(varargin{:});
            
            if ischar(units), units = {units}; end
            
            % do regexp matching
            if p.Results.regexp
                allUnits = td.listSpikeChannels();
                matches = cellvec(numel(units));
                for iU = 1:numel(units)
                    matchRes = regexp(allUnits, units{iU}, 'match');
                    matches{iU} = cat(1, matchRes{:});
                end
                units = cat(1, matches{:});
            end
            
            % get waveforms for each unit
            [waveByU, trialByU, timeByU] = deal(cellvec(numel(units)));
            for iU = 1:numel(units)
                [w, ~, t] = td.getSpikeWaveforms(units{iU});
                [waveByU{iU}, trialByU{iU}] = TensorUtils.catWhich(1, w{:});
                timeByU{iU} = cat(1, t{:});
            end
            
            [wavesMat, whichUnit] = TensorUtils.catWhich(1, waveByU{:});
            trialIdx = cat(1, trialByU{:});
            timeWithinTrial = cat(1, timeByU{:});
        end
    end
    
    % Plotting Spike data
    methods
        function plotPSTH(td, unitName, varargin)
            import TrialDataUtilities.Data.nanMeanSemMinCount;
            p = inputParser();
            p.addParameter('minTrials', [], @(x) isempty(x) || isscalar(x)); % minimum trial count to average
            p.addParameter('minTrialFraction', [], @(x) isempty(x) || isscalar(x)); % minimum trial fraction to average
            p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.addParameter('errorType', '', @(s) ismember(s, {'sem', 'std', ''}));
            p.addParameter('removeZeroSpikeTrials', true, @islogical);
            p.addParameter('axisMarginLeft', 2.5, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            axh = td.getRequestedPlotAxis(p.Unmatched);

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, tvecCell{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getSpikeRateFilteredMeanByGroup(unitName, ...
                    'minTrials', p.Results.minTrials, 'minTrialFraction', p.Results.minTrialFraction, ...
                    'timeDelta', p.Results.timeDelta, 'spikeFilter', p.Results.spikeFilter, ...
                    'removeZeroSpikeTrials', p.Results.removeZeroSpikeTrials);
            
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
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
            
            maskEmpty = cellfun(@isempty, tvecCell);
            if any(maskEmpty)
                error('No valid time window found for alignment. Perhaps there are not enough trials with spikes in all of the conditions? Try lowering minTrials or setting removeZeroSpikeTrials to false.');
            end
            
            td.plotProvidedAnalogDataGroupMeans(1, 'time', tvecCell, ...
                'data', meanMat, 'dataError', errorMat, 'axh', axh, ...
                'axisInfoX', 'time', 'axisInfoY', struct('name', 'Firing Rate', 'units', 'spikes/sec'), ...
                p.Unmatched);
            
            TrialDataUtilities.Plotting.setTitleIfBlank(axh, '%s : Unit %s', td.datasetName, unitName);
            axis(axh, 'tight');
            
            ylabel('spikes/sec');
            au = AutoAxis(axh);
            %au.addAutoAxisY();
            au.axisMarginLeft = p.Results.axisMarginLeft;
            au.update();
            
            hold(axh, 'off');
        end
        
        function plotRaster(td, unitName, varargin)
            p = inputParser();
            p.addParameter('conditionIdx', 1:td.nConditions, @isvector);
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
            p.addParameter('colorSpikes', false, @islogical);
            p.addParameter('timeAxisStyle', 'marker', @ischar);
            p.addParameter('intervalAlpha', 0.5, @isscalar);
            p.addParameter('intervalMinWidth', NaN, @isscalar); % if specified, draws intervals at least this wide to ensure visibility
            p.addParameter('gapBetweenConditions', [], @(x) isempty(x) || isscalar(x));
            
            % for plotting interval and marks above each condition
            p.addParameter('annotateAboveEachCondition', false, @islogical);
            p.addParameter('annotationHeight', 2, @islogical);
            p.addParameter('annotateUsingFirstTrialEachCondition', false, @islogical);
            
            % make room for labels using AutoAxis
            p.addParameter('axisMarginLeft', 2.5, @isscalar);
           
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            if td.nTrialsValid == 0
                error('No valid trials found');
            end
            
            axh = td.getRequestedPlotAxis(p.Unmatched);
            
            conditionIdx = p.Results.conditionIdx;
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            nConditionsUsed = numel(conditionIdx);
            alignIdx = p.Results.alignIdx;
            nAlignUsed = numel(alignIdx);
            
            timesByAlign = cell(nAlignUsed, nConditionsUsed);
            
            % compute x-axis offsets for each align
            timePointsCell = cell(nAlignUsed, 1);
            [startData, stopData] = deal(cell(nAlignUsed, nConditionsUsed));
            for iAlign = 1:nAlignUsed
                idxAlign = alignIdx(iAlign);
                
                % get grouped spike times by alignment
                thisC = td.useAlign(idxAlign).getSpikeTimesGrouped(unitName);
                timesByAlign(iAlign, :) = thisC(:);
                
                % figure out time validity window for this alignment
                % TODO might want to update this for the selected
                % conditions only
                [start, stop] = td.alignInfoSet{idxAlign}.getStartStopRelativeToZeroByTrial();
                timePointsCell{iAlign} = [start; stop];
                
                % and store start/stop by trial in each align/condition
                % cell
                for iCond = 1:nConditionsUsed
                    idxCond = conditionIdx(iCond);
                    startData{iAlign, iCond} = start(td.listByCondition{idxCond});
                    stopData{iAlign, iCond} = stop(td.listByCondition{idxCond});
                end
            end
            tOffsetByAlign = td.getAlignPlottingTimeOffsets(timePointsCell);
           
            % default gap depends on presence of annotations between conditions
            if isempty(p.Results.gapBetweenConditions)
                if p.Results.annotateAboveEachCondition
                    gap = p.Results.annotationHeight + 3; % 2 above, 1 below
                else
                    gap = 3;
                end
            else
                gap = p.Results.gapBetweenConditions;
            end
                    
            % compute y-axis offsets for each condition
            trialCounts = cellfun(@numel, td.listByCondition(conditionIdx));
            yOffsetByCondition = zeros(nConditionsUsed, 1);
            yLimsByCondition = nan(2, nConditionsUsed);
            currentOffset = 0;
            lastTrialCount = 0;
            for iC = 1:nConditionsUsed
                yOffsetByCondition(iC) = currentOffset;
                yLimsByCondition(:, iC) = [currentOffset; currentOffset-trialCounts(iC)];
                if trialCounts(iC) > 0
                    currentOffset = currentOffset - trialCounts(iC) - gap;
                    lastTrialCount = trialCounts(iC);
                end
            end
            % shift to lie entirely above y = 0
            delta = - min(yOffsetByCondition) + lastTrialCount;
            yOffsetByCondition = yOffsetByCondition + delta;
            yLimsByCondition = yLimsByCondition + delta;
           
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
                            if isempty(td.listByCondition{conditionIdx(iCond)}), continue; end
                            % use first trial to draw marks and intervals
                            % on
                            td.alignInfoSet{idxAlign}.drawOnRasterByTrial('startByTrial', startData{iAlign, iCond}(1), ...
                                'stopByTrial', stopData{iAlign, iCond}(1), ...
                                'trialIdx', td.listByCondition{conditionIdx(iCond)}(1), ...
                                'showInLegend', iCond == 1, 'tOffsetZero', tOffsetByAlign(iAlign), ...
                                'yOffsetTop', yOffsetByCondition(iCond) + p.Results.annotationHeight + 1, ...
                                'tickHeight', p.Results.annotationHeight, ...
                                'intervalMinWidth', p.Results.intervalMinWidth, ...
                                'axh', axh, 'intervalAlpha', p.Results.intervalAlpha);
                        end
                    end
                else
                    error('Not yet implemented');
                end
            else
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    for iCond = 1:nConditionsUsed
                        if isempty(td.listByCondition{conditionIdx(iCond)}), continue; end
                        td.alignInfoSet{idxAlign}.drawOnRasterByTrial('startByTrial', startData{iAlign, iCond}, ...
                            'stopByTrial', stopData{iAlign, iCond}, ...
                            'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                            'showInLegend', iCond == 1, 'tOffsetZero', tOffsetByAlign(iAlign), ...
                            'yOffsetTop', yOffsetByCondition(iCond), ...
                            'intervalMinWidth', p.Results.intervalMinWidth, ...
                            'axh', axh, 'intervalAlpha', p.Results.intervalAlpha);
                    end
                end
            end
            %set(gcf, 'Renderer', 'painters');
            
            % draw tick rasters in a grid pattern (conditions vertically,
            % aligns horizontally)
            for iAlign = 1:nAlignUsed
                for iC = 1:nConditionsUsed
                    app = td.conditionAppearances(conditionIdx(iC));
                    if p.Results.colorSpikes
                        color = app.Color;
                    else
                        color = 'k';
                    end
                        
                    TrialDataUtilities.Plotting.drawTickRaster(timesByAlign{iAlign, iC}, ...
                        'xOffset', tOffsetByAlign(iAlign), 'yOffset', yOffsetByCondition(iC), ...
                        'color', color);
                    hold(axh, 'on');
                end
            end
            
            % setup y axis condition labels
            colors = cat(1, td.conditionAppearances(conditionIdx).Color);
            conditionNames = td.conditionNamesMultiline(conditionIdx);
            
            % only include conditions with at least 1 trial
            mask = trialCounts(conditionIdx) > 0;
            au = AutoAxis(axh);
            au.addLabeledSpan('y', 'span', yLimsByCondition(:, mask), 'label', ...
                conditionNames(mask), 'color', colors(mask, :));
            
            % setup time axis markers
            if p.Results.annotateAboveEachCondition
                if iAlign == nAlignUsed
                    % just use horizontal scale bar
                    au = AutoAxis(axh);
                    au.xUnits = td.timeUnitName;
                    au.addAutoScaleBarX();
                    au.update();
                    au.installCallbacks();
                end
            else
                % use marks or tickBridges via AlignSummary
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    td.alignSummarySet{idxAlign}.setupTimeAutoAxis('which', 'x', ...
                        'style', p.Results.timeAxisStyle, 'tOffsetZero', tOffsetByAlign(iAlign));
                end
            end
            
            TrialDataUtilities.Plotting.setTitleIfBlank(axh, '%s Unit %s', td.datasetName, unitName);
            
            axis(axh, 'tight');
            au = AutoAxis(axh);
            au.axisMarginLeft = p.Results.axisMarginLeft; % make room for left hand side labels
            axis(axh, 'off');
            au.update();
            
            hold(axh, 'off');
        end
    end

    % Plotting Analog each trial
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
            p.addParameter('alignIdx', 1:td.nAlign, @isvector);
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
                elseif isempty(tvecCell{iAlign})
                    error('Time cell for alignment %d is empty', iAlign);
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
                currentOffset = currentOffset + maxs(iAlign-1) + ...
                    gaps(iAlign-1) - mins(iAlign);
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
            p.addParameter('time', [], @(x) isvector(x) || iscell(x)); % for D == 1,2,3 (for marking)
            p.addParameter('data', {}, @iscell); % for D == 1,2,3
            p.addParameter('yOffsetBetweenTrials', 0, @isscalar);
            p.addParameter('yOffsetBetweenConditions', [], @isscalar);
            
            p.addParameter('axisInfoX', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoZ', [], @(x) isempty(x) || ischar(x) || isa(x, 'ChannelDescriptor'));
            
            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @ischar);
            p.addParameter('axisStyleY', 'tickBridge', @ischar);
            
            p.addParameter('conditionIdx', 1:td.nConditions, @isnumeric);
            p.addParameter('alignIdx', [], @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 5, @isscalar);
            p.addParameter('timeAxisStyle', 'marker', @ischar);
            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('useTranslucentMark3d', false, @islogical);
            p.KeepUnmatched;
            p.parse(varargin{:});
            
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
                    zlevel = iCond / (2*nConditionsUsed);
                    
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
                            dataC = bsxfun(@plus, dataC, yOffsets);
                            if p.Results.alpha < 1
                               hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(timeC' + tOffset, dataC', ...
                                   ones(size(timeC)) * zlevel, ...
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                            else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, timeC', dataC', '-', ...
                                    plotArgs{:}, p.Results.plotOptions{:});
                            end

                        elseif D == 2
                            if p.Results.alpha < 1
                               hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dataC(:, :, 1)', dataC(:, :, 2)', ...
                                   ones(size(dataC(:, :, 1)')) * zlevel, ...
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                            else
                                plotArgs = app(iCond).getPlotArgs();
                                hData{iCond, iAlign} = plot(axh, dataC(:, :, 1)', dataC(:, :, 2)', '-', ...
                                    plotArgs{:}, p.Results.plotOptions{:});
                            end

                        elseif D == 3
                            if p.Results.alpha < 1
                                hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dataC(:, :, 1)', dataC(:, :, 2)', dataC(:, :, 3)', ... 
                                   'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                   'LineWidth', app(iCond).LineWidth, 'Parent', axh, p.Results.plotOptions{:});
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
                                yOffset = yOffsets(iTrial);
                                
                                if ischar(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                                    xvec = dvec;
                                    yvec = tvec;
                                else
                                    yvec = dvec;
                                    xvec = tvec;
                                end                         

                                if ~isempty(tvec) && ~isempty(dvec)
                                    if p.Results.alpha < 1
                                       hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(xvec + tOffset, yvec + yOffset, ...
                                           ones(size(yvec)) * zlevel, ...
                                           'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                                           'LineWidth', app(iCond).LineWidth, p.Results.plotOptions{:});
                                    else
                                        plotArgs = app(iCond).getPlotArgs();
                                        hData{iCond, iAlign}(iTrial) = plot(axh, xvec + tOffset, yvec + yOffset, '-', ...
                                            plotArgs{:}, p.Results.plotOptions{:});
                                    end
                                end

                            elseif D==2
                                dxvec = dataC{iTrial}(:, 1);
                                dyvec = dataC{iTrial}(:, 2);

                                if ~isempty(dxvec) && ~isempty(dyvec)
                                    if p.Results.alpha < 1
                                       hData{iCond, iAlign}(iTrial) = TrialDataUtilities.Plotting.patchline(dxvec, dyvec, ...
                                           ones(size(dxvec)) * zlevel, ...
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
                    
                    td.alignInfoSet{idxAlign}.drawOnDataByTrial('time', timeC, 'data', dataC, ...
                        'trialIdx', td.listByCondition{conditionIdx(iCond)}, ...
                        'showInLegend', iCond == 1, 'tOffsetZero', timeOffsetByAlign(iAlign), ...
                        'axh', axh, 'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize);
                end

                % setup time axis for this align
                if D == 1
                    td.alignSummarySet{idxAlign}.setupTimeAutoAxis('which', 'x', 'style', p.Results.timeAxisStyle, 'tOffsetZero', timeOffsetByAlign(iAlign));
                end
            end
            
            % setup non-time axes
            if D == 1
                if ischar(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
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

            box(axh, 'off');
            axis(axh, 'tight');
            
            au = AutoAxis(axh);
            if strcmp(axisStyleY, 'scaleBar')
                au.axisMarginLeft = 0.8;
            end
%             if strcmp(axisStyleX, 'scaleBar')
%                 au.axisMarginBottom = 1;
%             end
            au.update();
            au.installCallbacks();
            hold(axh, 'off');
        end
        
        function plotAnalogGroupedEachTrial(td, name, varargin) 
            % grab raw data (for marking) and grouped data (for plotting)
            [dataByGroup, timeByGroup] = td.getAnalogGroupedEachAlign(name);     
            td.plotProvidedAnalogDataGroupedEachTrial(1, 'time', timeByGroup, ...
                'data', dataByGroup, 'axisInfoY', td.channelDescriptorsByName.(name), ...
                varargin{:});
        end
        
        function plotAnalogGroupedEachTrial2D(td, name1, name2, varargin) 
            [dataCell, tvec] = td.getMultiAnalogAsMatrixGrouped({name1, name2});
            td.plotProvidedAnalogDataGroupedEachTrial(2, ...
                'time', tvec, 'data', dataCell(:), ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), varargin{:});
        end
        
        function plotAnalogGroupedEachTrial2DvsTime(td, name1, name2, varargin)
            [dataCell, tvec] = td.getMultiAnalogAsMatrixGrouped({name1, name2});
            td.plotProvidedAnalogDataGroupedEachTrial(2, ...
                'time', tvec, 'data', dataCell(:), ...
                'axisInfoZ', 'time', ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), varargin{:});
        end
        
        function plotAnalogGroupedEachTrial3D(td, name1, name2, name3, varargin) 
            [dataCell, tvec] = td.getMultiAnalogAsMatrixGrouped({name1, name2, name3});
            td.plotProvidedAnalogDataGroupedEachTrial(3, ...
                'time', tvec, 'data', dataCell(:), ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), ...
                'axisInfoZ', td.channelDescriptorsByName.(name3), varargin{:});
        end
    end
    
    % Plotting Analog means by group
    methods
        function plotProvidedAnalogDataGroupMeans(td, D, varargin)
            p = inputParser();
            p.addParameter('time', [], @(x) true);
            p.addParameter('alignIdx', 1:td.nAlign, @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            alignTimeOffsets = td.getAlignPlottingTimeOffsets(p.Results.time, 'alignIdx', p.Results.alignIdx);
            TrialDataConditionAlign.plotConditionAlignedAnalogDataGroupMeans(D, varargin{:}, ...
                'conditionDescriptor', td.conditionInfo, ...
                'alignSummarySet', td.alignSummarySet, ...
                'alignInfoActiveIdx', td.alignInfoActiveIdx, ...
                'alignTimeOffsets', alignTimeOffsets);
        end
        
        function plotAnalogGroupMeans(td, name, varargin) 
            % plot the mean and sem for an analog channel vs. time within
            % each condition
            
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('errorType', 'sem', @(s) ismember(s, {'sem', 'std', ''}));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getAnalogGroupMeans(name, 'minTrials', p.Results.minTrials);     
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
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
            
            td.plotProvidedAnalogDataGroupMeans(1, 'time', tvecCell, ...
                'data', meanMat, 'dataError', errorMat, p.Unmatched, ...
                'axisInfoX', 'time', 'axisInfoY', td.channelDescriptorsByName.(name), p.Unmatched);
        end
        
        function plotAnalogGroupMeans2D(td, name1, name2, varargin) 
            % plot the mean and sem for an analog channel vs. time within
            % each condition
            
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('errorType', 'sem', @(s) ismember(s, {'sem', 'std', ''}));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getMultiAnalogGroupMeans({name1, name2}, 'minTrials', p.Results.minTrials);     
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, ...
                    semMat{iAlign}, stdMat{iAlign});
            end
            td.plotProvidedAnalogDataGroupMeans(2, 'time', tvecCell, ...
                'data', meanMat, p.Unmatched, ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), p.Unmatched);
        end
        
        function plotAnalogGroupMeans3D(td, name1, name2, name3, varargin) 
            % plot the mean and sem for an analog channel vs. time within
            % each condition
            
            p = inputParser();
            p.addParameter('minTrials', 1, @isscalar);
            p.addParameter('errorType', 'sem', @(s) ismember(s, {'sem', 'std', ''}));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % loop over alignments and gather mean data
            % and slice each in time to capture only the non-nan region
            [meanMat, semMat, tvecCell, stdMat] = deal(cell(td.nAlign, 1));
            for iAlign = 1:td.nAlign
                [meanMat{iAlign}, semMat{iAlign}, tvecCell{iAlign}, stdMat{iAlign}] = ...
                    td.useAlign(iAlign).getMultiAnalogGroupMeans({name1, name2, name3}, 'minTrials', p.Results.minTrials);     
                [tvecCell{iAlign}, meanMat{iAlign}, semMat{iAlign}, stdMat{iAlign}] = ...
                    TrialDataUtilities.Data.sliceValidNonNaNTimeRegion(tvecCell{iAlign}, meanMat{iAlign}, ...
                    semMat{iAlign}, stdMat{iAlign});
            end
            
            td.plotProvidedAnalogDataGroupMeans(3, 'time', tvecCell, ...
                'data', meanMat, p.Unmatched, ...
                'axisInfoX', td.channelDescriptorsByName.(name1), ...
                'axisInfoY', td.channelDescriptorsByName.(name2), ...
                'axisInfoZ', td.channelDescriptorsByName.(name3));
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
    
    methods(Static) % Utility drawing methods
        function plotConditionAlignedAnalogDataGroupMeans(D, varargin)
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
            %  data, dataErrorY is either
            %     nConditions x T x D matrix 
            %     nAlign x 1 cell of nConditions x T x D matrices
            %
            p = inputParser();
            p.addParameter('time', [], @(x) isvector(x) || iscell(x)); % for D == 1,2,3 (for marking)
            p.addParameter('data', {}, @(x) isnumeric(x) || iscell(x)); % for D == 1,2,3
            p.addParameter('dataError', {}, @(x) isnumeric(x) || iscell(x)); 
            
            p.addParameter('conditionDescriptor', [], @(x) isa(x, 'ConditionDescriptor'));
            p.addParameter('alignSummarySet', [], @iscell);
            p.addParameter('alignInfoActiveIdx', 1, @isscalar);
            p.addParameter('alignTimeOffsets', [], @isvector);
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('zOffset', 0, @isscalar);
            
            p.addParameter('axisInfoX', [], @(x) isempty(x) || ischar(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoY', [], @(x) isempty(x) || ischar(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));
            p.addParameter('axisInfoZ', [], @(x) isempty(x) || ischar(x) || isstruct(x) || isa(x, 'ChannelDescriptor'));
            
            p.addParameter('scaleBars', false, @islogical);
            p.addParameter('axisStyleX', 'tickBridge', @ischar);
            p.addParameter('axisStyleY', 'tickBridge', @ischar);
            
            p.addParameter('lineWidth', 0.5, @isscalar);
            p.addParameter('conditionIdx', [], @isnumeric);
            p.addParameter('alignIdx', [], @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            
            p.addParameter('alpha', 1, @isscalar); % alpha for main traces
            p.addParameter('errorAlpha', 0.5, @isscalar); % alpha for surrounding error fills
            
            p.addParameter('markShowOnData', true, @islogical);
            p.addParameter('markShowOnAxis', true, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 2, @isscalar);
            
            p.addParameter('intervalShowOnData', true, @islogical);
            p.addParameter('intervalShowOnAxis', true, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);
            
            p.addParameter('showRangesOnData', true, @islogical); % show ranges for marks on traces
            p.addParameter('showRangesOnAxis', true, @islogical); % show ranges for marks below axis
            
            p.addParameter('timeAxisStyle', 'marker', @ischar);
            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('useTranslucentMark3d', true, @islogical);
            
            p.addParameter('axh', [], @(x) true); % pass thru to getRequestedPlotAxis
            
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            zOffset = p.Results.zOffset;
            
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
                assert(~iscell(time) && ~iscell(data) && (isEmpty(dataError) || ~iscell(dataError)), ...
                    'Arguments must be nAlign x 1 cell or C x T x D matrices');
                nAlignUsed = 1;
                
                time = {time};
                data = {data};
                dataError = {dataError};
            end
              
            % check size of time/data cell contents and mask conditions if
            % needed
            for iA = 1:nAlignUsed
                if size(time{iA}, 1) == nConditionsUsed
                    % okay as is
                elseif size(time{iA}, 1) == cd.nConditions
                    % needs to be masekd
                    time{iA} = time{iA}(conditionIdx, :);
                end
                if size(data{iA}, 1) == nConditionsUsed
                    % okay as is
                elseif size(data{iA}, 1) == cd.nConditions
                    % needs to be masekd
                    data{iA} = data{iA}(conditionIdx, :);
                end
                if ~isempty(dataError)
                    if size(dataError{iA}, 1) == nConditionsUsed
                        % okay as is
                    elseif size(dataError{iA}, 1) == cd.nConditions
                        % needs to be masekd
                        dataError{iA} = dataError{iA}(conditionIdx, :);
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
            assert(numel(alignIdx) == nAlignUsed, 'numel(time) must match numel(alignIdx)');
            
            % compute time offsets between successive alignments
            if D == 1
                timeOffsetByAlign = p.Results.alignTimeOffsets;
            else
                timeOffsetByAlign = zeros(nAlignUsed, 1);
            end

            % store handles as we go
            hData = cell(nConditionsUsed, nAlignUsed);
            
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
                        if plotErrorY
                            hShade = errorshade(tvec + tOffset + xOffset, dmat + yOffset, ...                   
                                errmat, app(iCond).Color, 'axh', axh, ...
                                'alpha', p.Results.errorAlpha, 'z', 1, 'showLine', false);
                            TrialDataUtilities.Plotting.hideInLegend(hShade);
                        end
                        if p.Results.alpha < 1
                            hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(tvec + tOffset + xOffset, dmat + yOffset, ...
                               'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                               'LineWidth', p.Results.lineWidth, 'Parent', axh, p.Results.plotOptions{:});
                        else
                            plotArgs = app(iCond).getPlotArgs();
                            hData{iCond, iAlign} = plot(axh, tvec + tOffset + xOffset, dmat + yOffset, '-', ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                plotArgs{:}, p.Results.plotOptions{:});
                        end

                    elseif D == 2
                        if p.Results.alpha < 1
                           hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, ...
                               'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                               'LineWidth', p.Results.lineWidth, 'Parent', axh, p.Results.plotOptions{:});
                        else
                            plotArgs = app(iCond).getPlotArgs();
                            hData{iCond, iAlign} = plot(axh, dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, '-', ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ... 
                                plotArgs{:}, p.Results.plotOptions{:});
                        end

                    elseif D == 3
                        if p.Results.alpha < 1
                            hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, dmat(:, 3) + zOffset, ... 
                               'EdgeColor', app(iCond).Color, 'EdgeAlpha', p.Results.alpha, ...
                               'LineWidth', p.Results.lineWidth,  'Parent', axh, p.Results.plotOptions{:});
                        else
                            plotArgs = app(iCond).getPlotArgs();
                            hData{iCond, iAlign} = plot3(axh, dmat(:, 1) + xOffset, dmat(:, 2) + yOffset, dmat(:, 3) + zOffset, '-', ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ... 
                                plotArgs{:}, p.Results.plotOptions{:});
                        end
                    end  
                    
                    % update name for inclusion in legend
                    if iAlign == 1
                        TrialDataUtilities.Plotting.showFirstInLegend(hData{iCond, iAlign}, cd.names{iCond});
                    end
                end
                
                box(axh, 'off');
                axis(axh, 'tight');
            end
            
            % draw marks and intervals on data
            if ~isempty(alignSummarySet)
                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    alignSummarySet{idxAlign}.drawOnDataByCondition(time{iAlign}, ...
                        permute(data{iAlign}, [2 3 1]), ...  % data needs to be T x D x C x N, currently C x T x D
                        'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
                        'xOffset', xOffset, 'yOffset', yOffset, 'zOffset', zOffset, ...
                        'tOffsetZero', timeOffsetByAlign(iAlign), 'alpha', p.Results.alpha, ...
                        'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize, ...
                        'intervalAlpha', p.Results.intervalAlpha, ...
                        'showRanges', p.Results.showRangesOnData, ...
                        'tMin', min(time{iAlign}), 'tMax', max(time{iAlign}));
                end
            end
            
            % setup time axes for each alignment
            if ~isempty(alignSummarySet)
                if D == 1
                    for iAlign = 1:nAlignUsed
                        idxAlign = alignIdx(iAlign);

                        % setup x axis 
                        alignSummarySet{idxAlign}.setupTimeAutoAxis('axh', axh, 'tOffsetZero', timeOffsetByAlign(iAlign) + xOffset, ...
                            'tMin', min(time{iAlign}), 'tMax', max(time{iAlign}), ...
                            'style', p.Results.timeAxisStyle, 'showRanges', p.Results.showRangesOnAxis);
                    end
                end
            end
                
            % setup non-time axes
            if D == 1
                if ischar(p.Results.axisInfoY) && strcmp(p.Results.axisInfoY, 'time')
                    % x is data, y is time
                    if ~isempty(p.Results.axisInfoX)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                            'which', 'x', 'style', axisStyleX);
                    end
                else
                    % y is data, x is time
                    if ~isempty(p.Results.axisInfoY)
                        TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                            'which', 'y', 'style', axisStyleY);
                    end
                end

            elseif D == 2
                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                        'which', 'x', 'style', axisStyleX);
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                        'which', 'y', 'style', axisStyleY);
                end

            elseif D == 3
                if ~isempty(p.Results.axisInfoX)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoX, ...
                        'which', 'x', 'useAutoAxis', false);
                end
                if ~isempty(p.Results.axisInfoY)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoY, ...
                        'which', 'y', 'useAutoAxis', false);
                end
                if ~isempty(p.Results.axisInfoZ)
                    TrialDataUtilities.Plotting.setupAxisForChannel(p.Results.axisInfoZ, ...
                        'which', 'z', 'useAutoAxis', false);
                end

                if p.Results.useThreeVector
                    ThreeVector(axh);
                    axis(axh, 'off');
                    axis(axh, 'vis3d');
                end                   
            end

            box(axh, 'off');
            axis(axh, 'off');
            axis(axh, 'tight');
            
            au = AutoAxis(axh);
            au.update();
            au.installCallbacks();
            hold(axh, 'off');
            set(axh, 'SortMethod', 'childorder');
        end
    end
        
        
end
