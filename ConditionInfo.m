% ConditionInfo may be optionally bound to a specific set of trial data
% which will be passed to getAttributeValueFn() when requesting the values
% of each attribute for each trial. If not bound to anything, this function will receive
% [] as its first argument, with the assumption being that the function handle in this
% case is already bound to a specific set of trial data.
classdef ConditionInfo < handle & matlab.mixin.Copyable & ConditionDescriptor

    properties(Transient)
        % trialData from which all attribute values are derived, not
        % altered in any way to prevent copying
        trialData
    end
    
    properties
        % function with signature:
        % valuesByAttribute = getAttributeValueFn(trialData, attributeNames)
        % 
        % - trialData: typically a struct array or TrialData instance
        % - attributeNames: cellstr of attribute names (from the
        %     "requestAs" list)
        % - valuesByAttribute : struct array where v(iTrial).attributeName = attribute value on this trial
        % 
        % 
        getAttributeValueFn = @ConditionInfo.defaultGetAttributeFn;
        
        usingManualAttributeValues = false; 
        
        bound = false; 
    end
    
    properties(SetAccess=protected)
        % T is number of trials
        % A is number of attributes
        values = {}; % T x A cell array : value of attribute a on trial t

        % a mask over trials (T x 1). A trial is valid if all of its attribute values are in the
        % value lists for those attributes, AND manualInvalid(i) == 0
        manualInvalid % used by markInvalid
        
        % which condition does each trial belong to
        conditionIdx % T x 1 array of linear index into conditions for each trials
        
        % T x A matrix of which condition each trial belongs to as a row vector of subscript indices
        conditionSubs 
        
        % nConditions x 1 cell array of idx in each condition
        listByCondition
    end

    properties(Dependent)
        nTrials
        
        countByCondition
        
        nConditionsNonEmpty
        
        % a mask over trials (T x 1). A trial is valid if all of its attribute values are in the
        % value lists for those attributes, AND manualInvalid(i) == 0
        % logical mask indicating which trials are valid to include when returning groups
        % this mask does not affect any other functions for grabbing attribute values / unique attributes, etc.
        valid  

        nValid
    end

    methods 
        function ci = ConditionInfo(varargin)
            p = inputParser;
            p.addOptional('trialData', [], @(x) true);
            p.parse(varargin{:});

            if ~isempty(p.Results.trialData)
                ci.bindTrialData(p.Results.trialData);
            end
        end
        
        function bind(ci, td)
            ci.bindTrialData(td);
        end

        function bindTrialData(ci, td)
            ci.bound = true;
            ci.trialData = td;
            ci.updateFromTrialData();
        end

        % allows 'attributes' paramValue to be specified in which only specific
        % attributes have been changed
        function updateTrialData(ci, td, varargin)
            ci.trialData = td;
            ci.updateFromTrialData(varargin);
        end
    end
    
    methods % Superclass overrides
        function ci = updateCache(ci)
            % flush the current conditions / appearances / names buffers and list
            % of trials in each condition. These will be refreshed upon request.
            
            % handle this on demand instead of precomputing
            ci.conditions = [];
            ci.appearances = {};
            ci.names = {};
            ci.conditionIdx = [];
            ci.conditionSubs = [];
            ci.listByCondition = [];
            
            % the above accomplishes the same effect as in ConditionDescriptor
            % except that ConditionDescriptor recomputes everything then since it 
            % is a value class
            %ci = updateCache@ConditionDescriptor(ci); 
        end
        
        function ci = freezeAppearances(ci)
            % freeze current appearance information, but only store
            % conditions that have a trial in them now (which saves a bunch
            % of searching time)
            mask = ci.countByCondition > 0;
            ci.frozenAppearanceConditions = ci.conditions(mask);
            ci.frozenAppearanceData = ci.appearances(mask);
            ci.appearanceFn = @ConditionDescriptor.frozenAppearanceFn;
        end
    end
    
    methods(Access=protected)
        function ci = maskAttributes(ci, mask)
            ci = maskAttributes@ConditionDescriptor(ci, mask);
            ci.values = ci.values(:, mask);
        end 
    end

    methods
        function set.usingManualAttributeValues(ci, value)
            if value
                % dump trialData
                ci.trialData = [];
            end
        end

        function conditionIdx = get.conditionIdx(ci)
            if isempty(ci.conditionIdx)
                if ci.nTrials > 0
                    ci.conditionIdx = TensorUtils.subMat2Ind(ci.conditionsSize, ci.conditionSubs);
                else
                    conditionIdx = [];
                    return;
                end
            end
            conditionIdx = ci.conditionIdx;
        end
        
        function initializeWithNTrials(ci, N)
            % build empty arrays for N trials
            ci.manualInvalid = false(N, 1);
            ci.values = cell(N, ci.nAttributes);
            ci.conditionIdx = nan(N, 1);
        end

        function filterValidTrials(ci)
            % drop all invalid trials
            ci.resampleTrials(ci.valid);
        end

        function resampleTrials(ci, selector)
            assert(isvector(selector), 'Selector must be vector of indices or vector mask');
            ci.values = ci.values(selector, :);
            ci.manualInvalid = false(nnz(selector), 1);
            ci.conditionIdx = ci.conditionIdx(selector);
            ci.conditionSubs = ci.conditionSubs(selector, :);
            % auto-updates on request:
            ci.listByCondition = [];
        end
        
        function updateFromTrialData(ci, varargin)
            % entirely refresh the internal attribute value cache (and number of trials)
            % from .trialData. Optionally refresh only specific attributes 
            p = inputParser;
            p.addOptional('attributes', {}, @iscellstr);
            p.parse(varargin{:});
            
            attributes = p.Results.attributes;
            if isempty(attributes)
                partialUpdate = false;
                attributes = ci.attributeNames;
            else
                partialUpdate = true;
            end
            
            if ci.usingManualAttributeValues
                error('This ConditionInfo has .usingManualAttributeValues == true');
            end
            
            if isempty(ci.trialData)
                warning('Cannot update from trialData because .trialData is empty');
                return;
            end
            
            % update trialCount to match length(trialData)
            nTrials = numel(ci.trialData);
            if ~partialUpdate
                ci.initializeWithNTrials(nTrials);
            else
                assert(nTrials == ci.nTrials, 'Number of trials must match');
            end

            if ci.nAttributes > 0 && ci.nTrials > 0
                attrInds = find(strcmp(attributes, ci.attributeNames));
                
                % fetch valuesByAttribute using callback function
                valueStruct = ci.requestAttributeValues(attributes);
                valueCell = struct2cell(valueStruct)';
                
                % store in .values cell
                ci.values(:, attrInds) = valueCell;
            end
            
            ci.updateCache();
        end
        
        function addAttribute(ci, name, varargin)
            % calls ConditionDescriptor.addAttribute but also accepts a vector of
            % attribute values as paramValue 'values', values, ...
            %
            % The unique set of values in this list, excluding
            % NaN and '', will constitute the valueList for that attribute, unless
            % 'valueList', valueList is manually specified as a paramValue
            p = inputParser;
            p.addRequired('name', @(x) ischar(x) || iscell(x));
            p.addParamValue('requestAs', '', @(x) ischar(x) || iscell(x)); % if specified, the attribute value will be requested as this instead of `name`
            p.addParamValue('values', {}, @(x) true);
            p.addParamValue('valueList', {}, @(x) true); % cell of values for this attribute, or cell of cells for multiple attributes
            p.addParamValue('matrix', [], @(x) islogical(x));
            p.KeepUnmatched = true;
            p.parse(name, varargin{:});
            values = p.Results.values;
            valueList = p.Results.valueList;
            requestAs = p.Results.requestAs;
            matrix = p.Results.matrix;
            if isempty(requestAs)
                requestAs = name;
            end

            if ~iscell(name)
                name = {name};
                % make valueList a cell of cells
                if ~iscell(valueList)
                    valueList = num2cell(valueList);
                end
                if ~iscell(valueList)
                    valueList = {valueList};
                end
            end
            if ~iscell(values) 
                values = num2cell(values);
            end
            if ~iscell(requestAs)
                requestAs = {requestAs};
            end
            
            nAttrAdd = numel(name);
            assert(numel(name) == numel(requestAs), 'List of name and requestAs must be same length');
            
            % check whether values was provided iff .usingManualAttributeValues
            if isempty(values)
                if ci.usingManualAttributeValues
                    error('.usingManualAttributeValues == true, please provide ''values'' param with the values for each attribute');
                end
                nTrialsAdd = numel(ci.trialData);
            else
                if ~ci.usingManualAttributeValues
                    error('.usingManualAttributeValues == false, ''values'' param should not be provided');
                end
                nTrialsAdd = size(values, 1);
                % check size of provided values 
                if ~isempty(ci.values)
                    % if the number of trials is already determined
                    assert(nTrialAdd == ci.nTrials, 'Values must be nTrials x number of new attributes cell array');
                end
                assert(size(values, 2) ==  nAttrAdd, 'Values must be nTrials x number of new attributes cell array');
            end
            
            % request the values first before making any modifications
            % in case something goes wrong
            if ~ci.usingManualAttributeValues
                valueStruct = ci.requestAttributeValues(name, requestAs);
                valueStruct = orderfields(valueStruct, requestAs);
                values = struct2cell(valueStruct);
            end

            % initialize/expand the ci.values array to the correct size
            nAttrOld = ci.nAttributes;
            if isempty(ci.values)
                % initialize the ci.values array
                % nTrialsAdd may be 0
                ci.initializeWithNTrials(nTrialsAdd);
                ci.values = cell(ci.nTrials, nAttrAdd);
            else
                % expand the existing ci.values array
                ci.values(:, nAttrOld+1:nAttrOld+nAttrAdd) = cell(nTrialsAdd, nAttrAdd);
            end
            
            ci.values(:, nAttrOld+1:nAttrOld+nAttrAdd) = values;
            
            % auto compute unique value lists?
            if isempty(valueList)
                valueList = cellvec(nAttrAdd);
                if ci.nTrials > 0
                    for iA = 1:nAttrAdd
                        valuesThis = values(iA, :);
                        % use different behavior depending on whether all values are scalars
                        [tf mat] = isScalarCell(valuesThis);
                        if tf 
                            valueUnique = unique(removenan(mat));
                        else
                            valueUnique = setdiff(unique(valuesThis), {''});
                        end
                        valueList{iA} = makecol(valueUnique);
                    end
                end
            else
                assert(numel(valueList) == nAttrAdd, 'ValueList must be a cell of cells with length nAttrAdd');
            end
            
            % defer to ConditionDescriptor to add the attributes to the remaining lists
            for i = 1:numel(name)
                ci = addAttribute@ConditionDescriptor(ci, name{i}, 'valueList', valueList{i}, 'requestAs', requestAs{i}, p.Unmatched);
            end
        end
        
        function valueStruct = requestAttributeValues(ci, attrNames, requestAs)
            % lookup requestAs name if not specified
            if nargin < 3
                inds = find(strcmp(ci.attributeNames, attrNames));
                requestAs = ci.attributeRequestAs(inds);
            end
                
            % translate into request as names
            if numel(ci.trialData) == 0
                valueStruct = struct();
            end
            
            
            valueStructRequestAs = ci.getAttributeValueFn(ci.trialData, requestAs);
            
            % check the returned size and field names
            assert(numel(valueStructRequestAs) == numel(ci.trialData), 'Number of elements returned by getAttributeFn must match nTrials');
            assert(all(isfield(valueStructRequestAs, requestAs)), 'Number of elements returned by getAttributeFn must match nTrials');
            
            % translate back into attribute names
            valueStruct = mvfield(valueStructRequestAs, requestAs, attrNames);
            valueStruct = orderfields(valueStruct, attrNames);
            valueStruct = makecol(valueStruct);
        end

        function subsMat = get.conditionSubs(ci)
            %subsMat = TensorUtils.ind2subAsMat(ci.conditionsSize, ci.conditionIdx);
            if isempty(ci.conditionSubs)
                if ci.nConditions > 0 && ci.nTrials > 0
                    if ci.nAttributesGroupBy == 0
                        % not grouping on anything, only 1 condition
                        assert(ci.nConditions == 1);
                        
                        % it's in condition 1 if it matches all of the
                        % non-grouped attributes, NaN otherwise
                        % filter by all other attributes having values in the
                        % valueList
                        subsMat = ones(ci.nTrials, 1);
                        for iA = 1:ci.nAttributes
                            attr = ci.attributeNames{iA};
                            values = ci.values(:, iA);
                            invalid = ~ismemberCell(values, ci.attributeValueList{iA});
                            subsMat(invalid, :) = NaN;
                        end
                    else
                        % build a matrix of subscripts subs{iAttr}(iTrial) is the index into 
                        % valueList for attribute iAttr of trial iTrial's attribute iAttr value
                        subsMat = nan(ci.nTrials, ci.nAttributesGroupBy);
                        for i = 1:ci.nAttributesGroupBy
                            % translate groupByList idx into attribute idx
                            iAttr = ci.groupByListAttributeIdx(i); 
                            values = ci.values(:, iAttr);
                            [~, subsMat(:, i)] = ismemberCell(values, ci.attributeValueList{iAttr});
                        end

                        % filter by all other attributes having values in the
                        % valueList
                        for iA = 1:ci.nAttributes
                            attr = ci.attributeNames{iA};
                            if ~ismember(attr, ci.groupByList)
                                values = ci.values(:, iA);
                                invalid = ~ismember(values, ci.attributeValueList{iA});
                                subsMat(invalid, :) = NaN;
                            end
                        end
                    end
                else
                    subsMat = [];
                    return;
                end
                ci.conditionSubs = subsMat;
                
                ci.conditionSubs(any(ci.conditionSubs == 0, 2), :) = NaN;
                ci.conditionSubs(ci.manualInvalid, :) = NaN;
            end
            
            subsMat = ci.conditionSubs;
        end
        
        function list = get.listByCondition(ci)
            if isempty(ci.listByCondition)
                list = cell(ci.conditionsSize);
                for iC = 1:ci.nConditions
                    list{iC} = makecol(find(ci.conditionIdx == iC));
                    if isempty(list{iC})
                        % ensure it can be concatenated into a column
                        % vector using cell2mat
                        list{iC} = nan(0, 1);
                    end
                end
                ci.listByCondition = list;
            else
                list = ci.listByCondition;
            end
        end

        function counts = get.countByCondition(ci)
            counts = cellfun(@length, ci.listByCondition);
        end
        
        function nConditions = get.nConditionsNonEmpty(ci)
            nConditions = nnz(~cellfun(@isempty, ci.listByCondition));
        end

        function nt = get.nTrials(ci)
            nt = size(ci.values, 1);
        end

        function markInvalid(ci, invalid)
            ci.manualInvalid(invalid) = true;
        end

        function valid = get.valid(ci)
            % return a mask which takes into account having a valid value for each attribute
            % specified, as well as the markInvalid function which stores its results in .manualInvalid
            valid = ~ci.manualInvalid & ci.getIsTrialInSomeGroup();
        end

        function nValid = get.nValid(ci)
            nValid = nnz(ci.valid);
        end
        
        function mask = getIsTrialInSomeGroup(ci)
            mask = ~isnan(ci.conditionIdx);
        end
    end

    methods % Resampling and shuffling
        function listByCondition = getListByConditionResampled(ci, varargin)
            % randomly resamples with replacement from the trials within each condition. 
            % The total number of idx in each condition will remain
            % the same but will typically contain duplicates and omissions of the
            % original set of trial idx
            %
            % listByCondition is in a format identical to ci.listByCondition
            % i.e. a tensor cell of size .conditionsSize containing vectors
            % of trialIdx belonging to that condition in that resample.

            listByCondition = TensorUtils.map(@ci.sampleWithReplacement, ci.listByCondition);
        end

        function listByCondition = getListByConditionResampledFromSingleAttributeValue(ci, attr, value)
            % randomly resamples with replacement from the trials within each condition. 
            % The total number of idx in each condition will remain
            % the same but will typically contain duplicates and omissions of the
            % original set of trial idx
            %
            % listByCondition is in a format identical to ci.listByCondition
            % i.e. a tensor cell of size .conditionsSize containing vectors

            attrIdx = ci.getAttributeIdxInGroupByList(attr);
            valueIdx = ci.getAttributeValueIdx(attr, value);

            listByCondition = TensorUtils.mapSlicesInPlace(@resampleFn, attrIdx, ci.listByCondition);

            function listAlongResampled = resampleFn(listAlong)
                % listAlong is a cell array of cell arrays containing trial idx 
                % be sure to preserve the orientation of listAlong 
                sz = size(listAlong);
                nPerCell = makecol(num2cell(cellfun(@length, listAlong)));
                sampleFromList = listAlong{valueIdx};

                listAlongResampled = cellfun(@(n) ci.sampleWithReplacement(sampleFromList, n), ...
                    nPerCell, 'UniformOutput', false); 
                listAlongResampled = reshape(listAlongResampled, sz);
            end
        end

        function listByCondition = getListByConditionShuffledAlong(ci, compareAlong)
            attrIdx = ci.getAttributeIdxInGroupByList(compareAlong);
            listByCondition = TensorUtils.mapSlices(@shuffleFn, attrIdx, ci.listByCondition);

            function listAlongShuffled = shuffleFn(listAlong)
                % listAlong is a cell array of cell arrays containing trial idx 
                % be sure to preserve the orientation of listAlong 
                sz = size(listAlong);
                nPer = cellfun(@length, listAlong);
                combinedIdx = cell2mat(makecol(squeeze(listAlong)));
                scrambledIdx = combinedIdx(randperm(length(combinedIdx)));
                splitIdx = mat2cell(scrambledIdx, nPer, 1);
                listAlongShuffled = reshape(splitIdx, sz);
            end
        end

        % utility function for linking the functions above (which reorder trial
        % memberships) to a new ConditionInfo built according to those trial 
        % memberships)
        function [ciCopy idxListNew] = buildCopyFromListByCondition(ci, listByCondition)
            % we have to build up a new list of attributes and a new list of trial idx
            
            attrInGroupByMask = ci.isAttributeInGroupByList;
            groupByListAttributeIdx = ci.groupByListAttributeIdx; 
            
            nPerCondition = cellfun(@length, listByCondition);

            nTrialsNew = sum(nPerCondition(:));
            newInvalid = false(nTrialsNew, 1);
            newValues = cell(nTrialsNew, ci.nAttributes);
            newListByCondition = cell(ci.conditionsSize);
            newConditionIdx = nan(nTrialsNew, 1);
            idxListNew = nan(nTrialsNew, 1);
            
            iTrialStart = 1;
            for iC = 1:ci.nConditions
                % set the all idx in this cell of the shuffle to the corresponding attribute value
                idx = listByCondition{iC};
                trialList = iTrialStart:(iTrialStart + nPerCondition(iC) - 1);

                idxListNew(trialList) = idx;
                newConditionIdx(trialList) = iC;
                newListByCondition{iC} = trialList;

                % store original attr values for attributes that aren't being grouped 
                newValues(trialList, ~attrInGroupByMask) = ci.values(idx, ~attrInGroupByMask);

                % overwrite original attr values for attributes that ARE being grouped
                % so that these trials actually look like they belong in this condition
                condInfo = ci.conditions(iC);
                for iAttrGroupBy = 1:ci.nAttributesGroupBy
                    [newValues{trialList, groupByListAttributeIdx(iAttrGroupBy)}] = ...
                        deal(condInfo.(ci.groupByList{iAttrGroupBy}));
                end
                
                iTrialStart = iTrialStart + nPerCondition(iC);
            end

            % build new ciCopy with the new resampling
            ciCopy = ci.copy();
            ciCopy.values = newValues;
            ciCopy.manualInvalid = newInvalid;

            % be sure this is all you need to do when updating the values
            ciCopy.conditionIdx = newConditionIdx;
            ciCopy.listByCondition = newListByCondition;
            % auto updates
            ciCopy.conditionSubs = [];
        end
         
        % sample trials with
        % replacement within each condition. Also returns a selection idx list
        % newIdxList which is nTrials x 1 which you can use to resample any associated
        % data to align that data with the resampled conditionInfo attributes
        function [ciCopy idxListNew] = buildResampled(ci, varargin)
            listByCondition = ci.getListByConditionResampled();
            [ciCopy idxListNew] = ci.buildCopyFromListByCondition(listByCondition);
        end

        function [ciCopy idxListNew] = buildResampledFromSingleAttributeValue(ci, attr, value, varargin)
            listByCondition = ci.getListByConditionResampledFromSingleAttributeValue(attr, value);
            [ciCopy idxListNew] = ci.buildCopyFromListByCondition(listByCondition);
        end

        function [ciCopy idxListNew] = buildShuffledAlong(ci, compareAlong, varargin)
            listByCondition = ci.getListByConditionShuffledAlong(compareAlong);
            [ciCopy idxListNew] = ci.buildCopyFromListByCondition(listByCondition);
        end

        function idxResample = sampleWithReplacement(ci, idx, nSample)
            assert(isempty(idx) || isvector(idx), 'Idx must be a vector');
            if nargin == 2
                nSample = length(idx);
            end
       
            idxResample = idx(randi([1 length(idx)], [nSample 1]));
        end
    end
    
    methods % Grouping and comparison axis shortcuts (largely legacy functions that wrap dependent properties)
        % builds a multidimensional cell array of indices, each of which contains a list
        % of indices that match the particular combination of attribute values for each attribute
        % in the groupBy list. The name of that group (generated by concatenating the attribute values) is
        % found in the corresponding cell of nameCell. infoStructArray contains the value of each attribute
        % in the field by the same name in the corresponding position.
        %
        % The list of trial indices excludes trials with an empty or nan value for any of the attributes, 
        % and excludes any trials marked as invalid in ci.valid (see ci.markInvalid)
        %
        % For example, if groupBy is {'target', 'stim'} and there are 4 targets and 2 stim classes,
        % then idxCell will be 4 x 2, with idxCell{i,j} containing trials for target i and stim class j.
        function [idxCell names appearanceProperties attributeValuesByCondition] = getGroups(ci)
            idxCell = ci.listByCondition; 
            names = ci.names;    
            appearanceProperties = ci.appearances;
            attributeValuesByCondition = ci.conditions;
        end

        % this is an expanded version of compareSlices / compareAlong in ConditionDescriptor
        function [trialIdxCompare, nameCompare, appearanceCompare, ...
                conditionIdxCell, conditionDescriptorOuter, conditionDescriptorInnerCell, conditionDescriptorInnerCommon] = ...
                getComparisonAxis(ci, compareAlong, varargin)

            if ~exist('compareAlong', 'var')
                error('Usage: getComparisonAxis(compareAcross attribute name)');
            end

            [conditionIdxCell, conditionDescriptorOuter, conditionDescriptorInnerCell, conditionDescriptorInnerCommon] = ...
                ci.compareSlices(compareAlong, varargin{:});

            nameCompare = conditionDescriptorOuter.names;
            appearanceCompare = conditionDescriptorOuter.appearances;
            trialIdxCompare = cellfun(@(inner) inner.listByCondition, conditionDescriptorInnerCell, 'UniformOutput', false);
        end

        % randomly shuffles the attribute values along attribute compareAlong 
        % of trials with each label
        function idxCompareShuffles = getConditionIdxShuffledAlong(ci, compareAlong, varargin)
            [idxCompare] = ci.getComparisonAxis(compareAlong);

            % nCompare is the number of unique attribute value tuples for attributes other than compareAlong
            nCompare = numel(idxCompare);
            [idxListByComparison nIdxPerByComparison nIdxTotalByComparison] = deal(cell(nCompare, 1));

            for iCompare = 1:nCompare
                % build a list of all trial idx across the comparison attribute for easy permutation
                idxListByComparison{iCompare} = cat(1, idxCompare{iCompare}{:});

                % and a list of how many entries are in each category
                nIdxByComparison{iCompare} = cellfun(@length, idxCompare{iCompare});

                % total number of trials in each category
                nIdxTotalByComparison{iCompare} = sum(nIdxByComparison{iCompare});
            end

            idxCompareShuffles = cell(nShuffles,1); 
            for iShuffle = 1:nShuffles
                idxCompareShuffles{iShuffle} = cell(nCompare, 1);
                for iCompare = 1:nCompare
                    % reorder the concatenated list of trial idx
                    reorderedList = idxListByComparison{iCompare}(randperm(nIdxTotalByComparison{iCompare}));

                    if isempty(reorderedList)
                        continue;
                    end
                    % and re-split it into the appropriate number of elements per condition along the compare axis
                    idxCompareShuffles{iShuffle}{iCompare} = mat2cell(reorderedList, nIdxByComparison{iCompare}, 1);
                end
            end
        end

        % similar to getComparisonAxis, except randomly relabels each trial with any 
        % label along the comparison axis while preserving the original proprortion
        % of trials with each label
        function [idxCompareShuffles nameCompare appearanceCompare] = ...
                getComparisonAxisShuffled(ci, compareAcross, varargin)
            nShuffles = 100;
            assignargs(varargin);

            [idxCompare nameCompare appearanceCompare] = ci.getComparisonAxis(compareAcross);

            % nCompare is the number of unique attribute value tuples for attributes other than compareAcross
            nCompare = numel(idxCompare);
            [idxListByComparison nIdxPerByComparison nIdxTotalByComparison] = deal(cell(nCompare, 1));

            for iCompare = 1:nCompare
                % build a list of all trial idx across the comparison attribute for easy permutation
                idxListByComparison{iCompare} = cat(1, idxCompare{iCompare}{:});

                % and a list of how many entries are in each category
                nIdxByComparison{iCompare} = cellfun(@length, idxCompare{iCompare});

                % total number of trials in each category
                nIdxTotalByComparison{iCompare} = sum(nIdxByComparison{iCompare});
            end

            idxCompareShuffles = cell(nShuffles,1); 
            for iShuffle = 1:nShuffles
                idxCompareShuffles{iShuffle} = cell(nCompare, 1);
                for iCompare = 1:nCompare
                    % reorder the concatenated list of trial idx
                    reorderedList = idxListByComparison{iCompare}(randperm(nIdxTotalByComparison{iCompare}));

                    if isempty(reorderedList)
                        continue;
                    end
                    % and re-split it into the appropriate number of elements per condition along the compare axis
                    idxCompareShuffles{iShuffle}{iCompare} = mat2cell(reorderedList, nIdxByComparison{iCompare}, 1);
                end
            end
        end

        % return a new ConditionInfo instance where trials identities have been shuffled among 
        % different attribute values of attribute `compareAcross`, preserving the trial
        % counts within each bin, and preserving all other attribute values
        function ciCell = buildShuffledAlongComparisonAxis(ci, compareAcross, varargin)
            nShuffles = 1;
            assignargs(varargin);

            [idxCompareShuffles] = ci.getComparisonAxisShuffled(compareAcross, 'nShuffles', nShuffles);

            nShuffles = length(idxCompareShuffles);
            attrIdx = ci.getAttributeIdxInGroupByList(compareAcross);
            attrValues = ci.getAttributeValueList(compareAcross);

            ciCell = cell(nShuffles, 1);
            for iShuffle = 1:nShuffles
                ciCopy = ci.copy();
                idxCompare = idxCompareShuffles{iShuffle};
                for iAlongAxis = 1:length(idxCompare)
                    % within means over the values of the compareAcross attribute
                    % which we here want to shuffle
                    attrValue = attrValues{iAlongAxis};
                    for iWithinAxis = 1:length(idxCompare{iAlongAxis})
                        % set the all idx in this cell of the shuffle to the corresponding attribute value
                        idx = idxCompare{iAlongAxis}{iWithinAxis};
                        [ciCopy.values{idx, attrIdx}] = deal(attrValue);
                    end
                end

                ciCell{iShuffle} = ciCopy;
            end
        end

        % similar to getComparisonAxis, except randomly resamples with replacement
        % trials from each group. The total number of idx in each set will remain
        % the same but will typically contain duplicates and omissions of the
        % original set of idx
        % DRAWS ALL SAMPLES FROM THE FIRST CONDITION ALONG THE COMPARE AXIS
        % for every value of the compare attribute
        function [idxCompareResampled nameCompare appearanceCompare] = ...
                getComparisonAxisResampledFromSame(ci, compareAcross, varargin)
            if ~exist('compareAcross', 'var')
                error('Usage: getComparisonAxisResampled(compareAcross attribute name)');
            end

            nResample = 100;
            resampleFrom = 1; % all samples will be drawn from this
            assignargs(varargin);

            [idxCompare nameCompare appearanceCompare] = getComparisonAxis(ci, compareAcross);

            % nCompare is the number of unique attribute value tuples for attributes other than compareAcross
            nCompare = numel(idxCompare);
            [idxListByComparison nIdxPerByComparison nIdxTotalByComparison] = deal(cell(nCompare, 1));
            
            for iCompare = 1:nCompare
                % build a list of all trial idx across the comparison attribute for easy permutation
                idxListByComparison{iCompare} = cat(1, idxCompare{iCompare}{:});

                % and a list of how many entries are in each category
                nIdxByComparison{iCompare} = cellfun(@length, idxCompare{iCompare});

                % total number of trials in each category
                nIdxTotalByComparison{iCompare} = sum(nIdxByComparison{iCompare});
            end

            % nPerCompare is the number of unique attribute values for compareAcross
            nPerCompare = numel(idxCompare{1});

            idxCompareResampled = cell(nResample,1); 
            for iResample = 1:nResample
                idxCompareResample{iResample} = cell(nCompare, 1);

                % loop over the unique non-compareAcross-attribute value tuples
                for iCompare = 1:nCompare
                    idxCompareResample{iResample}{iCompare} = cell(nPerCompare, 1);

                    % loop over the unique values of attr compareAcross
                    for iWithinCompare = 1:nPerCompare
                        % resample with replacement the list of idx from the resampleFrom'th list
                        % but take the same amount of trials as originally in this bin
                        idxCompareResampled{iResample}{iCompare}{iWithinCompare} = ...
                            ci.sampleWithReplacement(idxCompare{iCompare}{resampleFrom}, ...
                            nIdxByComparison{iCompare}(iWithinCompare));
                    end
                end
            end
        end

    end

    methods % Convenience attribute accessor methods
        
        function values = getAttribute(ci, name)
            idx = ci.getAttributeIdx(name);
            values = ci.values(:, idx);
            if ci.attributeNumeric(idx)
                values = cell2mat(values);
            end
            
            % exclude any values not found in the valueList, if specified
            valueList = ci.getAttributeValueList(name);
            if ~isempty(valueList)
                invalidValues = ~ismember(values, valueList);
                if iscellstr(values)
                    [values{invalidValues}] = deal('');
                elseif iscell(values)
                    [values{invalidValues}] = deal(NaN);
                else
                    values(invalidValues) = NaN;
                end
            end
        end

        function values = getAttributeUnique(ci, name)
            % if a value list is specified, we simply return that
            % otherwise, return the unique list of values
            valueList = ci.getAttributeValueList(name);
            if ~isempty(valueList)
                values = valueList; 
            else
                values = unique(ci.getAttribute(name));

                % remove empty and nan values
                if isnumeric(values) || islogical(values)
                    remove = isnan(values) | isempty(values);
                    values = num2cell(values);
                else
                    remove = cellfun(@isempty, values);
                end
                values(remove) = [];
            end
        end

        function valueCell = getMultipleAttributeUnique(ci, names)
           valueCell = cellfun(@(name) ci.getAttributeUnique(name), names, ...
                'UniformOutput', false);
        end
    
        function valueCell = getAllUnique(ci)
            valueCell = ci.getMultipleAttributeUnique(ci.attributeNames);
        end

        function idx = getIdxWithAttributeValue(ci, name, value)
            values = ci.getAttribute(name);
            if isnumeric(values) || islogical(values)
                values = num2cell(values);
            end
            match = cellfun(@(x) isequal(x, value), values);
            idx = find(match);
        end

        function idxCell = getIdxEachAttributeValue(ci, name)
            values = ci.getAttributeUnique(name);
            nValues = length(values);
            idxCell = cell(nValues, 1);

            if isnumeric(values) || islogical(values)
                values = num2cell(values);
            end
            for iValue = 1:nValues
                idxCell{iValue} = ci.getIdxWithAttributeValue(name, values{iValue});
            end
        end
        
        function idxList = getAttributeAsIdxUnique(ci, name)
            % return a numeric vector of the attribute value for trial i
            % like getAttribute, but instead of the raw value return an
            % index into getAttributeUnique(name)
            
            idxList = nan(ci.nTrials, 1);
            idxCell = ci.getIdxEachAttributeValue(name);
            for i = 1:length(idxCell);
                idxList(idxCell{i}) = i;
            end
        end
    end

    methods % ConditionDescriptor builder
        % build a static ConditionDescriptor for the current groupByList
        function cd = getConditionDescriptor(ci, varargin)
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end
    end

    methods(Static) 
        % Building from a condition descriptor with an accessor method
        function ci = fromConditionDescriptor(cd, varargin)
            p = inputParser;
            p.addOptional('trialData', [], @(x) true);
            p.addParamValue('getAttributeFn', @ConditionInfo.defaultGetAttributeFn, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            
            % build up the condition info
            ci = ConditionInfo();
            % Have conditionDescriptor copy over the important details
            ci = ConditionDescriptor.fromConditionDescriptor(cd, ci);
            % bind the trialData
            if ~isempty(p.Results.trialData
                ci.bind(p.Results.trialData);
            end
        end

        % return a scalar struct with one field for each attribute containing the attribute values
        % as a cell array or numeric vector
        function values = defaultGetAttributeFn(data, attributeNames, varargin)
            assert(isstruct(data) || isa(data, 'TrialData'), 'Please provide getAttributeFn if data is not struct or TrialData');
            valuesByAttribute = struct();
            for iAttr = 1:length(attributeNames)
                attr = attributeNames{iAttr};
                
                if isstruct(data)
                    valuesByAttribute.(attr) = {data.(attr)};
                elseif isa(data, 'TrialData')
                    valuesByAttribute.(attr) = data.getParam(attr);
                end
            end
            values = structOfArraysToStructArray(valuesByAttribute);
        end
        
        % same as ConditionDescriptor, except skips conditions with no
        % trials so that the colors stay maximally separated
        function a = defaultAppearanceFn(ci, varargin)
            % returns a struct specifying the default set of appearance properties 
            % for the given group. indsGroup is a length(ci.groupByList) x 1 array
            % of the inds where this group is located in the high-d array, and dimsGroup
            % gives the full dimensions of the list of groups.
            %
            % We vary color along all axes simultaneously, using the linear
            % inds.
            %
            % Alternatively, if no arguments are passed, simply return a set of defaults

            conditionsSize = ci.conditionsSize;
            nConditions = ci.nConditions;
            nConditionsNonEmpty = ci.nConditionsNonEmpty;
            
            a = emptyStructArray(ci.conditionsSize, {'color', 'lineWidth'});

            if nConditionsNonEmpty == 1
                cmap = [0.3 0.3 1];
            else
                cmap = distinguishable_colors(nConditionsNonEmpty);
            end
             
            colorInd = 1;
            for iC = 1:nConditions
                 if ci.countByCondition(iC) == 0
                     a(iC).lineWidth = 1;
                     a(iC).color = 'k';
                 else
                     a(iC).lineWidth = 2;
                     a(iC).color = cmap(colorInd, :);
                     colorInd = colorInd + 1;
                 end
            end
        end

    end

end
