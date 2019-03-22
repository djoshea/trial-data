% ConditionInfo may be optionally bound to a specific set of trial data
% which will be passed to getAttributeValueFn() when requesting the values
% of each attribute for each trial. If not bound to anything, this function will receive
% [] as its first argument, with the assumption being that the function handle in this
% case is already bound to a specific set of trial data.
%
% NOTE: shuffling and resampling along axes affects listByCondition, but not conditionSubs
classdef ConditionInfo < ConditionDescriptor

    properties(Hidden)
        % function with signature:
        % valuesByAttribute = getAttributeValueFn(trialData, attributeNames)
        %
        % - trialData: typically a struct array or TrialData instance
        % - attributeNames: cellstr of attribute names
        % - valuesByAttribute : struct array where v(iTrial).attributeName = attribute value on this trial
        getAttributeValueFn = @ConditionInfo.defaultGetAttributeFn;

        % function with signature:
        % nTrials = getNTrialsFn(trialData)
        getNTrialsFn = @ConditionInfo.defaultGetNTrialsFn;
    end

    properties(SetAccess=protected)
        % has apply to trial data already been called?
        applied = false;

        % T is number of trials
        % A is number of attributes
        values % T x A cell array : values{iTrial, iAttr} = value

        % a mask over trials (T x 1). A trial is valid if all of its attribute values are in the
        % value lists for those attributes, AND manualInvalid(i) == 0
        manualInvalid % used by markInvalid

        % a special mode in which all condition membership lists are set manually to this value
        % set via setConditionMembershipManual to allow manual condition definitions and weighting
        % T x C matrix either logical or numeric
        conditionMembershipRawManual
         % T x 1 cellstr of manually stored reasons why a given trial is excluded from all conditions
         % beyond it not satisfying the attribute filters
        conditionMembershipRawManualReasonInvalid

        % alternative manual mode to conditionMembershipRawManual. function signature looks like:
        % [membership, reasonInvalid] = ci.conditionMembershipManualFn(ci, membership, reasonInvalid, matchesFilters, varargin);
        % varargin is to support additional arguments in the future
        conditionMembershipManualFn
    end

    %%% End of properties saved to disk

    % Properties which wrap eponymous properties inside odc (on-demand cache)
    properties(Dependent, Transient, SetAccess=protected)
        % which condition does each trial belong to
        conditionIdx % T x 1 array of linear index into conditions for each trials

        % T x A matrix of which condition each trial belongs to as a row vector of subscript indices
        conditionSubsRaw

        % nTrials x 1 cellstr of reasons why a trial is not valid for this condition descriptor, without
        % worrying about manual invalid
        invalidCause

        % T x A matrix of which condition each trial belongs to as a row
        % vector of subscript indices, except invalid trials will have all
        % NaNs in their row
        conditionSubs

        % T x conditionSize logical table indicating whether each trial
        % belongs in each condition.
        conditionMembership % reflects randomization
        conditionMembershipRaw % does not reflect randomization

        % nConditions x 1 cell array of idx in each condition
        listByConditionRaw
        % nConditions x 1 cell array of corresponding weights for each included trial in that condition
        listByConditionWeightsRaw

        % nConditions x 1 cell array of idx in each condition
        % listByCondition IS affected by axis randomization, whereas listByConditionRaw is not
        listByCondition
        listByConditionWeights
    end

    properties(Dependent, Transient)
        nTrials

        countByCondition

        nConditionsNonEmpty

        % a mask over trials (T x 1). A trial is valid if all of its attribute values are in the
        % value lists for those attributes, AND manualInvalid(i) == 0
        % logical mask indicating which trials are valid to include when returning groups
        % this mask does not affect any other functions for grabbing attribute values / unique attributes, etc.
        valid

        computedValid

        nValid
    end

    methods % constructor, odc build
        function ci = ConditionInfo()
            ci = ci@ConditionDescriptor(); % calls buildOdc
            ci.odc = ConditionInfoOnDemandCache();
        end
    end

    methods % get / set data stored inside odc
        function v = get.conditionIdx(ci)
            v = ci.odc.conditionIdx;
            if isempty(v)
                ci.odc.conditionIdx = ci.buildConditionIdx();
                v = ci.odc.conditionIdx;
            end
        end

        function ci = set.conditionIdx(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionIdx = v;
        end

        function v = get.conditionSubsRaw(ci)
            v = ci.odc.conditionSubsRaw;
            if isempty(v)
                ci.odc.conditionSubsRaw = ci.buildConditionSubsRaw();
                v = ci.odc.conditionSubsRaw;
            end
        end

        function ci = set.conditionSubsRaw(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionSubsRaw = v;
        end

        function v = get.conditionMembershipRaw(ci)
            v = ci.odc.conditionMembershipRaw;
            if isempty(v)
                [ci.odc.conditionMembershipRaw,  ci.odc.invalidCause] = ci.buildConditionMembershipRaw();
                v = ci.odc.conditionMembershipRaw;
            end
        end

        function ci = set.conditionMembershipRaw(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionMembershipRaw = v;
        end

        function v = get.invalidCause(ci)
            v = ci.odc.invalidCause;
            if isempty(v)
                [ci.odc.conditionMembershipRaw,  ci.odc.invalidCause] = ci.buildConditionMembershipRaw();
                v = ci.odc.invalidCause;
            end
        end

        function ci = set.invalidCause(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.invalidCause = v;
        end

        function v = get.conditionSubs(ci)
            v = ci.odc.conditionSubs;
            if isempty(v)
                ci.odc.conditionSubs = ci.buildConditionSubs();
                v = ci.odc.conditionSubs;
            end
        end

        function ci = set.conditionSubs(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionSubs = v;
        end

        function v = get.conditionMembership(ci)
            v = ci.odc.conditionMembership;
            if isempty(v)
                ci.odc.conditionMembership = ci.buildConditionMembership();
                v = ci.odc.conditionMembership;
            end
        end

        function ci = set.conditionMembership(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionMembership = v;
        end

        function v = get.listByConditionRaw(ci)
            v = ci.odc.listByConditionRaw;
            if isempty(v)
                [ci.odc.listByConditionRaw, ci.odc.listByConditionWeightsRaw] = ci.buildListByConditionRaw();
                v = ci.odc.listByConditionRaw;
            end
        end

        function ci = set.listByConditionRaw(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.listByConditionRaw = v;
        end

        function v = get.listByCondition(ci)
            v = ci.odc.listByCondition;
            if isempty(v)
                [ci.odc.listByCondition, ci.odc.listByConditionWeights] = ci.buildListByCondition();
                v = ci.odc.listByCondition;
            end
        end

        function ci = set.listByCondition(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.listByCondition = v;
        end

        function v = get.listByConditionWeightsRaw(ci)
            v = ci.odc.listByConditionWeightsRaw;
            if isempty(v)
                [ci.odc.listByConditionRaw, ci.odc.listByConditionWeightsRaw] = ci.buildListByConditionRaw();
                v = ci.odc.listByConditionWeightsRaw;
            end
        end

        function ci = set.listByConditionWeightsRaw(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.listByConditionWeightsRaw = v;
        end

        function v = get.listByConditionWeights(ci)
            v = ci.odc.listByConditionWeights;
            if isempty(v)
                [ci.odc.listByCondition, ci.odc.listByConditionWeights] = ci.buildListByCondition();
                v = ci.odc.listByConditionWeights;
            end
        end

        function ci = set.listByConditionWeights(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.listByConditionWeights = v;
        end
    end

    methods % Build data stored inside odc
        function ci = notifyConditionsChanged(ci)
            % when the condition tensor's shape or size changes, we need to
            % reset the condition include mask
            ci.warnIfNoArgOut(nargout);
            ci = notifyConditionsChanged@ConditionDescriptor(ci);

            if ~isempty(ci.conditionMembershipRawManual)
                % flush this with a warning
                warning('TrialData:ConditionDescriptor:FlushConditionMembershipRawManual', 'Flushing conditionMembershipRawManual as conditions have changed');
                ci.conditionMembershipRawManual = [];
                ci.conditionMembershipRawManualReasonInvalid = {};
            end
        end

        function conditionIdx = buildConditionIdx(ci)
            if ci.nTrials > 0
                conditionIdx = TensorUtils.subMat2Ind(ci.conditionsSize, ci.conditionSubs);
            else
                conditionIdx = [];
                return;
            end
        end

        function [matchesFilters, reasonInvalid] = computeTrialsMatchingAttributeFilters(ci)
            valueList = ci.buildAttributeFilterValueListStruct();
            attrFilter = fieldnames(valueList);
            [matchesFilters, whichField] = ci.getAttributeMatchesOverTrials(valueList);

            reasonInvalid = cell(ci.nTrials, 1);
            reasonInvalid(:) = {''};
            explained = false(ci.nTrials, 1);

            for iF = 1:numel(attrFilter)
                reasonInvalid(whichField == iF & ~explained) = {sprintf('invalid attribute value for %s', attrFilter{iF})};
                explained = explained | whichField == iF;
            end
        end

        function [membership, reasonInvalid] = buildConditionMembershipRaw(ci)
            [matchesFilters, reasonInvalid] = ci.computeTrialsMatchingAttributeFilters();
            explained = ~matchesFilters;

            if ~isempty(ci.conditionMembershipRawManual)
                % use the manually specified condition membership raw
                membership = ci.conditionMembershipRawManual;
                reasonInvalid(~explained) = ci.conditionMembershipRawManualReasonInvalid(~explained);
                ci.checkConditionMembershipRaw(membership, matchesFilters);

            else
                % compute the normal way by checking the axis attributes
                if ci.nAxes == 0
                    membership = truevec(ci.nTrials);
                    assert(ci.nConditions == 1);
                    membership(~matchesFilters, :) = false;

                elseif ci.nConditions > 0 && ci.nTrials > 0
                    membership = true([ci.nTrials, ci.conditionsSize]);
                    membership = TensorUtils.assignIntoTensorAlongDimension(membership, false, 1, ~matchesFilters);
                    for iX = 1:ci.nAxes
                        % accept the first axis value that matches
                        % this will be nTrials x nValues that condition
                        matchMatrix = ci.getAttributeMatchesOverTrials(ci.axisValueLists{iX});
                        % mark false any conditions in membership table that
                        % aren't valid by this matchMatrix along axis iX
                        membership = TensorUtils.assignIntoTensorAlongMultipleDimensionsByMask(membership, false, [1 iX+1], ~matchMatrix);
                        tf = ~any(matchMatrix, 2);
                        reasonInvalid(~tf & ~explained) = {sprintf('invalid attributes for axis %s', ci.axisNames{iX})};
                        explained = explained | ~tf;
                    end

                else
                    membership = false(ci.nTrials, ci.nAxes);
                end

                % allow the user function to override values manually
                if ~isempty(ci.conditionMembershipManualFn)
                    [membership, reasonInvalid] = ci.conditionMembershipManualFn(ci, membership, reasonInvalid, matchesFilters);
                    ci.checkConditionMembershipRaw(membership, matchesFilters);
                end
            end

            % remove trials where conditionIncludeMask is false
            if ci.nConditions > 0
                membership(:, ~ci.conditionIncludeMask(:)) = false;
                removeMask = ~any(membership(:, :), 2) & ~explained;
                reasonInvalid(removeMask) = {'condition marked ignored by conditionIncludeMask'};
            end
        end

        function checkConditionMembershipRaw(ci, conditionMembership, matchesFilters)
            % check size
            assert(ismatrix(conditionMembership) && ...
                    size(conditionMembership, 1) == ci.nTrials && ...
                    size(conditionMembership, 2) == ci.nConditions, 'conditionMembership must be nTrials x nConditions matrix');

            % ensure that the condition membership does not include any trials that do not match the attributes
            trialsShouldBeExcluded = find(any(conditionMembership, 2) & ~matchesFilters);

            if ~isempty(trialsShouldBeExcluded)
                error('%d trials are included in conditionMembership matrix but do not pass the attribute filters: %s', ...
                    numel(trialsShouldBeExcluded), strjoin(trialsShouldBeExcluded));
            end
        end

        % mark manual invalid trials as NaNs in all conditionSubs
        % DO NOT PERFORM AXIS RANDOMIZATION HERE
        function membership = buildConditionMembership(ci)
            membership  = ci.conditionMembershipRaw;
            membership(ci.manualInvalid, :) = false;
        end

        function subsMat = buildConditionSubsRaw(ci)
            membership = ci.conditionMembershipRaw;
            someConditionMask = any(membership(:, :) ~= 0, 2);

            if ci.nAxes == 0
                subsMat = onesvec(ci.nTrials);
                assert(ci.nConditions == 1);
                subsMat(~someConditionMask, :) = NaN;

            elseif ci.nConditions > 0 && ci.nTrials > 0
                subsMat = nan(ci.nTrials, ci.nAxes);
                sz = ci.conditionsSizeNoExpand;
                subs = cellvec(numel(sz));
                for iT = 1:ci.nTrials
                    if ~someConditionMask(iT), continue; end
                    % old method, take first. new method, take max since membership can be real-valued now
                    % ind = find(membership(iT, :), 1, 'first');
                    [~, ind] = max(membership(iT, :)); % zero already ruled out above
                    [subs{:}] = ind2sub(sz, ind);
                    subsMat(iT, :) = cell2mat(subs);
                end
            else
                subsMat = nan(ci.nTrials, ci.nAxes);
            end
        end

        %         % compute which condition each trial falls into, without writing
        %         % NaNs for manualInvalid marked trials and without applying randomization along each axis
        %         % do remove trials where conditionIncludeMask == false
        %         function [subsMat, reasonInvalid] = buildConditionSubsRaw(ci)
        %             % filter out any that don't have a valid attribute value
        %             % along the other non-axis attributes which act as a filter
        %             % (i.e. have a manual value list as well). Reason invalid will
        %             % be nTrials x 1 cellstr explaining why a given trial does not
        %             % have a valid set of subs
        %
        %             reasonInvalid = cell(ci.nTrials, 1);
        %             reasonInvalid(:) = {''};
        %             explained = false(ci.nTrials, 1);
        %
        %             valueList = ci.buildAttributeFilterValueListStruct();
        %             attrFilter = fieldnames(valueList);
        %             [matchesFilters, whichField] = ci.getAttributeMatchesOverTrials(valueList);
        %
        %             for iF = 1:numel(attrFilter)
        %                 reasonInvalid(whichField == iF & ~explained) = {sprintf('invalid attribute value for %s', attrFilter{iF})};
        %                 explained = explained | whichField == iF;
        %             end
        %
        %             if ci.nAxes == 0
        %                 subsMat = onesvec(ci.nTrials);
        %                 assert(ci.nConditions == 1);
        %                 subsMat(~matchesFilters, :) = NaN;
        %
        %             elseif ci.nConditions > 0 && ci.nTrials > 0
        %                 subsMat = nan(ci.nTrials, ci.nAxes);
        %                 for iX = 1:ci.nAxes
        %                     % accept the first axis value that matches
        %                     matchMatrix = ci.getAttributeMatchesOverTrials(ci.axisValueLists{iX});
        %                     [tf, match] = max(matchMatrix, [], 2);
        %                     subsMat(tf, iX) = match(tf);
        %
        %                     reasonInvalid(~tf & ~explained) = {sprintf('invalid attributes for axis %s', ci.axisNames{iX})};
        %                     explained = explained | ~tf;
        %                 end
        %
        %                 % mark as NaN if it doesn't match for every attribute
        %                 subsMat(any(subsMat == 0 | isnan(subsMat), 2), :) = NaN;
        %
        %                 subsMat(~matchesFilters, :) = NaN;
        %             else
        %                 subsMat = nan(ci.nTrials, ci.nAxes);
        %             end
        %
        %             % remove trials where conditionIncludeMask is false
        %             if ci.nConditions > 0
        %                 conditionIdx = TensorUtils.subMat2Ind(ci.conditionsSize, subsMat); %#ok<*PROP>
        %                 removeMask = ~isnan(conditionIdx); % only consider trials still valid
        %                 removeMask(~isnan(conditionIdx)) = ~ci.conditionIncludeMask(conditionIdx(~isnan(conditionIdx)));
        %                 subsMat(removeMask, :) = NaN;
        %                 reasonInvalid(removeMask) = {'condition marked ignored by conditionIncludeMask'};
        %             end
        %         end

        function valueList = buildAttributeFilterValueListStruct(ci)
            mask = ci.attributeActsAsFilter;
            names = ci.attributeNames(mask);
            vals = ci.attributeValueLists(mask);
            valueList = struct();
            for iA = 1:numel(names)
                valueList.(names{iA}) = vals{iA};
            end
        end

        % mark manual invalid trials as NaNs in all conditionSubs
        % DO NOT PERFORM AXIS RANDOMIZATION HERE
        function subsMat = buildConditionSubs(ci)
            subsMat = ci.conditionSubsRaw;
            subsMat(ci.manualInvalid, :) = NaN;
        end

        function [list, weights] = buildListByConditionRaw(ci)
            % DO NOT PERFORM AXIS RANDOMIZATION OR SORTING HERE
            [list, weights] = deal(cell(ci.conditionsSize));
            for iC = 1:ci.nConditions
                % this allows for a condition to belong to multiple
                % conditions
                list{iC} = find(ci.conditionMembership(:, iC));
                weights{iC} = ci.conditionMembership(list{iC}, iC);
                if isempty(list{iC})
                    % ensure it can be concatenated into a column
                    % vector using cell2mat
                    list{iC} = nan(0, 1);
                    weights{iC} = nan(0, 1);
                end
            end
        end

        function list = generateSingleRandomizedListByCondition(ci, list, seed)
            if all(ci.axisRandomizeModes == ci.AxisOriginal) && ~ci.isResampledWithinConditions
                % no randomization to be done
                return;
            end

            % seed the random number generator exactly once
            ci.seedRandStream(seed);

            for iA = 1:ci.nAxes
                switch ci.axisRandomizeModes(iA)
                    case ci.AxisOriginal
                        continue;
                    case ci.AxisShuffled
                        replace = ci.axisRandomizeWithReplacement(iA);
                        list = TensorUtils.listShuffleAlongDimension(list, iA, replace);
                    case ci.AxisResampledFromSpecified

                        % build lookup index list in value list based on
                        % struct match
                        nValues = ci.nValuesAlongAxes(iA);
                        idxResampleFromList = cellvec(nValues);
                        for iV = 1:nValues
                            idxResampleFromList{iV} = ci.axisLookupValueInValueList(iA, ci.axisRandomizeResampleFromList{iA}{iV});
                        end

                        replace = true; % replace must be true since we will be sampling from a different condition with a different trial count
                        list = TensorUtils.listResampleFromSpecifiedAlongDimension(list, idxResampleFromList, iA, replace);
                    otherwise
                        error('Unknown randomize mode for axis %d', iA);
                end
            end

            % and finally, if resampleWithinConditions is true,
            % resampleFromSame everything
            if ci.isResampledWithinConditions
                list = TensorUtils.listResampleFromSame(list);
            end
        end


        function [list, weights] = buildListByCondition(ci)
            % Take .listByConditionRaw and perform axis randomization /
            % resampling and then sort lists by attributes
            % Perform axis randomization to listByConditionRaw
            % successively along each randomized axis

            listWithWeights = TensorUtils.listHorzCatContents(ci.listByConditionRaw, ci.listByConditionWeightsRaw);

            listWithWeights = ci.generateSingleRandomizedListByCondition(listWithWeights, ci.randomSeed);

            [list, weights] = TensorUtils.listSplitContents(listWithWeights);
            [list, weights] = ci.sortListByConditionByAttributes(list, weights);
        end

        function listCell = generateMultipleRandomizedListByCondition(ci, varargin)
            % When axis randomization is applied, generate a specific number of
            % listByCondition cells (i.e. cell tenors containing lists of trial indexes)
            % by using successive integer random seeds. This is useful when
            % performing statistical tests using axisRandomization techniques.
            % listCell is nConditions x nSamples
            p = inputParser();
            p.addOptional('n', 100, @isscalar);
            p.addParameter('initialSeed', ci.randomSeed, @(x) isscalar(x));
            p.addParameter('showProgress', true, @islogical);
            p.parse(varargin{:});

            assert(ci.hasRandomization, 'CondtionInfo needs some kind of randomization applied');

            N = p.Results.n;
            initialSeed = p.Results.initialSeed;
            showProgress = p.Results.showProgress;

            if showProgress
                prog = ProgressBar(N, 'Building n=%d randomized trial lists', N);
            end

            listOriginal = ci.listByConditionRaw;

            listCell = cell(ci.nConditions, N);
            for i = 1:N
                if showProgress
                    prog.update(i);
                end
                listTensor = ci.generateSingleRandomizedListByCondition(listOriginal, initialSeed+i-1);
                listCell(:, i) = listTensor(:);
            end
            if showProgress
                prog.finish();
            end
        end

        function [listCell, weightCell] = sortListByConditionByAttributes(ci, listCell, weightCell)
            sortByList = ci.attributeSortByList;
            if isempty(sortByList)
                return;
            end

            nSortAttr = numel(sortByList);
            attr = cell(nSortAttr, 1);
            rev = false(nSortAttr, 1);
            for i = 1:nSortAttr
                if strncmp(sortByList{i}, '-', 1)
                    attr{i} = sortByList{i}(2:end);
                    rev(i) = true;
                else
                    attr{i} = sortByList{i};
                    rev(i) = false;
                end
            end

            % sort over attributes in reverse order so that first takes preference
            for iA = nSortAttr:-1:1
                attrValsFull = ci.getAttributeValues(attr{iA});
                if rev(i)
                    mode = 'descend';
                else
                    mode = 'ascend';
                end

                for iC = 1:numel(listCell)
                    attrVals = attrValsFull(listCell{iC});
                    if iscellstr(attrVals) %#ok<*ISCLSTR>
                        [~, idx] = sortStrings(attrVals, mode);
                    elseif isnumeric(attrVals) && isvector(attrVals)
                        [~, idx] = sort(makecol(attrVals), 1, mode);
                    end

                    listCell{iC} = listCell{iC}(idx);
                    weightCell{iC} = weightCell{iC}(idx);
                end
            end
        end
    end

    methods % ConditionDescriptor overrides and utilities for auto list generation
        function printOneLineDescription(ci)
            if ci.nAxes == 0
                axisStr = 'no grouping axes';
            else
                axisStr = strjoin(ci.axisDescriptions, ', ');
            end

            attrFilter = ci.attributeNames(ci.attributeActsAsFilter);
            if isempty(attrFilter)
                filterStr = 'no filtering';
            else
                filterStr = sprintf('filtering by %s', strjoin(attrFilter));
            end

            validStr = sprintf('(%d valid)', nnz(ci.computedValid));

            tcprintf('inline', '{yellow}%s: {none}%s, %s %s\n', ...
                class(ci), axisStr, filterStr, validStr);
        end

        function valueList = buildAttributeValueLists(ci)
            if ~ci.applied
                % act like ConditionDescriptor before applied to trial data
                valueList = buildAttributeValueLists@ConditionDescriptor(ci);
                return;
            end

            % figure out the automatic value lists
            modes = ci.attributeValueModes;
            valueList = cellvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                switch modes(i)
                    case ci.AttributeValueListAuto
                        % compute unique bins
                        valueList{i} = ci.computeAutoListForAttribute(i);

                    case ci.AttributeValueListManual
                        % use manual list
                        valueList{i} = ci.attributeValueListsManual{i};

                    case ci.AttributeValueBinsManual
                        % use specified bins
                        valueList{i} = ci.attributeValueBinsManual{i};

                    case ci.AttributeValueBinsAutoUniform
                        % compute bin boundaries
                        valueList{i} = ci.computeAutoUniformBinsForAttribute(i);

                    case ci.AttributeValueBinsAutoQuantiles
                        % compute bin boundaries
                        valueList{i} = ci.computeAutoQuantileBinsForAttribute(i);
                end
                valueList{i} = makecol(valueList{i});
            end
        end

        function valueList = computeAutoListForAttribute(ci, attrIdx)
            % use tolerance to ensure we don't split up conditions due to
            % floating point imprecision

            vals = ci.getAttributeValues(attrIdx);
            if ci.attributeNumeric(attrIdx)
                if all(isnan(vals))
                    valueList = NaN;
                elseif islogical(vals)
                    valueList = unique(vals);
                else
                    valueList = TrialDataUtilities.Data.uniquetol(removenan(vals));
                    % include NaN in the list if one is found
                    if any(isnan(vals))
                        valueList(end+1) = NaN;
                    end
                end
            elseif iscellstr(vals) || iscategorical(vals) || isstring(vals)
                valueList = unique(vals);
            else
                valueList = TrialDataUtilities.Data.uniqueCellTol(vals);
            end
        end

        function bins = computeAutoUniformBinsForAttribute(ci, attrIdx)
            vals = cell2mat(ci.values(:, attrIdx));
            nBins = ci.attributeValueBinsAutoCount(attrIdx);
            minV = nanmin(vals);
            maxV = nanmax(vals);

            if isnan(minV) || isnan(maxV) || isnan(nBins)
                bins = [NaN, NaN];
            else
                binEdges = makecol(linspace(minV, maxV, nBins + 1));
                bins = [ binEdges(1:end-1), binEdges(2:end) ];
            end

            bins = mat2cell(bins, ones(size(bins, 1), 1), 2);
        end

        function bins = computeAutoQuantileBinsForAttribute(ci, attrIdx)
            vals = removenan(cell2mat(ci.values(:, attrIdx)));
            nBins = ci.attributeValueBinsAutoCount(attrIdx);

            if isempty(vals)
                bins = [NaN, NaN];
            else
                binEdges = makecol(quantile(vals, linspace(0, 1, nBins+1)));
                bins = [ binEdges(1:end-1), binEdges(2:end) ];
            end

            bins = mat2cell(bins, ones(size(bins, 1), 1), 2);
        end

        function valueListAsStrings = buildAttributeValueListsAsStrings(ci)
            if ~ci.applied
                % act like ConditionDescriptor before applied to trial data
                valueListAsStrings = buildAttributeValueListsAsStrings@ConditionDescriptor(ci);
                return;
            end

            % rely on ConditionDescriptor's implementation, substitute
            % where necessary
            modes = ci.attributeValueModes;
            valueListAsStrings = buildAttributeValueListsAsStrings@ConditionDescriptor(ci);
            valueListsDisplayAs = ci.attributeValueListsAsStringsManual;
            valueList = ci.attributeValueLists;

            for i = 1:ci.nAttributes
                if isempty(ci.attributeUnits{i})
                    unitsStr = '';
                else
                    unitsStr = [' ', ci.attributeUnits{i}];
                end
                if modes(i) == ci.AttributeValueListAuto
                    if ~isempty(valueListsDisplayAs{i})
                        displayAs = makecol(valueListsDisplayAs{i});
                        assert(numel(valueList{i}) == numel(displayAs), 'attributeValueListsAsStringsManual for attribute %s has the wrong number of entries');
                        valueListAsStrings{i} = displayAs;
                    else
                        % convert populated list to cellstr
                        if ~iscell(valueList{i})
                            if isstring(valueList{i}) || iscategorical(valueList{i})
                            valueListAsStrings{i} = arrayfun(@char, valueList{i}, 'UniformOutput', false);
                            else
                                valueListAsStrings{i} = arrayfun(@(x) [num2str(x), unitsStr], valueList{i}, 'UniformOutput', false);
                            %                             valueListAsStrings{i} = arrayfun(@(x) num2str(x), valueList{i}, 'UniformOutput', false);
                            end
                        elseif iscellstr(valueList{i})
                            valueListAsStrings{i} = valueList{i}; % cellfun(@(x) [x, unitsStr], valueList{i}, 'UniformOutput', false);
                            %                             valueListAsStrings{i} = [valueList{i}, unitsStr];
                        else
                            error('Not sure how to convert attribute values to string');
                        end
                    end
                end
                valueListAsStrings{i} = makecol(valueListAsStrings{i});
            end
        end

        function valueListByAxes = buildAxisValueLists(ci)
            % uses ConditionDescriptor's implementation but deals with
            valueListByAxes = buildAxisValueLists@ConditionDescriptor(ci);
            if ~ci.applied
                return;
            end

            % exclude trials which would be invalid because of attribute value list filters
            validTrials = ~ci.manualInvalid;
            valueList = ci.buildAttributeFilterValueListStruct();
            matchesFilters = ci.getAttributeMatchesOverTrials(valueList);
            validTrials(~matchesFilters) = false;

            % and exclude trials with invalid trials on the *manually*
            % specified axes
            for iX = 1:ci.nAxes
                switch ci.axisValueListModes(iX)
                    case ci.AxisValueListManual
                        matchMatrix = ci.getAttributeMatchesOverTrials(ci.axisValueListsManual{iX});
                        validTrials(~any(matchMatrix, 2)) = false;
                end
            end

            for iX = 1:ci.nAxes
                % build a cellstr of descriptions of the values along this axis
                switch ci.axisValueListModes(iX)
                    case ci.AxisValueListAutoOccupied
                        % need to filter by which values are actually
                        % occupied by at least one trial
                        maskTrialsByValues = ci.getAttributeMatchesOverTrials(valueListByAxes{iX});
                        % exclude invalid or soon-to-be invalid trials as
                        % marked above
                        maskTrialsByValues(~validTrials, :) = false;
                        keepValues = any(maskTrialsByValues, 1);

                        valueListByAxes{iX} = makecol(valueListByAxes{iX}(keepValues));
                end
            end
        end

        function [mask, whichField] = getAttributeMatchesOverTrials(ci, valueStruct)
            % valueStruct is a struct where .attribute = [vals] or {vals}
            % matches trials where attribute takes a value in vals
            % return a logical mask nTrials x 1 indicating these matches
            % if valueStruct is a length nValues struct vector, mask will
            % be nTrials x nValues. whichField will be nTrials x nValues
            % indiciating which field invalidated a given trial (or NaN)

            if ci.nTrials == 0
                mask = false(0, 1);
                whichField = nan(0, 1);
                return;
            end

            nValues = numel(valueStruct);
            mask = true(ci.nTrials, nValues);
            whichField = nan(ci.nTrials, nValues);

            if nValues == 0
                % no valid values provided, often because it's
                % an auto-occupied axis but no trials are valid
                mask = false(ci.nTrials, 1);
                return;
            end

            fields = fieldnames(valueStruct);
            attrIdx = ci.assertHasAttribute(fields);

            for iF = 1:numel(fields) % loop over attributes to match
                attrVals = ci.getAttributeValues(attrIdx(iF));
                switch ci.attributeValueModes(attrIdx(iF))

                    case {ci.AttributeValueListAuto, ci.AttributeValueListManual}
                        % match against value lists
                        for iV = 1:nValues % loop over each value in value list
                            valsThis = valueStruct(iV).(fields{iF});

                            if isempty(valsThis)
                                % empty is a wildcard match
                                continue;
                            end

                            % check whether value list has sublists within
                            % and flatten them if so
                            if ci.attributeNumeric(attrIdx(iF))
                                if iscell(valsThis)
                                    valsThis = [valsThis{:}];
                                    % groups of values per each element
                                end
                            else
                                % non-numeric
                                if iscell(valsThis) && ~iscellstr(valsThis) && ~ischar(valsThis)
                                    valsThis = [valsThis{:}];
                                end
                            end

                            matchesThis = TrialDataUtilities.Data.ismembertol(attrVals, valsThis);
                            mask(:, iV) = mask(:, iV) & matchesThis;

                            whichField(isnan(whichField(:, iV)) & ~matchesThis) = iF;
                        end

                    case {ci.AttributeValueBinsManual, ci.AttributeValueBinsAutoUniform, ...
                            ci.AttributeValueBinsAutoQuantiles}
                        % match against bins. valueStruct.attr is nBins x 2 bin edges
                        for iV = 1:nValues
                            valsThis = valueStruct(iV).(fields{iF});

                            if isempty(valsThis)
                                % empty is a wildcard match
                                continue;
                            end
                            matchesThis = matchAgainstBins(attrVals, valsThis);
                            mask(:, iV) = mask(:, iV) & matchesThis;
                            whichField(isnan(whichField(:, iV)) & ~matchesThis) = iF;
                        end
                end
            end

            function binAccept = matchAgainstBins(vals, bins)
                bins = cell2mat(bins);
                binAccept = any(bsxfun(@ge, vals, bins(:, 1)') & bsxfun(@le, vals, bins(:, 2)'), 2);
            end
        end

        function values = getAttributeValues(ci, name)
            idx = ci.getAttributeIdx(name);
            values = ci.values(:, idx);
            if ci.attributeAsVector(idx)
                values = cat(1, values{:});
            end
        end

        function ci = maskAttributes(ci, mask)
            ci.warnIfNoArgOut(nargout);
            ci.values = ci.values(:, mask);
            ci = maskAttributes@ConditionDescriptor(ci, mask);
        end
    end

    methods % Trial utilities and dependent properties
        function counts = get.countByCondition(ci)
            counts = cellfun(@length, ci.listByCondition);
        end

        function nConditions = get.nConditionsNonEmpty(ci)
            nConditions = nnz(~cellfun(@isempty, ci.listByCondition));
        end

        function nt = get.nTrials(ci)
            nt = size(ci.values, 1);
        end

        % mark additional trials invalid
        function ci = markInvalid(ci, invalid)
            ci.warnIfNoArgOut(nargout);
            ci.manualInvalid(invalid) = true;
            ci = ci.invalidateCache();
        end

        % overwrite manualInvalid with invalid, ignoring what was already
        % marked invalid
        function ci = setManualValidTo(ci, valid)
            % only invalidate if changing
            ci.warnIfNoArgOut(nargout);
            assert(isvector(valid) & numel(valid) == ci.nTrials, 'Size mismatch');
            if isempty(valid)
                return;
            end
            if any(ci.manualInvalid ~= ~valid)
                ci.manualInvalid = makecol(~valid);
                ci = ci.invalidateCache();
            end
        end

        function valid = get.valid(ci)
            % return a mask which takes into account having a valid value for each attribute
            % specified, as well as the markInvalid function which stores its results in .manualInvalid
            valid = ~ci.manualInvalid & ci.computedValid;
        end

        function computedValid = get.computedValid(ci)
            if ci.nTrials > 0
                computedValid = all(~isnan(ci.conditionSubsRaw), 2);
            else
                computedValid = false(0, 1);
            end
        end

        function nValid = get.nValid(ci)
            nValid = nnz(ci.valid);
        end

        function mask = getIsTrialInSomeGroup(ci)
            mask = ~isnan(ci.conditionIdx);
        end

        function ci = selectTrials(ci, selector)
            ci.warnIfNoArgOut(nargout);

            assert(isvector(selector), 'Selector must be vector of indices or vector mask');
            % cache everything ahead of time because some are dynamically
            % computed from the others

            ci.manualInvalid = ci.manualInvalid(selector);
            ci.values = ci.values(selector, :);
            ci = ci.invalidateCache();
        end

        % given a cellvec or nmeric vector, group its elements
        function varargout = groupElements(ci, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == ci.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) data(idx,:), ci.listByCondition, ...
                    'UniformOutput', false);
            end
        end

        % like groupElements, but returns flattened nConditions x 1 cells
        % rather than tensors in the shape of conditionsSize
        function varargout = groupElementsFlattened(ci, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == ci.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) data(idx,:), ci.listByCondition, ...
                    'UniformOutput', false);
                varargout{i} = varargout{i}(:);
            end
        end
    end

    methods % Apply to trial data
        function ci = initializeWithNTrials(ci, N)
            ci.warnIfNoArgOut(nargout);
            % build empty arrays for N trials
            ci.manualInvalid = false(N, 1);
            ci.values = cell(N, ci.nAttributes);
        end

        function ci = applyToTrialData(ci, td)
            % build the internal attribute value list (and number of trials)
            % from td.
            ci.warnIfNoArgOut(nargout);

            % set trialCount to match length(trialData)
            nTrials = ci.getNTrialsFn(td); %#ok<*PROPLC>
            ci = ci.initializeWithNTrials(nTrials);

            if ci.nAttributes > 0 && ci.nTrials > 0
                % fetch valuesByAttribute using callback function
                valueStruct = ci.requestAttributeValues(td, ci.attributeNames);
                valueCell = struct2cell(valueStruct)';

                % store in .values cell
                ci.values = valueCell;

                % fetch other details as well, units, and numeric status
                for iA = 1:ci.nAttributes
                    ci.attributeUnits{iA} = td.getChannelUnitsPrimary(ci.attributeNames{iA});
                    ci = ci.setAttributeNumeric(iA, td.isChannelNumericScalar(ci.attributeNames{iA}));
                    % TODO change this back to isChannelScalar
                    ci = ci.setAttributeAsVector(iA, td.isChannelScalar(ci.attributeNames{iA}));
%                     ci.attributeNumeric(iA) = td.isChannelScalar(ci.attributeNames{iA});
                end

                ci = ci.fixAttributeValues();
            end

            ci.applied = true;
            ci = ci.invalidateCache();
        end

        function ci = fixAttributeValues(ci, attrIdx)
            ci.warnIfNoArgOut(nargout);
            if ci.nAttributes == 0 || ci.nTrials == 0
                return;
            end

            if nargin < 2
                % go over all attributes if not specified
                attrIdx = 1:ci.nAttributes;
            end

            for iList = 1:numel(attrIdx)
                i = attrIdx(iList);
                vals = ci.values(:, i);

                % check for numeric, replace empty with NaN
                emptyMask = cellfun(@isempty, vals);
                vals(emptyMask) = {NaN};

                % moved this to applyToTrialData
                %                 isNumeric = cellfun(@(x) isscalar(x) && (islogical(x) || isnumeric(x)), vals);

                %                 if all(isNumeric)
                if ci.attributeNumeric(i)
                    %                     ci.attributeNumeric(i) = true;
                    isboolean = cellfun(@(x) isscalar(x) && islogical(x), vals);
                    if all(isboolean)
                        mat = cell2mat(vals);
                    else
                        % convert to double type
                        mat = cellfun(@double, vals);
                    end
                    assert(numel(vals) == numel(mat));
                    ci.values(:, i) = num2cell(mat);
                    %                     ci.attributeNumeric(i) = true;
                elseif ci.attributeAsVector(i)
                    % nothing to do for categorical or strings
                else
                    % replace empty and NaN with '' (NaN for strings)
                    nanMask = cellfun(@(x) any(ismissing(x)), vals);
                    vals(nanMask) = {''};

                    %                     ci.attributeNumeric(i) = false;

                    % check for cellstr
                    %if iscellstr(vals)
                    ci.values(:, i) = vals;
                    %                         ci.attributeNumeric(i) = false;
                    %else
                    %   error('Attribute %s values were neither uniformly scalar nor strings', ci.attributeNames{i});
                    %end
                end
            end
        end

        function ci = setAttributeValueData(ci, name, dataCell)
            if ~iscell(dataCell)
                dataCell = num2cell(dataCell);
            end
            assert(numel(dataCell) == ci.nTrials, 'Data must be size nTrials');

            idx = ci.assertHasAttribute(name);
            ci.values(:, idx) = dataCell;

            ci = ci.fixAttributeValues(idx);

            ci.invalidateCache();
        end

        function assertNotApplied(ci)
            if ci.applied
                error('You must unbind this ConditionInfo before adding attributes');
            end
        end

        function ci = addAttribute(ci, name, varargin)
            ci.warnIfNoArgOut(nargout);

            if ci.applied
                % ensure values are specified if already applied
                % since we won't be requesting them
                p = inputParser;
                p.KeepUnmatched = true;
                p.addParameter('values', {}, @(x) islogical(x) || isnumeric(x) || iscell(x) || iscategorical(x) || isstring(x));
                p.parse(varargin{:});

                if ismember('values', p.UsingDefaults)
                    error('This ConditionInfo has already been applied to data. values must be specified when adding new attributes');
                end

                % add via ConditionDescriptor
                ci = addAttribute@ConditionDescriptor(ci, name, p.Unmatched);

                % set the values in my .values cell array
                vals = p.Results.values;
                assert(numel(vals) == ci.nTrials, ...
                    'Values provided for attribute must have numel == nTrials');

                iAttr = ci.nAttributes;
                % critical to update attribute numeric here!

                % TODO fix this to deal with string arrays correctly
                if isstring(vals)
                    vals = cellstr(vals);
                end
                if iscell(vals)
                    ci.attributeNumeric(iAttr) = false;
                    ci.attributeAsVector(iAttr) = false;
                    ci.values(:, iAttr) = vals;
                else
                    ci.attributeAsVector(iAttr) = true;
                    if iscategorical(vals) || isstring(vals)
                        ci.attributeNumeric(iAttr) = false;
                    else
                        ci.attributeNumeric(iAttr) = true;
                    end
                    ci.values(:, iAttr) = num2cell(vals);
                end

                % fix everything up and rebuild the caches
                ci = ci.fixAttributeValues();
                ci = ci.invalidateCache();
            else
                % if not applied, no need to do anything special
                ci = addAttribute@ConditionDescriptor(ci, name, varargin{:});
            end
        end

        function valueStruct = requestAttributeValues(ci, td, attrNames)
            % lookup requestAs name if not specified

            % translate into request as names
            if ci.getNTrialsFn(td) == 0
                valueStruct = struct();
            else
                valueStructRequestAs = ci.getAttributeValueFn(td, attrNames);

                % check the returned size and field names
                assert(numel(valueStructRequestAs) == ci.nTrials, 'Number of elements returned by getAttributeFn must match nTrials');
                if ~all(isfield(valueStructRequestAs, attrNames))
                    missingFields = attrNames(~isfield(valueStructRequestAs, attrNames));
                    error('TrialData missing fields required by ConditionInfo: %s', strjoin(missingFields, ', '));
                end

                % translate back into attribute names
                %                 valueStruct = mvfield(valueStructRequestAs, requestAs, attrNames);
                %                 valueStruct = orderfields(valueStruct, attrNames);
                valueStruct = makecol(valueStructRequestAs);
            end
        end
    end

    methods % Manual conditionMembership specification and weighted trials within conditions
        function ci = setConditionMembershipManual(ci, conditionMembership, reasonInvalid)
            ci.assertAllAxisValueListsManual();

            assert(ismatrix(conditionMembership) && ...
                    size(conditionMembership, 1) == ci.nTrials && ...
                    size(conditionMembership, 2) == ci.nConditions, 'conditionMembership must be nTrials x nConditions matrix');

            if nargin < 2 || isempty(reasonInvalid)
                reasonInvalid = 'excluded via conditionMembershipManual';
            end

            if ischar(reasonInvalid)
                reasonInvalid = repmat({reasonInvalid}, ci.nTrials, 1);
            end
            assert(iscellstr(reasonInvalid), 'Must be cellstr or string');
            reasonInvalid = makecol(reasonInvalid);

            % clear out reasons when actually in a condition
            maskInvalidated = ~any(conditionMembership ~= 0, 2);
            reasonInvalid(~maskInvalidated) = {''};

            ci.conditionMembershipRawManual = conditionMembership;
            ci.conditionMembershipRawManualReasonInvalid = reasonInvalid;
            ci = ci.invalidateCache();

            % get conditionMembership to force recomputation, which will do some error checking
            ci.conditionMembership;
        end

        function ci = setConditionMembershipManualFn(ci, conditionMembershipFn)
            ci.warnIfNoArgOut(nargout);

            assert(isa(conditionMembershipFn, 'function_handle'));
            ci.conditionMembershipManualFn = conditionMembershipFn;
            ci = ci.invalidateCache();
        end

        function [membership, reasonInvalid] = testConditionMembershipManualFn(ci)
            ci.warnIfNoArgOut(nargout);

            % same as build appearances but does not catch errors to allow
            assert(~isempty(ci.conditionMembershipManualFn), 'No conditionMembershipManualFn is set');
            [membership, reasonInvalid] = ci.buildConditionMembershipRaw();
        end
    end

    methods
        % same as ConditionDescriptor, except skips conditions with no
        % trials so that the colors stay maximally separated
        %         function a = defaultAppearanceFn(ci, varargin)
        %             % returns a struct specifying the default set of appearance properties
        %             % for the given group. indsGroup is a length(ci.groupByList) x 1 array
        %             % of the inds where this group is located in the high-d array, and dimsGroup
        %             % gives the full dimensions of the list of groups.
        %             %
        %             % We vary color along all axes simultaneously, using the linear
        %             % inds.
        %             %
        %             % Alternatively, if no arguments are passed, simply return a set of defaults
        %             nConditionsNonEmpty = ci.nConditionsNonEmpty;
        %             countsByCondition = ci.countByCondition;
        %
        %             nConditions = ci.nConditions;
        %
        %             a(ci.conditionsSize()) = AppearanceSpec();
        %
        %             if nConditions == 1
        %                 cmap = [0.3 0.3 1];
        %             else
        %                 if nConditions > 256
        %                     cmap = jet(nConditions);
        %                 else
        %                     cmap = distinguishable_colors(nConditionsNonEmpty);
        %                 end
        %             end
        %
        %             colorInd = 1;
        %             for iC = 1:nConditions
        %                 if countsByCondition(iC) > 0
        %                     a(iC).Color = cmap(colorInd, :);
        %                     colorInd = colorInd + 1;
        %                 end
        %             end
        %         end
    end

    methods % Convert value lists to manual
        function ci = fixAttributeValueList(ci, name)
            ci.warnIfNoArgOut(nargout);
            iAttr = ci.assertHasAttribute(name);

            % cache this since it will be reset
            ciMask = ci.conditionIncludeMaskManual;

            switch(ci.attributeValueModes(iAttr))
                case {ci.AttributeValueListManual, ci.AttributeValueBinsManual}
                    % it's manual already
                    return;
                case ci.AttributeValueListAuto
                    ci = ci.setAttributeValueList(name, ci.attributeValueLists{iAttr});
                case {ci.AttributeValueBinsAutoUniform, ci.AttributeValueBinsAutoQuantiles}
                    ci = ci.binAttribute(name, ci.attributeValueLists{iAttr});
                otherwise
                    error('Unknown attributeValueList mode');
            end

            ci = ci.invalidateCache();
            ci.conditionIncludeMaskManual = ciMask;
        end

        function ci = fixAllAttributeValueLists(ci)
            ci.warnIfNoArgOut(nargout);
            for iA = 1:ci.nAttributes
                ci = ci.fixAttributeValueList(iA);
            end
        end

        function ci = fixAxisValueList(ci, axisSpec)
            ci.warnIfNoArgOut(nargout);
            % cache this since it will be reset
            ciMask = ci.conditionIncludeMaskManual;

            idx = ci.axisLookupByAttributes(axisSpec);
            for i = 1:numel(idx)
                ci = ci.setAxisValueList(idx(i), ci.axisValueLists{idx(i)});
            end

            ci.conditionIncludeMaskManual = ciMask;
        end

        function ci = fixAllAxisValueLists(ci)
            ci.warnIfNoArgOut(nargout);
            for iA = 1:ci.nAxes
                ci = ci.fixAxisValueList(iA);
            end
        end

        function ci = fixAllValueLists(ci)
            ci.warnIfNoArgOut(nargout);
            ci = ci.fixAllAttributeValueLists();
            ci = ci.fixAllAxisValueLists();
        end
    end

    methods % Conversion to ConditionDescriptor

        % build a static ConditionDescriptor with the same specs as this
        % ConditionInfo
        function cd = getConditionDescriptor(ci, varargin)
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        % build a static ConditionDescriptor with the same specs as this
        % ConditionInfo WITH ALL VALUE LISTS AND BINS FIXED
        function cd = getFixedConditionDescriptor(ci)
            cd = ci.fixAllValueLists().getConditionDescriptor();
        end
    end

    methods(Static) % From condition descriptor and default callbacks
        % Building from a condition descriptor with an accessor method
        function ci = fromConditionDescriptor(cd, varargin)
            p = inputParser;
            p.addOptional('trialData', [], @(x) true);
            p.parse(varargin{:});

            % build up the condition info
            ci = ConditionInfo();

            % Have conditionDescriptor copy over the important details
            ci = ConditionDescriptor.fromConditionDescriptor(cd, ci);

            ci.odc = ConditionInfoOnDemandCache();

            % and then apply to the trialData
            if ~isempty(p.Results.trialData)
                ci = ci.applyToTrialData(p.Results.trialData);
            end
        end

        % construct condition descriptor from a struct of attribute values
        function cd = fromStruct(s)
            cd = ConditionInfo();
            cd = cd.addAttributes(fieldnames(s));
            cd = cd.applyToTrialData(s);
        end

        % return a scalar struct with one field for each attribute containing the attribute values
        % as a cell array or numeric vector
        function values = defaultGetAttributeFn(data, attributeNames, varargin)
            assert(isstruct(data) || isa(data, 'TrialData'), 'Please provide getAttributeFn if data is not struct or TrialData');

            if isstruct(data)
                % TODO implement request as renaming here
                values = keepfields(data, attributeNames);
                %                 for iAttr = 1:length(attributeNames)
                %                     attr = attributeNames{iAttr};
                %                     valuesByAttribute.(attr) = {data.(attr)};
                %                 end
            elseif isa(data, 'TrialData')
                values = data.getRawChannelDataAsStruct(attributeNames);
                %values = keepfields(data.getParamStruct, attributeNames);
            else
                error('Please provide .getAttributeFn to request attributes from this data type');
            end
        end

        function nTrials = defaultGetNTrialsFn(data, varargin)
            if isempty(data)
                nTrials = 0;
                return;
            end

            assert(isstruct(data) || isa(data, 'TrialData'), 'Please provide getNTrialsFn if data is not struct or TrialData');
            if isstruct(data)
                nTrials = numel(data);
            else
                nTrials = data.nTrials;
            end
        end
    end

end
