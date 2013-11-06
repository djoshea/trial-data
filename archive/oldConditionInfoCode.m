
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
        
        
            
%             % exclude any values not found in the valueList, if specified
%             valueList = ci.getAttributeValueList(name);
%             if ~isempty(valueList)
%                 invalidValues = ~ismember(values, valueList);
%                 if iscellstr(values)
%                     [values{invalidValues}] = deal('');
%                 elseif iscell(values)
%                     [values{invalidValues}] = deal(NaN);
%                 else
%                     values(invalidValues) = NaN;
%                 end
%             end
%     end
% 
%         function values = getAttributeUnique(ci, name)
%             % if a value list is specified, we simply return that
%             % otherwise, return the unique list of values
%             valueList = ci.getAttributeValueList(name);
%             if ~isempty(valueList)
%                 values = valueList; 
%             else
%                 values = unique(ci.getAttribute(name));
% 
%                 % remove empty and nan values
%                 if isnumeric(values) || islogical(values)
%                     remove = isnan(values) | isempty(values);
%                     values = num2cell(values);
%                 else
%                     remove = cellfun(@isempty, values);
%                 end
%                 values(remove) = [];
%             end
%         end
% 
%         function valueCell = getMultipleAttributeUnique(ci, names)
%            valueCell = cellfun(@(name) ci.getAttributeUnique(name), names, ...
%                 'UniformOutput', false);
%         end
%     
%         function valueCell = getAllUnique(ci)
%             valueCell = ci.getMultipleAttributeUnique(ci.attributeNames);
%         end
% 
%         function idx = getIdxWithAttributeValue(ci, name, value)
%             values = ci.getAttribute(name);
%             if isnumeric(values) || islogical(values)
%                 values = num2cell(values);
%             end
%             match = cellfun(@(x) isequal(x, value), values);
%             idx = find(match);
%         end
% 
%         function idxCell = getIdxEachAttributeValue(ci, name)
%             values = ci.getAttributeUnique(name);
%             nValues = length(values);
%             idxCell = cell(nValues, 1);
% 
%             if isnumeric(values) || islogical(values)
%                 values = num2cell(values);
%             end
%             for iValue = 1:nValues
%                 idxCell{iValue} = ci.getIdxWithAttributeValue(name, values{iValue});
%             end
%         end
%         
%         function idxList = getAttributeAsIdxUnique(ci, name)
%             % return a numeric vector of the attribute value for trial i
%             % like getAttribute, but instead of the raw value return an
%             % index into getAttributeUnique(name)
%             
%             idxList = nan(ci.nTrials, 1);
%             idxCell = ci.getIdxEachAttributeValue(name);
%             for i = 1:length(idxCell);
%                 idxList(idxCell{i}) = i;
%             end
%         end
    end