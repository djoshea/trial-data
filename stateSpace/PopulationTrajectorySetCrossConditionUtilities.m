classdef PopulationTrajectorySetCrossConditionUtilities
   
    % combining / weighting along condition axes
    methods(Static)
        function psetDiff = computeDifferenceAlongAxis(pset, axisName, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('autoNamesAlongAxis', true, @islogical);
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.addParameter('newNamesShortAlongAxis', {}, @iscellstr);
            p.addParameter('reverse', false, @islogical); % default is 2-1, reverse is 1-2
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            reverse = p.Results.reverse;
            nAlongAxis = pset.conditionsSize(aIdx);
            
            % generate new names from differences
            if p.Results.autoNamesAlongAxis
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                stringList = stringLists{aIdx};
                if ~reverse
                    newNamesAlongAxis = cellfun(@(v1, v2) [v2 ' - ' v1], stringList(1:end-1), ...
                        stringList(2:end), 'UniformOutput', false);
                else
                    newNamesAlongAxis = cellfun(@(v1, v2) [v1 ' - ' v2], stringList(1:end-1), ...
                        stringList(2:end), 'UniformOutput', false);
                end
                
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                stringList = stringLists{aIdx};
                if ~reverse
                    newNamesShortAlongAxis = cellfun(@(v1, v2) [v2 ' - ' v1], stringList(1:end-1), ...
                        stringList(2:end), 'UniformOutput', false);
                else
                    newNamesShortAlongAxis = cellfun(@(v1, v2) [v1 ' - ' v2], stringList(1:end-1), ...
                        stringList(2:end), 'UniformOutput', false);
                end
                
            else
                if isempty(p.Results.newNamesAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                    if ~reverse
                        newNamesAlongAxis = valueLists{aIdx}(2:end);
                    else
                        newNamesAlongAxis = valueLists{aIdx}(1:end-1);
                    end   
                else
                    newNamesAlongAxis = p.Results.newNamesAlongAxis;
                end
                
                if isempty(p.Results.newNamesShortAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                    if ~reverse
                        newNamesShortAlongAxis = valueLists{aIdx}(2:end);
                    else
                        newNamesShortAlongAxis = valueLists{aIdx}(1:end-1);
                    end  
                else
                    newNamesShortAlongAxis = p.Results.newNamesAlongAxis;
                end
            end
            
            newValueList = pset.conditionDescriptor.getAxisValueList(aIdx);
            if ~reverse
                newValueList = newValueList(2:end);
            else
                newValueList = newValueList(1:end-1);
            end
            
            % generate differencing matrix
            if ~reverse
                % compute conditions 2-1, 3-2, 4-3, etc.
                mat = diag(-onesvec(nAlongAxis)) + diag(onesvec(nAlongAxis-1), 1);
                wNbyO = mat(1:nAlongAxis-1, :);
            else
                % compute conditions 1-2, 2-3, 3-4, etc.
                mat = diag(onesvec(nAlongAxis)) + diag(-onesvec(nAlongAxis-1), 1);
                wNbyO = mat(1:nAlongAxis-1, :);
            end
            
            psetDiff = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newValueListAlongAxis', newValueList, ...
                'newNamesAlongAxis', newNamesAlongAxis, 'newNamesShortAlongAxis', newNamesShortAlongAxis, p.Unmatched);
        end
        
        function psetDiff = subtractOneConditionFromOthersAlongAxis(pset, axisName, conditionToSubtract, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('autoNamesAlongAxis', true, @islogical);
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.addParameter('newNamesShortAlongAxis', {}, @iscellstr);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            nAlongAxis = pset.conditionsSize(aIdx);
            
            idxKeep = setdiff(1:nAlongAxis, conditionToSubtract);
            
            % generate new names from differences
            if p.Results.autoNamesAlongAxis
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                stringList = stringLists{aIdx}(idxKeep);
                stringSubtract = stringLists{aIdx}{conditionToSubtract};
                newNamesAlongAxis = cellfun(@(v1) [v1 ' - ' stringSubtract], stringList, 'UniformOutput', false);
                
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                stringSubtract = stringLists{aIdx}{conditionToSubtract};
                newNamesShortAlongAxis = cellfun(@(v1) [v1 ' - ' stringSubtract], stringList, 'UniformOutput', false);
                
            else
                if isempty(p.Results.newNamesAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                    newNamesAlongAxis = valueLists{aIdx}(idxKeep);
                else
                    newNamesAlongAxis = p.Results.newNamesAlongAxis;
                end
                
                if isempty(p.Results.newNamesShortAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                    newNamesShortAlongAxis = valueLists{aIdx}(idxKeep);
                else
                    newNamesShortAlongAxis = p.Results.newNamesAlongAxis;
                end
            end
           
            % generate new value lists
            newValueList = pset.conditionDescriptor.getAxisValueList(aIdx);
            newValueList = newValueList(idxKeep);
            
            % generate differencing matrix
            mat = eye(nAlongAxis);
            mat(:, conditionToSubtract) = -1;
            wNbyO = mat(idxKeep, :);
            
            psetDiff = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newValueListAlongAxis', newValueList, ...
                'newNamesAlongAxis', newNamesAlongAxis, 'newNamesShortAlongAxis', newNamesShortAlongAxis, p.Unmatched);
        end
        
        function psetMean = computeMeanAlongAxis(pset, axisName, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('autoNamesAlongAxis', true, @islogical);
            p.addParameter('newNamesAlongAxis', {}, @iscell);
            p.addParameter('newNamesShortAlongAxis', {}, @iscellstr);
            p.addParameter('removeAxis', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            nAlongAxis = pset.conditionsSize(aIdx);
            axisName = pset.conditionDescriptor.axisNames{aIdx};
            
            % generate new names from differences
            if p.Results.autoNamesAlongAxis
                newNamesAlongAxis = { sprintf('Mean Over %s', axisName) }; 
                newNamesShortAlongAxis = { sprintf('Mean Over %s', axisName) }; 
            else
                newNamesAlongAxis = p.Results.newNamesAlongAxis;
                newNamesShortAlongAxis = p.Results.newNamesShortAlongAxis;
            end
            
            newValueList = pset.conditionDescriptor.getAxisValueList(aIdx);
            newValueList = newValueList(1);
            flds = fieldnames(newValueList);
            for iF = 1:numel(flds)
                newValueList.(flds{iF}) = sprintf('Mean Over %s', axisName);
            end
            
            % normalization is done by normalizeCoefficientsByNumNonNaN
            % below
            wNbyO = ones(1, nAlongAxis);
            
            psetMean = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                aIdx, wNbyO, 'removeAxis', p.Results.removeAxis, 'newValueListAlongAxis', newValueList, ...
                'newNamesAlongAxis', newNamesAlongAxis, 'newNamesShortAlongAxis', newNamesShortAlongAxis, ...
                'replaceNaNWithZero', true, 'normalizeCoefficientsByNumNonNaN', true, p.Unmatched);
        end
        
        function psetMean = subtractMeanAlongAxis(pset, axisName, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('autoNamesAlongAxis', true, @islogical);
            p.addParameter('newNamesAlongAxis', '', @iscell);
            p.addParameter('newNamesShortAlongAxis', {}, @iscellstr);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            axisName = pset.conditionDescriptor.axisNames{aIdx};
            nAlongAxis = pset.conditionsSize(aIdx);
            
            % generate new names from differences
            if p.Results.autoNamesAlongAxis
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                stringList = stringLists{aIdx};
                newNamesAlongAxis = cellfun(@(s) sprintf('%s mean-subtracted across %s', s, axisName), stringList, 'UniformOutput', false); 
                
                stringLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                stringList = stringLists{aIdx};
                newNamesShortAlongAxis = cellfun(@(s) sprintf('%s mean-subtracted across %s', s, axisName), stringList, 'UniformOutput', false); 
                
            else
                if isempty(p.Results.newNamesAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', false);
                    newNamesAlongAxis = valueLists{aIdx};
                else
                    newNamesAlongAxis = p.Results.newNamesAlongAxis;
                end
                
                if isempty(p.Results.newNamesShortAlongAxis)
                    % keep same names
                    valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings('short', true);
                    newNamesShortAlongAxis = valueLists{aIdx};
                else
                    newNamesShortAlongAxis = p.Results.newNamesAlongAxis;
                end
            end
            
            newValueList = pset.conditionDescriptor.getAxisValueList(aIdx);
            
            % normalization is done by 'normalizeCoefficientsByNumNonNaN' below
            % and the identity matrix is added in after normalization by
            % 'addToOriginal'
            wNbyO = -ones(nAlongAxis, nAlongAxis);
            
            psetMean = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newValueListAlongAxis', newValueList, ...
                'newNamesAlongAxis', newNamesAlongAxis, 'newNamesShortAlongAxis', newNamesShortAlongAxis, ...
                'replaceNaNWithZero', true, ...
                'normalizeCoefficientsByNumNonNaN', true, ...
                'addToOriginal', true, ...
                p.Unmatched);
        end
        
        function psetReweighted = applyLinearCombinationAlongConditionAxis(pset, axisName, weightsNewCByOldC, varargin)
            p = inputParser();
            p.addParameter('newValueListAlongAxis', [], @isstruct);
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.addParameter('newNamesShortAlongAxis', {}, @iscellstr);
            p.addParameter('removeAxis', false, @islogical);
            p.addParameter('conditionAppearanceFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            
            % if true, ignore NaNs by replacing them with zero. by enabling
            % this flag, you allow combined data to be valid if _any_ of
            % the conditions that contribute to the combination are valid.
            % If false, all conditions that contribute must be valid for
            % the combined data to be valid
            p.addParameter('replaceNaNWithZero', false, @islogical);
            
            % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis
            % this enables nanmean like computations
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical); 
            
            % requires that the matrix be square, the equivalent of adding
            % the identity matrix to the weight matrix, except that this
            % will be added after normalization
            p.addParameter('addToOriginal', false, @islogical);
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            wNbyO = weightsNewCByOldC;
            cOld = pset.conditionsSize(aIdx);
            assert(size(wNbyO, 2) == cOld, ...
                'Weighting matrix must have same column count as number of existing conditions along axis (%d)', cOld);
            
            cNewAlongAxis = size(wNbyO, 1);
            
            removeAxis = p.Results.removeAxis;
            
            if ~removeAxis
                newNamesAlongAxis = p.Results.newNamesAlongAxis;
                if isempty(newNamesAlongAxis)
                    newNamesAlongAxis = arrayfun(@(i) sprintf('Condition Combination %d', i), 1:cNewAlongAxis, 'UniformOutput', false);
                end
                assert(numel(newNamesAlongAxis) == cNewAlongAxis, 'newNamesAlongAxis must have numel == number of new conditions (%d)', cNewAlongAxis);

                newNamesShortAlongAxis = p.Results.newNamesShortAlongAxis;
                if ~isempty(newNamesShortAlongAxis)
                    assert(numel(newNamesShortAlongAxis) == cNewAlongAxis, 'newNamesShortAlongAxis must have numel == number of new conditions (%d)', cNewAlongAxis);
                end
                
                newValueListAlongAxis = p.Results.newValueListAlongAxis;
                assert(numel(newValueListAlongAxis) == cNewAlongAxis, 'newNamesAlongAxis must have numel == number of new conditions (%d)', cNewAlongAxis);
                
            end
            
            % pset.warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns();
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % setup new condition descriptor, optionally drop the axis
            % we're combining along, if the new size is 1
            if removeAxis
                assert(cNewAlongAxis == 1, 'New condition count along axis must be 1 in order to removeAxis');
                newCD = pset.conditionDescriptor.removeAxis(aIdx);
            else
                newCD = pset.conditionDescriptor.setAxisValueList(axisName, newValueListAlongAxis, ...
                    'asStrings', newNamesAlongAxis, 'asStringsShort', newNamesShortAlongAxis);
            end
            
            % update the condition appearance fn if specified
            if ~ismember('conditionAppearanceFn', p.UsingDefaults) % we dont just check isempty b/c the user may or may not want to set it to empty
                newCD.appearanceFn = p.Results.conditionAppearanceFn;
            end
            
            b.conditionDescriptor = newCD;
            
            nConditionsNew = b.conditionDescriptor.nConditions;
%             conditionsSizeNew = b.conditionDescriptor.conditionsSize;
            
            % adjust mean and sem to reflect difference between conditions
            [b.dataMean, b.dataSem] = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                % build N x TA x C1 x C2 x ...
                tensorMean = pset.arrangeNbyTAbyConditionAttributes('alignIdx', iAlign);
                tensorMeanReweighted = TensorUtils.linearCombinationAlongDimension(tensorMean, aIdx+2, wNbyO, ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal);
                % back to N x C x TA
                b.dataMean{iAlign} = permute(tensorMeanReweighted(:, :, :), [1 3 2]);
                
                % build N x TA x C1 x C2 x ...
                % use sd1+2 = sqrt(sd1^2 / n1 + sd2^2 / n2) formula
                % which here means semNew = sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                tensorSem = pset.arrangeNbyTAbyConditionAttributes('type', 'sem', 'alignIdx', iAlign);
                tensorSemReweighted = sqrt( TensorUtils.linearCombinationAlongDimension(tensorSem.^2, aIdx+2, abs(wNbyO), ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal) );
                % back to N x C x TA
                b.dataSem{iAlign} = permute(tensorSemReweighted(:, :, :), [1 3 2]);
            end
            
            % reshape pset.dataCachedSampledTrialsTensor
            % N x TA x C x Trials -> N x TA x size(conditions) x Trials
            if ~isempty(pset.dataCachedSampledTrialsTensor)
                cachedTrialsAttr = reshape(pset.dataCachedSampledTrialsTensor, ...
                    [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize), size(pset.dataCachedSampledTrialsTensor, 4)]);
                newTrialsAttr = TensorUtils.linearCombinationAlongDimension(cachedTrialsAttr, ...
                    aIdx+2, wNbyO, ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal);
                b.dataCachedSampledTrialsTensor = reshape(newTrialsAttr, ...
                    [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew, size(pset.dataCachedSampledTrialsTensor, 4)]);
                
                cachedTrialsAttr = reshape(pset.dataCachedMeanExcludingSampledTrialsTensor, ...
                    [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize), size(pset.dataCachedMeanExcludingSampledTrialsTensor, 4)]);
                newTrialsAttr = TensorUtils.linearCombinationAlongDimension(cachedTrialsAttr, ...
                    aIdx+2, wNbyO, ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal);
                b.dataCachedMeanExcludingSampledTrialsTensor = reshape(newTrialsAttr, ...
                    [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew, size(pset.dataCachedSampledTrialsTensor, 4)]);
                
                % N x C: need to compute min over number of trials for conditions included
                % in that new combined condition
                trialCountsNbyAttr = reshape(pset.dataCachedSampledTrialCounts, [pset.nBases, makerow(pset.conditionsSize)]);
                trialCountsNew = TensorUtils.linearCombinationApplyScalarFnAlongDimension(trialCountsNbyAttr, aIdx+1, wNbyO, @min); 
                b.dataCachedSampledTrialCounts = reshape(trialCountsNew, [pset.nBases, nConditionsNew]);
            end
            
            % update difference of trials scaled noise estimates so that we
            % can compute noise variance floors when projecting. since the
            % noise estimates are already scaled by 1/sqrt(2*nTrials), we
            % simply add them together to get the new scaled estimate
            if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimate)
                scaledNoiseEstimate_NbyTAbyC = pset.dataDifferenceOfTrialsScaledNoiseEstimate;
                scaledNoiseEstimate_NbyTAbyAttr = reshape(scaledNoiseEstimate_NbyTAbyC, [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize)]);
                newScaledNoiseEstimate_NbyTAbyAttr = TensorUtils.linearCombinationAlongDimension(...
                    scaledNoiseEstimate_NbyTAbyAttr, aIdx+2, abs(wNbyO), ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal);
                b.dataDifferenceOfTrialsScaledNoiseEstimate = reshape(newScaledNoiseEstimate_NbyTAbyAttr, ...
                    [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew]);
            end
            
            % diff randomized data if present, recompute intervals
            if ~isempty(pset.dataMeanRandomized)
                % setup new condition descriptor, optionally drop the axis
                % we're combining along, if the new size is 1
                newCD = pset.conditionDescriptorRandomized.setAxisValueList(axisName, newNamesAlongAxis);
                if p.Results.removeAxis
                    assert(cNewAlongAxis == 1, 'New condition count along axis must be 1 in order to removeAxis');
                    newCD = newCD.removeAxis(aIdx);
                end
                % update the condition appearance fn if specified
                if ~ismember('conditionAppearanceFn', p.UsingDefaults) % we dont just check isempty b/c the user may or may not want to set it to empty
                    newCD.appearanceFn = p.Results.conditionAppearanceFn;
                end
                b.conditionDescriptorRandomized = newCD;
                
                [b.dataMeanRandomized, b.dataSemRandomized] = deal(cell(pset.nAlign, 1));
                for iAlign = 1:pset.nAlign
                    % dataMeanRandomized is N x C x TA x R (where R is number of random samples)
                    % randTensor is N x C1 x C2 x ... x TA x R
                    meanTensor = reshape(pset.dataMeanRandomized{iAlign}, ...
                        [pset.nBases, makerow(pset.conditionsSize), pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    
                    meanTensorReweighted = TensorUtils.linearCombinationAlongDimension(meanTensor, aIdx+1, wNbyO, ...
                        'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                        'keepNaNIfAllNaNs', true, ...
                        'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                        'addToOriginal', p.Results.addToOriginal); 
                    
                    % back to N x C x TA x R
                    b.dataMeanRandomized{iAlign} = reshape(meanTensorReweighted, [pset.nBases, nConditionsNew, pset.nTimeDataMean(iAlign), pset.nRandomSamples]);

                    % ensure invalid bases remain invalid
                    b.dataMeanRandomized{iAlign}(~pset.basisValid, :, :, :) = NaN;
                    
                    % dataSemRandomized is N x C x TA x R (where R is number of random samples)
                    % randTensor is N x C1 x C2 x ... x TA x R
                    semTensor = reshape(pset.dataSemRandomized{iAlign}, ...
                        [pset.nBases, makerow(pset.conditionsSize), pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    
                    semTensorReweighted = TensorUtils.linearCombinationAlongDimension(semTensor, aIdx+1, wNbyO, ...
                        'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                        'keepNaNIfAllNaNs', true, ...
                        'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                        'addToOriginal', p.Results.addToOriginal); 

                    % back to N x C x TA x R
                    b.dataSemRandomized{iAlign} = reshape(semTensorReweighted, [pset.nBases, nConditionsNew, pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    b.dataSemRandomized{iAlign}(~pset.basisValid, :, :, :) = NaN;
                end
                
                % and scale randomized difference of trials
                scaledNoiseEstimate_NbyTAbyCbyS = pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized;
                scaledNoiseEstimate_NbyTAbyAttrbyS = reshape(scaledNoiseEstimate_NbyTAbyCbyS, ...
                    [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize), pset.nRandomSamples]);
                newScaledNoiseEstimate_NbyTAbyAttrbyS = TensorUtils.linearCombinationAlongDimension(...
                    scaledNoiseEstimate_NbyTAbyAttrbyS, aIdx+2, abs(wNbyO), ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', p.Results.normalizeCoefficientsByNumNonNaN, ...
                    'addToOriginal', p.Results.addToOriginal);
                b.dataDifferenceOfTrialsScaledNoiseEstimateRandomized = reshape(newScaledNoiseEstimate_NbyTAbyAttrbyS, ...
                [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew, pset.nRandomSamples]);
            end
            
            b.trialLists = {}; % no longer relevant
            
            % A x N x C
            % for data valid, we need all input conditions to be valid for
            % output conditions to be valid, so we change wNbyO such that
            % the linear combination will be 1 iff all/any bases that
            % contribute to that output are valid. All is if
            % replaceNanWithZero is false, any is if replaceNanWithZero is true
            wNbyO_forValid = bsxfun(@rdivide, wNbyO ~= 0, sum(wNbyO ~= 0, 2));
            [dataValidTensor, cdims] = TensorUtils.reshapeDimsInPlace(pset.dataValid, 3, pset.conditionsSize);
            
            if p.Results.replaceNaNWithZero
                b.dataValid = TensorUtils.flattenDimsInPlace(TensorUtils.linearCombinationAlongDimension(...
                    dataValidTensor, aIdx+2, wNbyO_forValid) == 1, cdims);
            else
                b.dataValid = TensorUtils.flattenDimsInPlace(TensorUtils.linearCombinationAlongDimension(...
                    dataValidTensor, aIdx+2, wNbyO_forValid) ~= 0, cdims);
            end
            
            % sum trials from all included conditions
            [dataNTrialsTensor, cdims] = TensorUtils.reshapeDimsInPlace(pset.dataNTrials, 3, pset.conditionsSize);
            b.dataNTrials = TensorUtils.flattenDimsInPlace(TensorUtils.linearCombinationAlongDimension(...
                dataNTrialsTensor, aIdx+2, wNbyO ~= 0, 'replaceNaNWithZero', p.Results.replaceNaNWithZero), cdims);

            % shrink the time windows over all considered conditions
            if p.Results.replaceNaNWithZero
                % don't worry about nans in the time windows since we're
                % taking data from any condition along the axis
                nanMode = 'omitnan';
            else
                % we must use includenan here otherwise we'll end up having
                % tMinValid where no data is present
                nanMode = 'includenan';
            end
            
            [tMinValidOld, cdims] = TensorUtils.reshapeDimsInPlace(pset.tMinValidByAlignBasisCondition, 3, pset.conditionsSize);
            
            tMinValidCellByNew = arrayfun(@(iNew) max(TensorUtils.selectAlongDimension(tMinValidOld, aIdx+2, wNbyO(iNew, :) ~= 0), [], aIdx+2, nanMode), ...
                1:cNewAlongAxis, 'UniformOutput', false);
            tMinValidNew = cat(aIdx+2, tMinValidCellByNew{:});
            b.tMinValidByAlignBasisCondition = TensorUtils.flattenDimsInPlace(tMinValidNew, cdims);
            
            [tMaxValidOld, cdims] = TensorUtils.reshapeDimsInPlace(pset.tMaxValidByAlignBasisCondition, 3, pset.conditionsSize);
            tMaxValidCellByNew = arrayfun(@(iNew) min(TensorUtils.selectAlongDimension(tMaxValidOld, aIdx+2, wNbyO(iNew, :) ~= 0), [], aIdx+2, nanMode), ...
                1:cNewAlongAxis, 'UniformOutput', false);
            tMaxValidNew = cat(aIdx+2, tMaxValidCellByNew{:});
            b.tMaxValidByAlignBasisCondition = TensorUtils.flattenDimsInPlace(tMaxValidNew, cdims);
            
            cIndsTensor = TensorUtils.containingLinearInds(pset.conditionDescriptor.conditionsSize);
            setsAlongAxis = arrayfun(@(iNew) find(wNbyO(iNew, :)), 1:cNewAlongAxis, 'UniformOutput', false);
            conditionIdxSetsTensor = TensorUtils.selectSetsAlongDimension(cIndsTensor, aIdx, setsAlongAxis);
            
            prog = ProgressBar(pset.nAlignSummaryData, 'Collecting conditions within AlignSummary data');
            for iSource = 1:pset.nAlignSummaryData
                prog.update(iSource);
                for iAlign = 1:pset.nAlign
                    if ~isempty(b.alignSummaryData{iSource, iAlign})
                        b.alignSummaryData{iSource, iAlign} = b.alignSummaryData{iSource, iAlign}.combineSetsOfConditions(...
                            b.conditionDescriptor, conditionIdxSetsTensor(:));
                    end
                end
            end
            prog.finish();
            
            psetReweighted = b.buildManualWithTrialAveragedData();
        end
    end
    
    % concatenating along a condition axis
    methods(Static)
        function psetCat = concatenateAlongNewConditionAxis(psetCell, axisName, axisValueList, varargin)
            p = inputParser();
            p.addParameter('aggregateMarks', true, @islogical);
            p.addParameter('aggregateIntervals', true, @islogical);
            p.addParameter('conditionAppearanceFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('equalizeTimeVectors', false, @islogical); % takes time but necessary if time vectors (min to max) are not matched exactly 
            p.parse(varargin{:});
            
            % check for equal conditions size
            sz = psetCell{1}.conditionsSize;
            for i = 2:numel(psetCell)
                assert(isequal(sz, psetCell{i}.conditionsSize), 'Condition Size for psetCell{%d} does not match psetCell{1}', i);
            end
                
            % expand the time vectors of each pset to the full common
            tMinByAlign = min(cell2mat(cellfun(@(p) p.tMinForDataMean, makerow(psetCell), 'UniformOutput', false)), [], 2);
            tMaxByAlign = max(cell2mat(cellfun(@(p) p.tMaxForDataMean, makerow(psetCell), 'UniformOutput', false)), [], 2);
            
            if p.Results.equalizeTimeVectors
                prog = ProgressBar(numel(psetCell), 'Equalizing time vectors across psets');
                for iP = 1:numel(psetCell)
                    prog.update(iP);
                    psetCell{iP} = psetCell{iP}.manualSliceOrExpandTimeWindow(tMinByAlign, tMaxByAlign);
                end
                prog.finish();
            end
            
%             PopulationTrajectorySet.assertBasesMatch(psetCell{:});
%             assert(PopulationTrajectorySet.checkSameBasesValid(psetCell{:}), ...
%                 'Psets being concatenated do not have the same .basisValid. Use PopulationTrajectorySet.equalizeBasesValid first');
            pset = psetCell{1};
            
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % take union of valid bases and make this the new basis valid
            basisValid = psetCell{1}.basisValid;
            for i = 2:numel(psetCell)
                basisValid = basisValid | psetCell{i}.basisValid;
            end
            b.basisValidManual = basisValid;
            if isempty(b.basisInvalidCauseManual)
                b.basisInvalidCauseManual = cell(numel(basisValid), 1);
                b.basisInvalidCauseManual(:) = {''};
            else
                b.basisInvalidCauseManual(~basisValid) = {''};
            end
            
            assert(numel(axisValueList) == numel(psetCell), 'Value list along new concatenation axis must have same length as psetCell');
            
            % update condition descriptor
            cd = pset.conditionDescriptor;
            cd = cd.addAttribute(axisName, 'valueList', axisValueList);
            cd = cd.addAxis(axisName, 'valueList', axisValueList);
            % update the condition appearance fn if specified
            if ~ismember('conditionAppearanceFn', p.UsingDefaults) % we dont just check isempty b/c the user may or may not want to set it to empty
                cd.appearanceFn = p.Results.conditionAppearanceFn;
            end
            b.conditionDescriptor = cd;
            
            cAxis = cd.nAxes;
            
            cSize = pset.conditionsSizeNoExpand;
            cSize(end+1) = 1; % make room for new axis
            
            if ~p.Results.aggregateMarks
                b.alignDescriptorSet = cellfun(@(ad) ad.clearMarks(), b.alignDescriptorSet, 'UniformOutput', false);
            end
            if ~p.Results.aggregateIntervals
                b.alignDescriptorSet = cellfun(@(ad) ad.clearIntervals(), b.alignDescriptorSet, 'UniformOutput', false);
            end
            
            if ~isempty(pset.conditionDescriptorRandomized)
                cdRand = pset.conditionDescriptorRandomized;
                cdRand = cdRand.addAttribute(axisName, 'valueList', axisValueList);
                cdRand = cdRand.addAxis(axisName, 'valueList', axisValueList);
                b.conditionDescriptorRandomized = cdRand;
            end
            
            % pre touch all the data means to ensure that buildDataMean
            % has run on each pset. this is time consuming and best done up
            % front so we get a nice progress bar
            prog = ProgressBar(numel(psetCell), 'Generating trial-averaged data');
            for i = 1:numel(psetCell)
                prog.update(i);
                temp = psetCell{i}.dataMean; %#ok<NASGU>
                temp = psetCell{i}.dataDifferenceOfTrialsScaledNoiseEstimate; %#ok<NASGU>
                temp = psetCell{i}.dataCachedSampledTrialsTensor; %#ok<NASGU>
            end
            prog.finish();

            debug('Concatenating trial-averaged data\n');
            for iAlign = 1:pset.nAlign
                % N x C x T
                b.dataMean{iAlign} = catConditionsFlat(psetCell, @(p) p.dataMean{iAlign}, 2, cSize, cAxis);
                b.dataSem{iAlign} = catConditionsFlat(psetCell, @(p) p.dataSem{iAlign}, 2, cSize, cAxis);
            end
            
            % A x N x C
            b.tMinValidByAlignBasisCondition = catConditionsFlat(psetCell, @(p) p.tMinValidByAlignBasisCondition, 3, cSize, cAxis);
            b.tMaxValidByAlignBasisCondition = catConditionsFlat(psetCell, @(p) p.tMaxValidByAlignBasisCondition, 3, cSize, cAxis);
            b.dataNTrials = catConditionsFlat(psetCell, @(p) p.dataNTrials, 3, cSize, cAxis);
            b.dataValid = catConditionsFlat(psetCell, @(p) p.dataValid, 3, cSize, cAxis);
            
            % N x C
            b.trialLists = catConditionsFlat(psetCell, @(p) p.trialLists, 2, cSize, cAxis);
            
            % N x TA x C
            b.dataDifferenceOfTrialsScaledNoiseEstimate = catConditionsFlat(psetCell, @(p) p.dataDifferenceOfTrialsScaledNoiseEstimate, 3, cSize, cAxis);
        
            % N x TA x C x Trials
            b.dataCachedSampledTrialsTensor = catConditionsFlat(psetCell, @(p) p.dataCachedSampledTrialsTensor, 3, cSize, cAxis);
            b.dataCachedMeanExcludingSampledTrialsTensor = catConditionsFlat(psetCell, @(p) p.dataCachedMeanExcludingSampledTrialsTensor, 3, cSize, cAxis);
            
            % N x C 
            b.dataCachedSampledTrialCounts = catConditionsFlat(psetCell, @(p) p.dataCachedSampledTrialCounts, 2, cSize, cAxis);

            % adjust alignSummaryData
            % N x A
            temp = pset.alignSummaryData; %#ok<NASGU> % request up front to trigger computation before progress bar
            prog = ProgressBar(pset.nBases, 'Aggregating AlignSummary data');
            for iSource = 1:pset.nAlignSummaryData
                prog.update(iSource);
                for iAlign = 1:pset.nAlign
                    emptyMask = cellfun(@(pset) isempty(pset.alignSummaryData{iSource, iAlign}), psetCell);
                    if all(emptyMask) || ~pset.basisValid(iSource)
                        % missing entirely, that's okay
                        b.alignSummaryData{iSource, iAlign} = [];
                        
                    elseif ~any(emptyMask)
                        % all found, combine
                        clear alignSummarySet;
                        for iP = numel(psetCell):-1:1
                            alignSummarySet(iP) = psetCell{iP}.alignSummaryData{iSource, iAlign};
                        end
                        b.alignSummaryData{iSource, iAlign} = ...
                            AlignSummary.aggregateByConcatenatingConditionsAlongNewAxis(alignSummarySet, cd, ...
                            'aggregateMarks', p.Results.aggregateMarks, 'aggregateIntervals', p.Results.aggregateIntervals);
                    else
                        % missing for some, that's a problem
                        error('AlignSummary missing for %d PopTrajSets for source %d align %d', nnz(emptyMask), iSource, iAlign);
                    end
                    
                end
            end
            prog.finish();

            hasDataRandomized = cellfun(@(p) p.hasDataRandomized, psetCell);
            if all(hasDataRandomized)
                % update condition descriptor
                cd = pset.conditionDescriptorRandomized;
                cd = cd.addAttribute(axisName, 'valueList', axisValueList);
                cd = cd.addAxis(axisName, 'valueList', axisValueList);
                % update the condition appearance fn if specified
                if ~ismember('conditionAppearanceFn', p.UsingDefaults) % we dont just check isempty b/c the user may or may not want to set it to empty
                    cd.appearanceFn = p.Results.conditionAppearanceFn;
                end
                b.conditionDescriptorRandomized = cd;
                
                debug('Concatenating data randomized\n');
                for iAlign = 1:pset.nAlign
                    b.dataMeanRandomized{iAlign} = catConditionsFlat(psetCell, @(p) p.dataMeanRandomized{iAlign}, 2, cSize, cAxis);
                    b.dataSemRandomized{iAlign} = catConditionsFlat(psetCell, @(p) p.dataSemRandomized{iAlign}, 2, cSize, cAxis);
                end
                
                % N x T x C by S
                b.dataDifferenceOfTrialsScaledNoiseEstimateRandomized = catConditionsFlat(psetCell, @(p) p.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, 3, cSize, cAxis);
            end
            
            psetCat = b.buildManualWithTrialAveragedData();

            function res = catConditionsFlat(objCell, accessFn, conditionDim, conditionsSize, conditionAxis)
%                 if ~exist('timeDim', 'var')
%                     timeDim = [];
%                     timeCell = {};
%                 else
%                     timeCell = cellfun(timeAccessFn, objCell, 'UniformOutput', false);
%                 end
                
                % grab from each obj in cell
                tensorCell = cellfun(accessFn, objCell, 'UniformOutput', false);
                
                % reshape to make conditions a tensor rather than flat
                sz = TensorUtils.expandSizeToNDims(size(tensorCell{1}), conditionDim);
                newSz = [sz(1:conditionDim-1), conditionsSize, sz(conditionDim+1:end)];
                tz = cellfun(@(x) reshape(x, newSz), tensorCell, 'UniformOutput', false);
                tensorCell = tz;
                % concatenate along the condition axis
                catAxis = conditionDim + conditionAxis-1;
                
%                 if isempty(timeDim)
                    catTensor = cat(catAxis, tensorCell{:});
%                 else
%                     [catTensor, tvec] = TensorUtils.concatenateAlignedToCommonTimeVector(timeCell, tensorCell, timeDim, catAxis);
%                     sz(timeDim) = numel(tvec);
%                 end
                
                % reshape the result to be flat again
                catSz = sz;
                catSz(conditionDim) = catSz(conditionDim) * numel(objCell);
                res = reshape(catTensor, catSz);
            end
            
        end
        
        function psetCat = concatenateAlongExistingConditionAxis(psetCell, axisName, varargin)
            p = inputParser();
            p.addParameter('aggregateMarks', true, @islogical);
            p.addParameter('aggregateIntervals', true, @islogical);
            p.addParameter('conditionAppearanceFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('equalizeTimeVectors', false, @islogical); % takes time but necessary if time vectors (min to max) are not matched exactly 
            p.parse(varargin{:});
            
            if p.Results.equalizeTimeVectors
                psetCell = PopulationTrajectorySetCrossConditionUtilities.equalizeTimeVectors(psetCell);
            end
            pset = psetCell{1};
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            % check for equal conditions size
            sz = psetCell{1}.conditionsSize;
            sz(aIdx) = 1;
            for i = 2:numel(psetCell)
                szThis = psetCell{i}.conditionsSize;
                szThis(aIdx) = 1;
                assert(isequal(sz, szThis), 'Condition Size for psetCell{%d} does not match psetCell{1}', i);
            end

            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % take union of valid bases and make this the new basis valid
            basisValid = psetCell{1}.basisValid;
            for i = 2:numel(psetCell)
                basisValid = basisValid | psetCell{i}.basisValid;
            end
            b.basisValidManual = basisValid;
            if isempty(b.basisInvalidCauseManual)
                b.basisInvalidCauseManual = cell(numel(basisValid), 1);
                b.basisInvalidCauseManual(:) = {''};
            else
                b.basisInvalidCauseManual(~basisValid) = {''};
            end
            
            % update condition descriptor
            cdCell = cellfun(@(pset) pset.conditionDescriptor, psetCell, 'UniformOutput', false);
            cd = ConditionDescriptor.concatenateAlongAxis(cdCell, aIdx);
            b.conditionDescriptor = cd;
            
            cAxis = aIdx;
            cSize = cellfun(@(pset) pset.conditionsSizeNoExpand, psetCell, 'UniformOutput', false);
            
            if ~p.Results.aggregateMarks
                b.alignDescriptorSet = cellfun(@(ad) ad.clearMarks(), b.alignDescriptorSet, 'UniformOutput', false);
            end
            if ~p.Results.aggregateIntervals
                b.alignDescriptorSet = cellfun(@(ad) ad.clearIntervals(), b.alignDescriptorSet, 'UniformOutput', false);
            end
            
            % pre touch all the data means to ensure that buildDataMean
            % has run on each pset. this is time consuming and best done up
            % front so we get a nice progress bar
            prog = ProgressBar(numel(psetCell), 'Generating trial-averaged data');
            for i = 1:numel(psetCell)
                prog.update(i);
                temp = psetCell{i}.dataMean; %#ok<NASGU>
                temp = psetCell{i}.dataDifferenceOfTrialsScaledNoiseEstimate; %#ok<NASGU>
                temp = psetCell{i}.dataCachedSampledTrialsTensor; %#ok<NASGU>
            end
            prog.finish();

            debug('Concatenating trial-averaged data\n');
            for iAlign = 1:pset.nAlign
                % N x C x T
                b.dataMean{iAlign} = catConditions(psetCell, @(p) p.dataMean{iAlign}, 2, cSize, cAxis);
                b.dataSem{iAlign} = catConditions(psetCell, @(p) p.dataSem{iAlign}, 2, cSize, cAxis);
            end
            
            % A x N x C
            b.tMinValidByAlignBasisCondition = catConditions(psetCell, @(p) p.tMinValidByAlignBasisCondition, 3, cSize, cAxis);
            b.tMaxValidByAlignBasisCondition = catConditions(psetCell, @(p) p.tMaxValidByAlignBasisCondition, 3, cSize, cAxis);
            b.dataNTrials = catConditions(psetCell, @(p) p.dataNTrials, 3, cSize, cAxis);
            b.dataValid = catConditions(psetCell, @(p) p.dataValid, 3, cSize, cAxis);
            
            % N x C
            b.trialLists = catConditions(psetCell, @(p) p.trialLists, 2, cSize, cAxis, true);
            
            % N x TA x C
            b.dataDifferenceOfTrialsScaledNoiseEstimate = catConditions(psetCell, @(p) p.dataDifferenceOfTrialsScaledNoiseEstimate, 3, cSize, cAxis, true);
        
            % N x TA x C x Trials
            b.dataCachedSampledTrialsTensor = catConditions(psetCell, @(p) p.dataCachedSampledTrialsTensor, 3, cSize, cAxis, true);
            b.dataCachedMeanExcludingSampledTrialsTensor = catConditions(psetCell, @(p) p.dataCachedMeanExcludingSampledTrialsTensor, 3, cSize, cAxis, true);
            
            % N x C 
            b.dataCachedSampledTrialCounts = catConditions(psetCell, @(p) p.dataCachedSampledTrialCounts, 2, cSize, cAxis, true);

            % adjust alignSummaryData
            % N x A
            temp = pset.alignSummaryData; %#ok<NASGU> % request up front to trigger computation before progress bar
            prog = ProgressBar(pset.nBases, 'Aggregating AlignSummary data');
            for iSource = 1:pset.nAlignSummaryData
                prog.update(iSource);
                for iAlign = 1:pset.nAlign
                    emptyMask = cellfun(@(pset) isempty(pset.alignSummaryData{iSource, iAlign}), psetCell);
                    if all(emptyMask) || ~pset.basisValid(iSource)
                        % missing entirely, that's okay
                        b.alignSummaryData{iSource, iAlign} = [];
                        
                    elseif ~any(emptyMask)
                        % all found, combine
                        clear alignSummarySet;
                        for iP = numel(psetCell):-1:1
                            alignSummarySet(iP) = psetCell{iP}.alignSummaryData{iSource, iAlign};
                        end
                        b.alignSummaryData{iSource, iAlign} = ...
                            AlignSummary.aggregateByConcatenatingConditionsAlongExistingAxis(alignSummarySet, cd, cAxis, ...
                            'aggregateMarks', p.Results.aggregateMarks, 'aggregateIntervals', p.Results.aggregateIntervals);
                    else
                        % missing for some, that's a problem
                        error('AlignSummary missing for %d PopTrajSets for source %d align %d', nnz(emptyMask), iSource, iAlign);
                    end
                    
                end
            end
            prog.finish();

           
            psetCat = b.buildManualWithTrialAveragedData();

            function res = catConditions(objCell, accessFn, conditionDim, conditionsSizeCell, catAlongAxis, emptyOkay)
                if nargin < 6
                    emptyOkay = false;
                end
                % grab from each obj in cell
                tensorCell = cellfun(accessFn, objCell, 'UniformOutput', false);
                
                if any(cellfun(@isempty, tensorCell))
                    if emptyOkay
                        res = {};
                        return;
                    else
                        error('Some values were empty');
                    end
                end
                
                % reshape to make conditions a tensor rather than flat
                for iO = 1:numel(tensorCell)
                    szt = TensorUtils.expandSizeToNDims(size(tensorCell{iO}), catAlongAxis);
                    newSz = [szt(1:conditionDim-1), conditionsSizeCell{iO}, szt(conditionDim+1:end)];
                    tensorCell{iO} = reshape(tensorCell{iO}, newSz);
                end
                
                % concatenate along the condition axis
                catAxis = conditionDim + catAlongAxis - 1;
                catTensor = cat(catAxis, tensorCell{:});
                
                % reshape the result to be flat again
                szMat = cat(1, conditionsSizeCell{:});
                catSz = newSz;
                catSz(catAxis) = sum(szMat(:, catAlongAxis));
                res = reshape(catTensor, catSz);
            end
        end
    end
   
    methods(Static)
       function psetSum = addToDataMean(pset, addToDataMean, varargin)
            p = inputParser();
            p.addParameter('gain1', 1, @isscalar);
            p.addParameter('gain2', 1, @isscalar);
            p.addParameter('sem', [], @(x) isempty(x) || iscell(x)); % will be treated as the S.E.M. of the quantity addToDataMean, and added to dataSem in the RMS sense
            p.addParameter('applyToSingleTrial', false, @islogical);
            p.parse(varargin{:});
            
            assert(~p.Results.applyToSingleTrial, 'Not yet supported');
            addSem = p.Results.sem;
            gain1 = p.Results.gain1;
            gain2 = p.Results.gain2;
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % adjust mean and sem to reflect difference between conditions
            [b.dataMean, b.dataSem] = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                b.dataMean{iAlign} = gain1 .* pset.dataMean{iAlign} + gain2 .* addToDataMean{iAlign};
                
                % build N x TA x C1 x C2 x ...
                % use sd1+2 = sqrt(sd1^2 / n1 + sd2^2 / n2) formula
                % which here means semNew = sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                if ~isempty(addSem)
                    b.dataSem{iAlign} = sqrt((gain1.*pset.dataSem{iAlign}).^2 + (gain2.*addSem{iAlign}).^2);
                end
            end
            
            b.trialLists = {}; % no longer relevant
            
            psetSum = b.buildManualWithTrialAveragedData();
       end 
        
       function psetCCM = crossConditionMean(pset, varargin)
            % similar in spirit to subtractMeanAlongAxis except subtracts
            % the global condition mean across all axes
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('conditionName', 'ccm', @ischar);
            p.parse(varargin{:});
           
            psetCCM = PopulationTrajectorySetCrossConditionUtilities.computeMeanAlongAxis(pset.flattenConditionAxes(), 1, 'removeAxis', true);
            if ~isempty(p.Results.conditionName)
                psetCCM = psetCCM.setConditionNames({p.Results.conditionName});
            end
       end
        
       function psetSum = addPset(pset, psetAdd, varargin)
           psetSum = PopulationTrajectorySetCrossConditionUtilities.addToDataMean(pset, psetAdd.dataMean, 'sem', psetAdd.dataSem, varargin{:});
       end
       
       function psetMinusCCM = subtractCrossConditionMean(pset, varargin)
           psetCCM = PopulationTrajectorySetCrossConditionUtilities.crossConditionMean(pset, varargin{:});
           psetMinusCCM = PopulationTrajectorySetCrossConditionUtilities.addPset(pset, psetCCM, 'gain2', -1);
       end
    end
    
    methods(Static) % Utilities
        function psetCell = equalizeTimeVectors(psetCell)
        % expand the time vectors of each pset to the full common
            tMinByAlign = min(cell2mat(cellfun(@(p) p.tMinForDataMean, makerow(psetCell), 'UniformOutput', false)), [], 2);
            tMaxByAlign = max(cell2mat(cellfun(@(p) p.tMaxForDataMean, makerow(psetCell), 'UniformOutput', false)), [], 2);
            
            prog = ProgressBar(numel(psetCell), 'Equalizing time vectors across psets');
            for iP = 1:numel(psetCell)
                prog.update(iP);
                psetCell{iP} = psetCell{iP}.manualSliceOrExpandTimeWindow(tMinByAlign, tMaxByAlign);
            end
            prog.finish();
        end
    end
end