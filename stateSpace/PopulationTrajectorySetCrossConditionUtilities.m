classdef PopulationTrajectorySetCrossConditionUtilities
   
    % combining / weighting along condition axes
    methods(Static)
        function psetDiff = computeDifferenceAlongAxis(pset, axisName, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.addParameter('reverse', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            reverse = p.Results.reverse;
            nAlongAxis = pset.conditionsSize(aIdx);
            
            % generate new names from differences
            if isempty(p.Results.newNamesAlongAxis)
                valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings(' ', true);
                valueList = valueLists{aIdx};
                if reverse
                    newNamesAlongAxis = cellfun(@(v1, v2) [v2 ' - ' v1], valueList(1:end-1), ...
                        valueList(2:end), 'UniformOutput', false);
                else
                    newNamesAlongAxis = cellfun(@(v1, v2) [v1 ' - ' v2], valueList(1:end-1), ...
                        valueList(2:end), 'UniformOutput', false);
                end
            else
                newNamesAlongAxis = p.Results.newNamesAlongAxis;
            end
            
            % generate differencing matrix
            if reverse
                % compute conditions 2-1, 3-2, 4-3, etc.
                mat = diag(-onesvec(nAlongAxis)) + diag(onesvec(nAlongAxis-1), 1);
                wNbyO = mat(1:nAlongAxis-1, :);
            else
                % compute conditions 1-2, 2-3, 3-4, etc.
                mat = diag(onesvec(nAlongAxis)) + diag(-onesvec(nAlongAxis-1), 1);
                wNbyO = mat(1:nAlongAxis-1, :);
            end
            
            psetDiff = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newNamesAlongAxis', newNamesAlongAxis, p.Unmatched);
        end
        
        function psetDiff = subtractOneConditionFromOthersAlongAxis(pset, axisName, conditionToSubtract, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            nAlongAxis = pset.conditionsSize(aIdx);
            
            idxKeep = setdiff(1:nAlongAxis, conditionToSubtract);
            
            % generate new names from differences
            if isempty(p.Results.newNamesAlongAxis)
                valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings(' ', true);
                valueList = valueLists{aIdx}(idxKeep);
                valueSubtract = valueLists{aIdx}{conditionToSubtract};
                
                newNamesAlongAxis = cellfun(@(v) [v ' - ' valueSubtract], valueList, ...
                    'UniformOutput', false);
            else
                newNamesAlongAxis = p.Results.newNamesAlongAxis;
            end
            
            % generate differencing matrix
            mat = eye(nAlongAxis);
            mat(:, conditionToSubtract) = -1;
            wNbyO = mat(idxKeep, :);
            
            psetDiff = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newNamesAlongAxis', newNamesAlongAxis, p.Unmatched);
        end
        
        function psetMean = computeMeanAlongAxis(pset, axisName, varargin)
            % more parameters available in applyLinearCombinationAlongConditionAxis
            p = inputParser();
            p.addParameter('newNameAlongAxis', '', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            nAlongAxis = pset.conditionsSize(aIdx);
            
            % generate new names from differences
            if isempty(p.Results.newNameAlongAxis)
                newNamesAlongAxis = {sprintf('Mean Over %s', pset.conditionDescriptor.axisNames{aIdx}) }; 
            else
                newNamesAlongAxis = {p.Results.newNameAlongAxis};
            end

            % normalization is done by normalizeCoeffientsByNumConditions
            % below
            wNbyO = ones(1, nAlongAxis);
            
            psetMean = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, ...
                axisName, wNbyO, 'newNamesAlongAxis', newNamesAlongAxis, ...
                'replaceNaNWithZero', true, 'normalizeCoeffientsByNumConditions', true, p.Unmatched);
        end
        
        function psetReweighted = applyLinearCombinationAlongConditionAxis(pset, axisName, weightsNewCByOldC, varargin)
            p = inputParser();
            p.addParameter('newNamesAlongAxis', {}, @iscellstr);
            p.addParameter('removeAxis', false, @islogical);
            p.addParameter('conditionAppearanceFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            
            % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis
            % this enables nanmean like computations
            p.addParameter('normalizeCoeffientsByNumConditions', false, @islogical); 
            p.parse(varargin{:});
            
            aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisName);
            
            wNbyO = weightsNewCByOldC;
            cOld = pset.conditionsSize(aIdx);
            assert(size(wNbyO, 2) == cOld, ...
                'Weighting matrix must have same column count as number of existing conditions along axis (%d)', cOld);
            
            cNewAlongAxis = size(wNbyO, 1);
            
            newNamesAlongAxis = p.Results.newNamesAlongAxis;
            if isempty(newNamesAlongAxis)
                newNamesAlongAxis = arrayfun(@(i) sprintf('Condition Combination %d', i), 1:cNewAlongAxis, 'UniformOutput', false);
            end
            assert(numel(newNamesAlongAxis) == cNewAlongAxis, 'newNamesAlongAxis must have numel == number of new conditions (%d)', cNewAlongAxis);
            
            pset.warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns();
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % setup new condition descriptor, optionally drop the axis
            % we're combining along, if the new size is 1
            newCD = pset.conditionDescriptor.setAxisValueList(axisName, newNamesAlongAxis);
            if p.Results.removeAxis
                assert(cNewAlongAxis == 1, 'New condition count along axis must be 1 in order to removeAxis');
                newCD = newCD.removeAxis(aIdx);
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
                tensorMean = pset.buildNbyTAbyConditionAttributes('alignIdx', iAlign);
                tensorMeanReweighted = TensorUtils.linearCombinationAlongDimension(tensorMean, aIdx+2, wNbyO, ...
                    'replaceNaNWithZero', p.Results.replaceNaNWithZero, ...
                    'normalizeCoeffientsByNumConditions', p.Results.normalizeCoeffientsByNumConditions);
                % back to N x C x TA
                b.dataMean{iAlign} = permute(tensorMeanReweighted(:, :, :), [1 3 2]);
                
                % build N x TA x C1 x C2 x ...
                % use sd1+2 = sqrt(sd1^2 / n1 + sd2^2 / n2) formula
                % which here means semNew = sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                tensorSem = pset.buildNbyTAbyConditionAttributes('type', 'sem', 'alignIdx', iAlign);
                tensorSemReweighted = sqrt( TensorUtils.linearCombinationAlongDimension(tensorSem.^2, aIdx+2, abs(wNbyO)) );
                % back to N x C x TA
                b.dataSem{iAlign} = permute(tensorSemReweighted(:, :, :), [1 3 2]);
            end
            
            % update difference of trials scaled noise estimates so that we
            % can compute noise variance floors when projecting. since the
            % noise estimates are already scaled by 1/sqrt(2*nTrials), we
            % simply add them together to get the new scaled estimate
            scaledNoiseEstimate_NbyTAbyC = pset.dataDifferenceOfTrialsScaledNoiseEstimate;
            scaledNoiseEstimate_NbyTAbyAttr = reshape(scaledNoiseEstimate_NbyTAbyC, [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize)]);
            newScaledNoiseEstimate_NbyTAbyAttr = TensorUtils.linearCombinationAlongDimension(scaledNoiseEstimate_NbyTAbyAttr, aIdx+2, abs(wNbyO));
            b.dataDifferenceOfTrialsScaledNoiseEstimate = reshape(newScaledNoiseEstimate_NbyTAbyAttr, ...
                [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew]);
            
            % diff randomized data if present, recompute intervals
            if ~isempty(pset.dataMeanRandomized)
                [b.dataMeanRandomized, b.dataSemRandomized] = deal(cell(pset.nAlign, 1));
                for iAlign = 1:pset.nAlign
                    % dataMeanRandomized is N x C x TA x R (where R is number of random samples)
                    % randTensor is N x C1 x C2 x ... x TA x R
                    meanTensor = reshape(pset.dataMeanRandomized{iAlign}, ...
                        [pset.nBases, makerow(pset.conditionsSize), pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    
                    meanTensorReweighted = TensorUtils.linearCombinationAlongDimension(meanTensor, aIdx+1, wNbyO); 
                    
                    % back to N x C x TA x R
                    b.dataMeanRandomized{iAlign} = reshape(meanTensorReweighted, [pset.nBases, nConditionsNew, pset.nTimeDataMean(iAlign), pset.nRandomSamples]);

                    % ensure invalid bases remain invalid
                    b.dataMeanRandomized{iAlign}(~pset.basisValid, :, :, :) = NaN;
                    
                    % dataSemRandomized is N x C x TA x R (where R is number of random samples)
                    % randTensor is N x C1 x C2 x ... x TA x R
                    semTensor = reshape(pset.dataSemRandomized{iAlign}, ...
                        [pset.nBases, makerow(pset.conditionsSize), pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    
                    semTensorReweighted = TensorUtils.linearCombinationAlongDimension(semTensor, aIdx+1, wNbyO); 

                    % back to N x C x TA x R
                    b.dataSemRandomized{iAlign} = reshape(semTensorReweighted, [pset.nBases, nConditionsNew, pset.nTimeDataMean(iAlign), pset.nRandomSamples]);
                    b.dataSemRandomized{iAlign}(~pset.basisValid, :, :, :) = NaN;
                end
            end
            
            b.trialLists = {}; % no longer relevant
            
            % A x N x C
            % for data valid, we need all input conditions to be valid for
            % output conditions to be valid, so we change wNbyO such that
            % the linear combination will be 1 iff all bases that
            % contribute to that output are valid.
            wNbyO_forValid = bsxfun(@rdivide, wNbyO ~= 0, sum(wNbyO ~= 0, 2));
            [dataValidTensor, cdims] = TensorUtils.reshapeDimsInPlace(pset.dataValid, 3, pset.conditionsSize);
            b.dataValid = TensorUtils.flattenDimsInPlace(TensorUtils.linearCombinationAlongDimension(...
                dataValidTensor, aIdx+2, wNbyO_forValid) == 1, cdims);
            
            % sum trials from all included ocnditions
            [dataNTrialsTensor, cdims] = TensorUtils.reshapeDimsInPlace(pset.dataNTrials, 3, pset.conditionsSize);
            b.dataNTrials = TensorUtils.flattenDimsInPlace(TensorUtils.linearCombinationAlongDimension(...
                dataNTrialsTensor, aIdx+2, wNbyO ~= 0), cdims);

            % shrink the time windows over all considered conditions
            [tMinValidOld, cdims] = TensorUtils.reshapeDimsInPlace(pset.tMinValidByAlignBasisCondition, 3, pset.conditionsSize);
            tMinValidCellByNew = arrayfun(@(iNew) max(TensorUtils.selectAlongDimension(tMinValidOld, aIdx+2, wNbyO(iNew, :) ~= 0), [], aIdx+2), ...
                1:cNewAlongAxis, 'UniformOutput', false);
            tMinValidNew = cat(aIdx+2, tMinValidCellByNew{:});
            b.tMinValidByAlignBasisCondition = TensorUtils.flattenDimsInPlace(tMinValidNew, cdims);
            
            [tMaxValidOld, cdims] = TensorUtils.reshapeDimsInPlace(pset.tMaxValidByAlignBasisCondition, 3, pset.conditionsSize);
            tMaxValidCellByNew = arrayfun(@(iNew) min(TensorUtils.selectAlongDimension(tMaxValidOld, aIdx+2, wNbyO(iNew, :) ~= 0), [], aIdx+2), ...
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
                    b.alignSummaryData{iSource, iAlign} = b.alignSummaryData{iSource, iAlign}.combineSetsOfConditions(...
                        b.conditionDescriptor, conditionIdxSetsTensor(:));
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
            
            p.parse(varargin{:});
            
            pset = psetCell{1};
            
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            cd = pset.conditionDescriptor;
            cd = cd.addAttribute(axisName, 'valueList', axisValueList);
            cd = cd.addAxis(axisName, 'valueList', axisValueList);
            
            % update the condition appearance fn if specified
            if ~ismember('conditionAppearanceFn', p.UsingDefaults) % we dont just check isempty b/c the user may or may not want to set it to empty
                cd.appearanceFn = p.Results.conditionAppearanceFn;
            end
            
            b.conditionDescriptor = cd;
            cAxis = cd.nAxes;
            
            cSize = pset.conditionsSize;
            cSize(end+1) = 1;
            
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

            debug('Concatenating trial-averaged data\n');
            for iAlign = 1:pset.nAlign
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
            
            % N x T x C
            b.dataDifferenceOfTrialsScaledNoiseEstimate = catConditionsFlat(psetCell, @(p) p.dataDifferenceOfTrialsScaledNoiseEstimate, 3, cSize, cAxis);
        
            % adjust alignSummaryData
            % N x A
            temp = pset.alignSummaryData; %#ok<NASGU> % request up front to trigger computation before progress bar
            prog = ProgressBar(pset.nBases, 'Aggregating AlignSummary by basis');
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                for iAlign = 1:pset.nAlign
                    clear alignSummarySet;
                    for iP = numel(psetCell):-1:1
                        alignSummarySet(iP) = psetCell{iP}.lookupAlignSummaryDataForBasisAlign(iBasis, iAlign);
                    end
                    b.alignSummaryData{iBasis, iAlign} = ...
                        AlignSummary.aggregateByConcatenatingConditionsAlongNewAxis(alignSummarySet, cd, ...
                        'aggregateMarks', p.Results.aggregateMarks, 'aggregateIntervals', p.Results.aggregateIntervals);
                end
            end
            prog.finish();

            hasDataRandomized = cellfun(@(p) p.hasDataRandomized, psetCell);
            if all(hasDataRandomized)
                debug('Concatenating data randomized\n');
                for iAlign = 1:pset.nAlign
                    b.dataMeanRandomized{iAlign} = catConditionsFlat(psetCell, @(p) p.dataMeanRandomized{iAlign}, 2, cSize, cAxis);
                    b.dataSemRandomized{iAlign} = catConditionsFlat(psetCell, @(p) p.dataSemRandomized{iAlign}, 2, cSize, cAxis);
                end
            end
            
            psetCat = b.buildManualWithTrialAveragedData();

            function res = catConditionsFlat(objCell, accessFn, conditionDim, conditionsSize, conditionAxis)
                % grab from each obj in cell
                tensorCell = cellfun(accessFn, objCell, 'UniformOutput', false);
                
                % reshape to make conditions a tensor rather than flat
                sz = TensorUtils.sizeNDims(tensorCell{1}, conditionDim);
                newSz = [sz(1:conditionDim-1), conditionsSize, sz(conditionDim+1:end)];
                tensorCell = cellfun(@(x) reshape(x, newSz), tensorCell, 'UniformOutput', false);
                
                % concatenate along the condition axis
                catAxis = conditionDim + conditionAxis-1;
                catTensor = cat(catAxis, tensorCell{:});
                
                % reshape the result to be flat again
                catSz = sz;
                catSz(conditionDim) = catSz(conditionDim) * numel(objCell);
                res = reshape(catTensor, catSz);
            end
            
        end
    end
    
end