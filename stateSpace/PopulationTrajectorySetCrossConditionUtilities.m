classdef PopulationTrajectorySetCrossConditionUtilities
   
    methods(Static)
        
        function psetDiff = computeDifferenceAlongAxis(pset, axisName, varargin)
            p = inputParser();
            p.addParameter('reverse', false, @islogical);
            p.parse(varargin{:});
            reverse = p.Results.reverse;
            
            cd = pset.conditionDescriptor;
            aIdx = cd.axisLookupByAttributes(axisName);
            
            pset.warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns()
            
            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            
            % adjust cd to reflect differences
            valueLists = pset.conditionDescriptor.generateAxisValueListsAsStrings(' ', true);
            valueList = valueLists{aIdx};
            if reverse
                diffValueList = cellfun(@(v1, v2) [v2 ' - ' v1], valueList(1:end-1), ...
                    valueList(2:end), 'UniformOutput', false);
            else
                diffValueList = cellfun(@(v1, v2) [v1 ' - ' v2], valueList(1:end-1), ...
                    valueList(2:end), 'UniformOutput', false);
            end
            b.conditionDescriptor = pset.conditionDescriptor.setAxisValueList(axisName, diffValueList);
            
            nConditionsNew = b.conditionDescriptor.nConditions;
            
            % adjust mean and sem to reflect difference between conditions
            [b.dataMean, b.dataSem] = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                % N x TA x C1 x C2 x ...
                dataMean = pset.buildNbyTAbyConditionAttributes('alignIdx', iAlign);
                
                diffTensor = diff(dataMean, 1, aIdx+2);
                if reverse
                    diffTensor = -diffTensor;
                end
                
                % N x TA x C
                diffReshape = reshape(diffTensor, [pset.nBases, pset.nTimeDataMean(iAlign), nConditionsNew]);
                
                % N x C x TA
                b.dataMean{iAlign} = permute(diffReshape, [1 3 2]);
                
                dataSem = pset.buildNbyTAbyConditionAttributes('type', 'sem', 'alignIdx', iAlign);
                
                % use sd1+2 = sqrt(sd1^2 / n1 + sd2^2 / n2) formula
                % which equals sem1+2 = sqrt(sem1^2 + sem2^2)
     
                % add squares of sem down the line and take sqrt
                mask1 = TensorUtils.maskByDimCellSelectAlongDimension(size(dataSem), aIdx+2, 1:(size(dataSem, aIdx+2)-1));
                mask2 = TensorUtils.maskByDimCellSelectAlongDimension(size(dataSem), aIdx+2, 2:size(dataSem, aIdx+2));
                
                semTensor = sqrt(dataSem(mask1{:}).^2 + dataSem(mask2{:}).^2);
                semReshape = reshape(semTensor, [pset.nBases, pset.nTimeDataMean(iAlign), nConditionsNew]);
                
                b.dataSem{iAlign} = permute(semReshape, [1 3 2]);
            end
            
            % update difference of trials scaled noise estimates so that we
            % can compute noise variance floors when projecting. since the
            % noise estimates are already scaled by 1/sqrt(2*nTrials), we
            % simply add them together to get the new scaled estimate
            scaledNoiseEstimate_NbyTAbyC = pset.dataDifferenceOfTrialsScaledNoiseEstimate;
            scaledNoiseEstimate_NbyTAbyAttr = reshape(scaledNoiseEstimate_NbyTAbyC, [pset.nBases, sum(pset.nTimeDataMean), makerow(pset.conditionsSize)]);
            
            % subtract adjacent noiseEstimates down the line (adding would
            % also be okay
            mask1 = TensorUtils.maskByDimCellSelectAlongDimension(size(scaledNoiseEstimate_NbyTAbyAttr), aIdx+2, 1:(size(scaledNoiseEstimate_NbyTAbyAttr, aIdx+2)-1));
            mask2 = TensorUtils.maskByDimCellSelectAlongDimension(size(scaledNoiseEstimate_NbyTAbyAttr), aIdx+2, 2:size(scaledNoiseEstimate_NbyTAbyAttr, aIdx+2));
           
            newScaledNoiseEstimate_NbyTAbyAttr = scaledNoiseEstimate_NbyTAbyAttr(mask1{:}) - scaledNoiseEstimate_NbyTAbyAttr(mask2{:});
            b.dataDifferenceOfTrialsScaledNoiseEstimate = reshape(newScaledNoiseEstimate_NbyTAbyAttr, ...
                [pset.nBases, sum(pset.nTimeDataMean), nConditionsNew]);
            
            % diff randomized data if present, recompute intervals
            if ~isempty(pset.dataMeanRandomized)
                [b.dataMeanRandomized, b.dataIntervalLow, b.dataIntervalHigh] = deal(cell(pset.nAlign, 1));
                
                % originally N x TA x C x R (where R is number of random samples)
                % randTensor is N x TA x C1 x C2 x ... x R
                randTensor = reshape(pset.dataMeanRandomized{iAlign}, ...
                    [pset.nBases, pset.nTimeDataMean(iAlign), makerow(pset.conditionsSize), pset.nRandomSamples]);

                diffTensor = diff(randTensor, 1, aIdx+2);
                if reverse
                    diffTensor = -diffTensor;
                end
                
                % back to N x TA x C x R
                b.dataMeanRandomized{iAlign} = reshape(diffTensor, [pset.nBases, pset.nTimeDataMean(iAlign), nConditionsNew, pset.nRandomSamples]);
                
                % recompute quantiles
                quantiles = quantile(b.dataMeanRandomized{iAlign}, ...
                    [pset.dataIntervalQuantileLow, pset.dataIntervalQuantileHigh], 4);
                b.dataIntervalQuantileLow{iAlign} = quantiles(:, :, :, 1);
                b.dataIntervalQuantileHigh{iAlign} = quantiles(:, :, :, 2);
            end
            
            psetDiff = b.buildManualWithTrialAveragedData();
        end
    end
    
end