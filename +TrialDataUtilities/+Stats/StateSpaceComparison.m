classdef StateSpaceComparison

    methods(Static)
        function [valueTensor, tvec, conditionDescriptorSansAxis] = dprime(pset, axisSpec, varargin)
            % compute (mu_each - mu_control) / sigma_control along each
            % axis for a specified idxControl (default 1)
            % valueTensor will match size of dataMean
            p = inputParser();
            p.addParameter('idxControl', 1, @islogical); % which position to use as the "control group"
            p.parse(varargin{:});

            axisIdx = pset.conditionDescriptor.axisLookupByAttributes(axisSpec);
            condSize = pset.conditionsSize;
            condSizeOther = TensorUtils.sizeOtherDims(pset.conditions, axisIdx);
            
            dataMean = pset.dataMean;
            dataStd = pset.computeDataStdRandomized
            pTensorByAlign = cell(pset.nAlign, 1);
            for iA = 1:pset.nAlign
                pTensorByAlign{iA} = nan(pset.nBases, pset.nTimeDataMean(iA), condSizeOther);

                if ~pset.alignValid(iA), continue; end
                prog = ProgressBar(pset.nBases, 'Computing %s for each basis for align %d', p.Results.functionName, iA);
                for iB = 1:pset.nBases
                    prog.update(iB);
                    if ~pset.basisValid(iB), continue; end

                    dataThis = reshape(dataByTrialGrouped(iB, iA, :), condSize);
                    [pTensor, ~, conditionDescriptorSansAxis] = fn(pset.getDataSourceForBasis(iB), ...
                        'axis', axisIdx, 'data', dataThis);
                    pTensorByAlign{iA}(iB, :) = pTensor(:);
                end
                prog.finish();
            end

            tvecByAlign = pset.tvecDataMean;
        end

    end
    
    methods(Static)
        
        %%%%%%%%%
        % kruskal wallis
        %%%%%%%%%

        function [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = kruskalWallisAlongAxisVsTimeEachAlign(pset, varargin)
            % [pTensorByAlign] = kruskalWallisAlongAxisVsTime(pset, varargin)
            % pTensorByAlign is nAlign x 1 cell, each contains N x T x
            % conditionDescriptorSansAxis.conditionsSize
            p = inputParser();
            p.addParamValue('axis', 1, @(x) true);
            p.addParamValue('removeZeroSpikeTrials', false, @islogical);
            p.parse(varargin{:});

            fn = @TrialDataUtilities.Stats.TimeseriesComparisonStatistics.kruskalWallisAlongAxisVsTime;
            [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.evaluateTimeseriesComparisonFnAlongAxisVsTimeEachAlign(pset, fn, ...
                'functionName', 'Kruskal Wallis', p.Results);
        end

        function [fractionSignificantByAlign, firstSignificantTimeByAlign, pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                fractionKruskalWallisSignificantAlongAxisVsTimeEachAlign(pset, varargin)
            % [fractionSignificantByAlign, firstSignificantTime, pTensorByAlign, conditionDescriptorSansAxis] = ...
            %    fractionKruskalWallisSignificantAlongAxisVsTime(tdca, varargin)
            %
            % let Cnew = conditionDescriptorSansAxis.conditionsSize;
            %
            % fractionSignificantByAlign: nAlign x 1 containing T x Cnew
            %
            % firstSignificantTimeAlign: nBases x nAlign. either zero or
            %   one value per row will be non-NaN, so that the time of first
            %   significance across all aligns will be present in the
            %   iAlign'th column where that basis first reaches significant
            p = inputParser();
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar);
            p.addParameter('nConsecutive', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.nConsecutiveDefault, @isscalar);
            p.addParameter('axis', 1, @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.kruskalWallisAlongAxisVsTimeEachAlign(pset, ...
                'axis', p.Results.axis);
            [fractionSignificantByAlign, firstSignificantTimeByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.findSignificantRunsInPValueTensor(pset, pTensorByAlign, ...
                'nConsecutive', p.Results.nConsecutive, 'alpha', p.Results.alpha);
        end

        %%%%%%%%%
        % anova
        %%%%%%%%%

        function [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = anovaAlongAxisVsTimeEachAlign(pset, varargin)
            % [pTensorByAlign] = kruskalWallisAlongAxisVsTime(pset, varargin)
            % pTensorByAlign is nAlign x 1 cell, each contains N x T x
            % conditionDescriptorSansAxis.conditionsSize
            p = inputParser();
            p.addParamValue('axis', 1, @(x) true);
            p.parse(varargin{:});

            fn = @TrialDataUtilities.Stats.TimeseriesComparisonStatistics.anovaAlongAxisVsTime;
            [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.evaluateTimeseriesComparisonFnAlongAxisVsTimeEachAlign(pset, fn, ...
                'functionName', 'anova', p.Results);
        end

        function [fractionSignificantByAlign, firstSignificantTimeByAlign, pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                fractionAnovaSignificantAlongAxisVsTimeEachAlign(pset, varargin)
            % [fractionSignificantByAlign, firstSignificantTime, pTensorByAlign, conditionDescriptorSansAxis] = ...
            %    fractionKruskalWallisSignificantAlongAxisVsTime(tdca, varargin)
            %
            % let Cnew = conditionDescriptorSansAxis.conditionsSize;
            %
            % fractionSignificantByAlign: nAlign x 1 containing T x Cnew
            %
            % firstSignificantTimeAlign: nBases x nAlign. either zero or
            %   one value per row will be non-NaN, so that the time of first
            %   significance across all aligns will be present in the
            %   iAlign'th column where that basis first reaches significant
            p = inputParser();
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar);
            p.addParameter('nConsecutive', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.nConsecutiveDefault, @isscalar);
            p.addParameter('axis', 1, @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.anovaAlongAxisVsTimeEachAlign(pset, ...
                'axis', p.Results.axis);
            [fractionSignificantByAlign, firstSignificantTimeByAlign] = ...
                TrialDataUtilities.Stats.StateSpaceComparison.findSignificantRunsInPValueTensor(pset, pTensorByAlign, ...
                'nConsecutive', p.Results.nConsecutive, 'alpha', p.Results.alpha);
        end
    end

    methods(Static) % Generic utilities for evaluating functions for all bases / aligns over time
        function [pTensorByAlign, conditionDescriptorSansAxis, tvecByAlign] = ...
                evaluateTimeseriesComparisonFnAlongAxisVsTimeEachAlign(pset, fn, varargin)
            % [pTensorByAlign] = kruskalWallisAlongAxisVsTime(pset, varargin)
            % pTensorByAlign is nAlign x 1 cell, each contains N x T x
            % conditionDescriptorSansAxis.conditionsSize
            p = inputParser();
            p.addParamValue('axis', 1, @(x) true);
            p.addParamValue('functionName', 'comparison function', @ischar);
            p.parse(varargin{:});

            axisIdx = pset.conditionDescriptor.axisLookupByAttributes(p.Results.axis);
            condSize = pset.conditionsSize;
            condSizeOther = TensorUtils.sizeOtherDims(pset.conditions, axisIdx);

            % nBases x nAlign
            dataByTrialGrouped = pset.buildDataByTrialGroupedWithDataMeanTimeVector();
            pTensorByAlign = cell(pset.nAlign, 1);
            for iA = 1:pset.nAlign
                pTensorByAlign{iA} = nan(pset.nBases, pset.nTimeDataMean(iA), condSizeOther);

                if ~pset.alignValid(iA), continue; end
                prog = ProgressBar(pset.nBases, 'Computing %s for each basis for align %d', p.Results.functionName, iA);
                for iB = 1:pset.nBases
                    prog.update(iB);
                    if ~pset.basisValid(iB), continue; end

                    dataThis = reshape(dataByTrialGrouped(iB, iA, :), condSize);
                    [pTensor, ~, conditionDescriptorSansAxis] = fn(pset.getDataSourceForBasis(iB), ...
                        'axis', axisIdx, 'data', dataThis);
                    pTensorByAlign{iA}(iB, :) = pTensor(:);
                end
                prog.finish();
            end

            tvecByAlign = pset.tvecDataMean;
        end

        function [fractionSignificantByAlign, firstSignificantTimeByAlign] = findSignificantRunsInPValueTensor(pset, pTensorByAlign, varargin)
            p = inputParser();
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar);
            p.addParameter('nConsecutive', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.nConsecutiveDefault, @isscalar);
            p.addParameter('axis', 1, @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % concatenate over alignments
            % N x TA x conditionsSizeRemaining
            [pTensorCat, whichAlign] = TensorUtils.catWhich(2, pTensorByAlign{:});
            tvecCat = cat(1, pset.tvecDataMean{:});
            sigTensor = pTensorCat <= p.Results.alpha;
            [changeIdx, inRunMask] = TrialDataUtilities.Data.findFirstConsecutiveRun(sigTensor, p.Results.nConsecutive, 2, 1);

            % find crossing times
            firstSignificantTimeByAlign = nan(pset.nBases, pset.nAlign);
            crossTimes = TensorUtils.selectAlongDimensionWithNaNs(tvecCat, 1, changeIdx);
            crossWhichAlign = TensorUtils.selectAlongDimensionWithNaNs(whichAlign, 1, changeIdx);
            mask = ~isnan(crossTimes);
            inds = sub2ind(size(firstSignificantTimeByAlign), find(mask), crossWhichAlign(mask));
            firstSignificantTimeByAlign(inds) = crossTimes(mask);

            % take mean over this over neurons and squeeze out first dim (size N)
            % this will be TA x Cnew
            fractionSignificantCat = TensorUtils.squeezeDims(mean(inRunMask, 1), 1);

            fractionSignificantByAlign = TensorUtils.splitAlongDimension(fractionSignificantCat, 1, pset.nTimeDataMean);
        end
    end
end
