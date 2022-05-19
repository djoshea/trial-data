classdef TimeseriesComparisonStatistics

    properties(Constant)
        nConsecutiveDefault = 5;
        alphaDefault = 0.05;
    end

    methods(Static) % Effect size methods
        function [dprimeTensor, dprimeCI, tvec, conditionDescriptorSansAxis] = dPrimeAlongAxisVsTime(tdca, varargin)
            import(getPackageImportString);
            p = inputParser;
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;

            [misc, dprimeTensor, dprimeHigh, dprimeLow] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.evaluateComparisonAlongAxisVsTime(...
                tdca, @(in) TrialDataUtilities.Stats.TimeseriesComparisonStatistics.dprimeFn(alpha, in), p.Unmatched);
            tvec = misc.tvec;

            dprimeCI = cat(1, shiftdim(dprimeLow, -1), shiftdim(dprimeHigh, -1));
            conditionDescriptorSansAxis = misc.conditionDescriptorSansAxis;
        end

        function [gTensor, gCI, tvec, conditionDescriptorSansAxis] = hedgesGAlongAxisVsTime(tdca, varargin)
            import(getPackageImportString);
            p = inputParser;
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;

            [misc, gTensor, gHigh, gLow] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.evaluateComparisonAlongAxisVsTime(...
                tdca, @(in) TrialDataUtilities.Stats.TimeseriesComparisonStatistics.hedgesGFn(alpha, in), p.Unmatched);
            tvec = misc.tvec;

            gCI = cat(1, shiftdim(gLow, -1), shiftdim(gHigh, -1));

            conditionDescriptorSansAxis = misc.conditionDescriptorSansAxis;
        end

        function [gTensor, gCI, tvec, conditionDescriptorSansAxis] = meanDifferenceAlongAxisVsTime(tdca, varargin)
            import(getPackageImportString);
            p = inputParser;
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;

            [misc, gTensor, gHigh, gLow] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.evaluateComparisonAlongAxisVsTime(...
                tdca, @(in) TrialDataUtilities.Stats.TimeseriesComparisonStatistics.meanDiffFn(alpha, in), p.Unmatched);
            tvec = misc.tvec;

            gCI = cat(1, shiftdim(gLow, -1), shiftdim(gHigh, -1));
            conditionDescriptorSansAxis = misc.conditionDescriptorSansAxis;
        end
    end

    methods(Static) % difference of means hypothesis testing
        function [pValTensor, tvec, conditionDescriptorSansAxis] = kruskalWallisAlongAxisVsTime(tdca, varargin)
            % pValTensor will be T x size(other condition axes) x neurons
            import(getPackageImportString);
            p = inputParser;
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [misc, pValTensor] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.evaluateComparisonAlongAxisVsTime(tdca, @TrialDataUtilities.Stats.TimeseriesComparisonStatistics.kwFun, p.Unmatched);
            tvec = misc.tvec;
            conditionDescriptorSansAxis = misc.conditionDescriptorSansAxis;
        end

        function [pValTensor, tvec, conditionDescriptorSansAxis] = anovaAlongAxisVsTime(tdca, varargin)
            % pValTensor will be T x size(other condition axes)
            import(getPackageImportString);
            p = inputParser;
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [misc, pValTensor] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.evaluateComparisonAlongAxisVsTime(tdca, @TrialDataUtilities.Stats.TimeseriesComparisonStatistics.anova1Fun, p.Unmatched);
            tvec = misc.tvec;
            conditionDescriptorSansAxis = misc.conditionDescriptorSansAxis;
        end
    end

    methods(Static) % find first time of effect size divergence above threshold
        function [crossTimes, gTensor, gCI, tvec, conditionDescriptorSansAxis] = findTimeHedgesGAboveThreshold(tdca, varargin)
            import(getPackageImportString);
            p = inputParser();
            p.addParameter('thresh', 1, @isscalar);
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar);
            p.addParameter('nConsecutive', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.nConsecutiveDefault, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % will be T x size(other condition axes)
            [gTensor, gCI, tvec, conditionDescriptorSansAxis] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.hedgesGAlongAxisVsTime(tdca, 'alpha', p.Results.alpha, p.Unmatched);
            if isempty(tvec)
                error('tvec parameter must be provided if data passed in as parameter');
            end

            isAboveThresh = TensorUtils.squeezeDims(min(abs(gCI), [], 1) >= p.Results.thresh, 1);
            changeIdx = TrialDataUtilities.Data.findFirstConsecutiveRun(isAboveThresh, p.Results.nConsecutive, 1);

            crossTimes = TensorUtils.selectAlongDimensionWithNaNs(makecol(tvec), 1, changeIdx);
            crossTimes = reshape(crossTimes, size(changeIdx));
        end

        function [crossTimes, pValTensor, tvec, conditionDescriptorSansAxis] = findTimeKruskalWallisSignificantAlongAxis(tdca, varargin)
            import(getPackageImportString);

            p = inputParser();
            p.addParameter('alpha', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.alphaDefault, @isscalar);
            p.addParameter('nConsecutive', TrialDataUtilities.Stats.TimeseriesComparisonStatistics.nConsecutiveDefault, @isscalar);
            p.addParameter('axis', 1, @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            % will be T x size(other condition axes)
            [pValTensor, tvec, conditionDescriptorSansAxis] = TrialDataUtilities.Stats.TimeseriesComparisonStatistics.kruskalWallisAlongAxisVsTime(...
                tdca, 'axis', p.Results.axis, p.Unmatched);
            if isempty(tvec)
                error('tvec parameter must be provided if data passed in as parameter');
            end

            isSignificant = TensorUtils.squeezeDims(pValTensor <= p.Results.alpha, 1);
            changeIdx = TrialDataUtilities.Data.findFirstConsecutiveRun(isSignificant, p.Results.nConsecutive, 1);

            crossTimes = TensorUtils.selectAlongDimensionWithNaNs(makecol(tvec), 1, changeIdx);
            crossTimes = reshape(crossTimes, size(changeIdx));
        end

    end

    methods(Static) % internal utility functions
        function [misc, varargout] = evaluateComparisonAlongAxisVsTime(tdca, fn, varargin)
            % at each time point, computes the pvalue of a ttest2 or anova1
            % for samples from analog or spike channel chName over the condition axis
            % axisSpec.
            p = inputParser();
            p.addParameter('axis', 1, @(x) true);
            p.addParameter('name', '', @ischar);
            p.addParameter('data', {}, @iscell);
            p.addParameter('tvec', {}, @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if ~isscalar(p.Results.axis) || ~isempty(tdca)
                axisIdx = tdca.conditionInfo.axisLookupByAttributes(p.Results.axis);
            else
                axisIdx = p.Results.axis;
            end

            % dataAxisFirst is permuted from conditionsSize to put selected
            % axis first
            if isempty(p.Results.data)
                cdesc = tdca.channelDescriptorsByName.(p.Results.name);
                if isa(cdesc, 'AnalogChannelDescriptor')
                    [data, misc.tvec] = tdca.getAnalogAsMatrixGrouped(p.Results.name, p.Unmatched);
                else
                    [data, misc.tvec] = tdca.getSpikeRateFilteredAsMatrixGrouped(p.Results.name, p.Unmatched);
                end
            else
                data = p.Results.data;
                misc.tvec = p.Results.tvec;
            end

            otherDims = TensorUtils.otherDims(size(data), axisIdx);
            dimPerm = [axisIdx, otherDims];
            dataAxisFirst = permute(data, dimPerm);
%             nAlongAxis = size(dataAxisFirst, 1);
            sizeOtherAxes = TensorUtils.sizeOtherDims(dataAxisFirst, 1);
            numelOtherAxes = prod(sizeOtherAxes);

            % in case tvec isn't specified, grab the first non-empty data
            % and determine the number of timepoints
            Tmat = cellfun(@(x) size(x, 2), dataAxisFirst);
            T = max(Tmat(:));
            C = size(dataAxisFirst{1}, 3);

            nArg = nargout - 1;
%             nArg = nargout(fn);
            [varargout{1:nArg}] = deal(nan([T, numelOtherAxes, C]));

            for iOtherAxes = 1:numelOtherAxes
                % nAlongAxis x 1 { tRials x Time x Neurons }
                dataThis = dataAxisFirst(:, iOtherAxes);
                [argThis{1:nArg}] = fn(dataThis); % --> Time x Neurons

                for iArg = 1:nArg
                    % this is wrong as the output is already T x C
%                     argPerm = ipermute(argThis{iArg}, dimPerm); % reorient the way it was
                    varargout{iArg}(:, iOtherAxes, :) = permute(argThis{iArg}, [1 3 2]);
                end
            end

            % re-expand other axes in place inside varargout
            for iArg = 1:nArg
                varargout{iArg} = reshape(varargout{iArg}, [T, sizeOtherAxes, C]);
            end

            if ~isempty(tdca)
                misc.conditionDescriptorSansAxis = tdca.conditionInfo.getConditionDescriptor.removeAxis(axisIdx);
            else
                misc.conditionDescriptorSansAxis = [];
            end
        end


    end

    methods(Static)
        function [dprime, dprimeHigh, dprimeLow] = dprimeFn(alpha, inCell)
            assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
            delta = nanmean(inCell{1}, 1) - nanmean(inCell{2}, 1);
            v1 = nanvar(inCell{1}, 0, 1);
            v2 = nanvar(inCell{2}, 0, 1);
            sd = sqrt(0.5 * (v1 + v2));
            dprime = delta ./ sd;
            n1 = sum(~isnan(inCell{1}), 1);
            n2 = sum(~isnan(inCell{2}), 1);
            df = n1+n2-2;
            se = sqrt(v1/n1 + v2/n2);
            ciFac = -tinv(alpha/2, df);
            dprimeHigh = dprime + se.*ciFac;
            dprimeLow = dprime - se.*ciFac;
        end

        function [g, gHigh, gLow] = hedgesGFn(alpha, inCell)
            assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
            stats = mes(inCell{1}, inCell{2}, 'hedgesg', 'confLevel', 1-alpha);
            g = stats.hedgesg;
            gHigh = stats.hedgesgCi(1, :);
            gLow = stats.hedgesgCi(2, :);
        end

        function [g, gHigh, gLow] = meanDiffFn(alpha, inCell)
            assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
            stats = mes(inCell{1}, inCell{2}, 'md', 'confLevel', 1-alpha);
            g = stats.md;
            gHigh = stats.mdCi(1, :);
            gLow = stats.mdCi(2, :);
        end

        function pValVsTime = kwFun(inCell)
            % supports channels on dim 3
            [dataCat, group] = TensorUtils.catWhich(1, inCell{:});
            T = size(dataCat, 2); % tRials x Time (x channels)
            C = size(dataCat, 3);
            pValVsTime = nan(T, C);
            for iT = 1:T
                for iC = 1:C
                    keep = ~isnan(dataCat(:, iT, iC));
                    pValVsTime(iT, iC) = kruskalwallis(dataCat(keep, iT, iC), group(keep), 'off');
                end
            end
        end

        function pValVsTime = anova1Fun(inCell)
            % supports channels on dim 3
            [dataCat, group] = TensorUtils.catWhich(1, inCell{:});
            T = size(dataCat, 2);
            C = size(dataCat, 3);
            pValVsTime = nan(T, C);
            for iT = 1:T
                for iC = 1:C
                    keep = ~isnan(dataCat(:, iT, iC));
                    pValVsTime(iT, iC) = anova1(dataCat(keep, iT, iC), group(keep), 'off');
                end
            end
        end
    end

end
