classdef TimeseriesStatistics
    
    methods(Static) % Effect size methods
        function [dprimeTensor, dprimeCI, tvec] = dPrimeAlongAxisVsTime(tdca, varargin)
            p = inputParser;
            p.addParameter('alpha', 0.05, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;
            
            import(getPackageImportString);
            [misc, dprimeTensor, dprimeHigh, dprimeLow] = TimeseriesStatistics.evaluateComparisonAlongAxisVsTime(tdca, @dprimeFn, p.Unmatched);
            tvec = misc.tvec;
            
            dprimeCI = cat(1, shiftdim(dprimeLow, -1), shiftdim(dprimeHigh, -1));
            
            function [dprime, dprimeHigh, dprimeLow] = dprimeFn(inCell)
                assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
                delta = nanmean(inCell{1}, 1) - nanmean(inCell{2}, 1);
                v1 = nanvar(inCell{1}, 0, 1);
                v2 = nanvar(inCell{2}, 0, 1);
                sd = sqrt(0.5* (v1 + v2));
                dprime = delta ./ sd;
                n1 = sum(~isnan(inCell{1}), 1);
                n2 = sum(~isnan(inCell{2}), 1);
                df = n1+n2-2;
                se = sqrt(v1/n1 + v2/n2);
                ciFac = -tinv(alpha/2, df);
                dprimeHigh = dprime + se.*ciFac;
                dprimeLow = dprime - se.*ciFac;
            end        
        end
        
        function [gTensor, gCI, tvec] = hedgesGAlongAxisVsTime(tdca, varargin)
            p = inputParser;
            p.addParameter('alpha', 0.05, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;
            
            import(getPackageImportString);
            [misc, gTensor, gHigh, gLow] = TimeseriesStatistics.evaluateComparisonAlongAxisVsTime(tdca, @hedgesGFn, p.Unmatched);
            tvec = misc.tvec;
            
            gCI = cat(1, shiftdim(gLow, -1), shiftdim(gHigh, -1));
            
            function [g, gHigh, gLow] = hedgesGFn(inCell)
                assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
                stats = mes(inCell{1}, inCell{2}, 'hedgesg', 'confLevel', 1-alpha);
                g = stats.hedgesg;
                gHigh = stats.hedgesgCi(1, :);
                gLow = stats.hedgesgCi(2, :);
            end        
        end
        
        function [gTensor, gCI, tvec] = meanDifferenceAlongAxisVsTime(tdca, varargin)
            p = inputParser;
            p.addParameter('alpha', 0.05, @isscalar); % 1-a defines confidence intervals
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            alpha = p.Results.alpha;
            
            import(getPackageImportString);
            [misc, gTensor, gHigh, gLow] = TimeseriesStatistics.evaluateComparisonAlongAxisVsTime(tdca, @hedgesGFn, p.Unmatched);
            tvec = misc.tvec;
            
            gCI = cat(1, shiftdim(gLow, -1), shiftdim(gHigh, -1));
            
            function [g, gHigh, gLow] = hedgesGFn(inCell)
                assert(numel(inCell) == 2, 'd'' only supported for axes with 2 conditions');
                stats = mes(inCell{1}, inCell{2}, 'md', 'confLevel', 1-alpha);
                g = stats.md;
                gHigh = stats.mdCi(1, :);
                gLow = stats.mdCi(2, :);
            end        
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
            
            axisIdx = tdca.conditionInfo.axisLookupByAttributes(p.Results.axis);
            
            % dataAxisFirst is permuted from conditionsSize to put selected
            % axis first
            if isempty(p.Results.data)
                [data, misc.tvec] = tdca.getAnalogAsMatrixGrouped(p.Results.name, p.Unmatched);
            else
                data = p.Results.data;
                misc.tvec = p.Results.tvec;
            end
            
            otherDims = TensorUtils.otherDims(size(data), axisIdx);
            dimPerm = [axisIdx, otherDims];
            dataAxisFirst = permute(data, dimPerm);
%             nAlongAxis = size(dataAxisFirst, 1);
            sizeOtherAxes = TensorUtils.sizeOtherDims(dataAxisFirst, 1);
            nOtherAxes = prod(sizeOtherAxes);
            
            % in case tvec isn't specified, grab the first non-empty data
            % and determine the number of timepoints
            Tmat = cellfun(@(x) size(x, 2), dataAxisFirst);
            T = max(Tmat(:));
            
            nArg = nargout(fn);
            [varargout{1:nArg}] = deal(nan([T, sizeOtherAxes]));
            
            for iOtherAxes = 1:nOtherAxes
                % nAlongAxis x 1
                dataThis = dataAxisFirst(:, iOtherAxes);
                [argThis{1:nArg}] = fn(dataThis);
                
                for iArg = 1:nArg
                    varargout{iArg}(:, iOtherAxes) = argThis{iArg};
                end
            end        
        end
        
        
    end
    
end

% for iOtherAxes = 1:nOtherAxes
%                 % nAlongAxis x 1
%                 dataThis = dataAxisFirst(:, iOtherAxes);
%                 
%                 % nTrialsTotal x T, nTrialsTotal x 1
%                 [dataCat, whichCondition] = TensorUtils.catWhich(1, dataThis{:});
%                
%                 for iT = 1:T
%                     keep = ~isnan(dataCat(:, iT));
%                     pValTensor(iT, iOtherAxes) = kruskalwallis(dataCat(keep, iT), whichCondition(keep), 'off');
%                 end
%             end