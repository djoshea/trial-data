classdef StateSpaceComparisons 

    methods(Static)
        function [pTensor] = kruskalWallisAlongAxisVsTime(pset, varargin) 
            % [pTensor] = kruskalWallisAlongAxisVsTime(pset, varargin) 
            % pTensor is N x C x TA 
            
            p = inputParser();
            p.addParamValue('axis', 1, @(x) true);
            p.parse(varargin{:});
           
            pTensor = nan(pset.nBases, pset.nConditions, sum(pset.nTimeDataMean));
            axisIdx = pset.conditionDescriptor.axisLookupByAttributes(p.Results.axis);

            prog = ProgressBar(pset.nBases, 'Evaluating for each basis');
            for iB = 1:pset.nBases
                prog.update(iB);
                if ~pset.basisValid(iB)
                    continue;
                end

                [taByCByTrials] = TensorUtils.squeezeDims(pset.buildDataNbyTAbyCbyTrials('basisIdx', iB), 1);
                pThis = TrialData.TimeseriesComparisonStatistics.kruskalWallisAlongAxisVsTime([], 'axis', axisIdx, 'data', jk:wq

                                
            end
            prog.finish();
       end
    end

end
