classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
        lambda = 1e-6;
        
        axisCombinations % see TrialDataUtilities.DPCA.dpca for dimListsToCombineList argument
        
        makeOrthogonal = false;
    end
    
    methods
        function proj = set.axisCombinations(proj, v)
            assert(iscell(v), 'Axis Combinations must be a cell of vectors');
            proj.axisCombinations = v;
        end
    end

    methods
        function proj = ProjDPCA_NonOrthogonal(varargin)
            proj = proj@StateSpaceProjection();
        end
        
        function setBuildFromConditionIdx(proj, varargin) %#ok<INUSD>
            error('Not supported');
        end

        function coeff = computeProjectionCoefficients(proj, pset, varargin) 
            NbyTAbyAttr = pset.buildNbyTAbyConditionAttributes();
            if isempty(proj.K)
                proj.K = 10;
            end
            
            % generate combinations to add together
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            dimIdx = find(dimMask);
            
            % merge each covariate with each covariate + time mixture
            combinedParams = TrialDataUtilities.DPCA.dpca_generateTimeCombinedParams(dimIdx, ...
                'combine', proj.axisCombinations, 'combineEachWithTime', true);
            coeff = TrialDataUtilities.DPCA.dpca(NbyTAbyAttr, proj.K, 'combinedParams', combinedParams, 'lambda', 1e-6);
            
            if proj.makeOrthogonal
                coeff = orth(coeff);
            end
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
