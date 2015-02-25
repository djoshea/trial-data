classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
        lambda = 1e-6;
        
        makeOrthogonal = false;
    end
    
    methods
        function proj = ProjDPCA_NonOrthogonal(varargin)
            proj = proj@StateSpaceProjection();
        end
        
        function setBuildFromConditionIdx(proj, varargin) %#ok<INUSD>
            error('Not supported');
        end

        function coeff = computeProjectionCoefficients(proj, pset, varargin) 
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            if isempty(proj.K)
                proj.K = min(10, pset.nBasesValid);
            end
            
            % generate combinations to add together
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            dimIdx = find(dimMask);
            
            % merge each covariate with each covariate + time mixture
            combinedParams = TrialDataUtilities.DPCA.dpca_generateTimeCombinedParams(dimIdx, ...
                'combine', proj.axisCombinations, 'combineEachWithTime', true);
            coeffValid = TrialDataUtilities.DPCA.dpca(NvbyTAbyAttr, proj.K, 'combinedParams', combinedParams, 'lambda', 1e-6);
            
            if proj.makeOrthogonal
                coeffValid = orth(coeffValid);
            end
            
            coeff = TensorUtils.inflateMaskedTensor(coeffValid, 1, pset.basisValid);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
