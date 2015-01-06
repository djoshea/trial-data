classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
        lambda = 1e-6;
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
                K = 10;
            else
                K = proj.K;
            end
            
            % generate combinations to add together
            combinedParams = dpca_generateTimeCombinedParams(pset.conditionDescriptor.nAxes);
            coeff = dpca(NbyTAbyAttr, K, 'combinedParams', combinedParams, 'lambda', 1e-6);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
