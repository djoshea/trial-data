classdef ProjDPCA < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
    end

    methods
        function proj = ProjDPCA(varargin)
            proj = proj@StateSpaceProjection();
        end
        
        function setBuildFromConditionIdx(proj, varargin) %#ok<INUSD>
            error('Not supported');
        end

        function coeff = computeProjectionCoefficients(proj, pset, varargin) 
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            coeffValid = dpca_nanSafe(NvbyTAbyAttr, proj.K, 10000, 1e-7);
            coeff = TensorUtils.inflateMaskedTensor(coeffValid, 1, pset.basisValid);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
