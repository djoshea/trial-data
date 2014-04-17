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
            NbyTAbyAttr = pset.buildNbyTAbyConditionAttributes();
            coeff = dpca_nanSafe(NbyTAbyAttr, proj.K, 10000, 1e-7);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
