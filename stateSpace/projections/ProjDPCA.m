classdef ProjDPCA < StateSpaceProjection

        
    methods(Static)
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjDPCA();
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjDPCA();
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, varargin{:});
        end
    end
    
    methods
        function proj = ProjDPCA(varargin)
            proj = proj@StateSpaceProjection();
        end

        function coeff = computeProjectionCoefficients(proj, pset, varargin) 
            p = inputParser;
            p.addParamValue('nBasesProj', pset.nBasesValid, @isscalar);
            p.parse(varargin{:});
            K = p.Results.nBasesProj;
            
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            coeffValid = dpca_nanSafe(NvbyTAbyAttr, K, 10000, 1e-7);
            coeff = TensorUtils.inflateMaskedTensor(coeffValid, 1, pset.basisValid);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
