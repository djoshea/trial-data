classdef ProjDPCA < StateSpaceProjection

    properties
        P = 10;
    end

    methods
        function proj = ProjDPCA(varargin)
            proj = proj@StateSpaceProjection();
        end
        
        function coeff = calculateBasisCoefficients(proj, pset, ctaByN, varargin)
            proj.useCommonValidTimeWindow = true;
            
            dataTensor = proj.extractDataForMarginalization(pset);
            dataTensor = proj.prepareBases(dataTensor, 1);
            
            coeff = dpca_nanSafe(dataTensor, proj.P, [], []);
        end

        function names = getBasisNames(proj, pset, ctaByN)
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
