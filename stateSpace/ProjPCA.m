classdef ProjPCA < StateSpaceProjection

    properties
        tsquared
        
        P = [];
    end

    methods
        function proj = ProjPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
    end

    methods
        function coeff = calculateBasisCoefficients(proj, pset, ctaByN, varargin)
            if exist('pca') == 2
                [coeff] = pca(ctaByN, 'Rows', 'complete');
            else
                [coeff] = princomp(ctaByN);
            end
            
            if ~isempty(proj.P)
                coeff = coeff(:, 1:proj.P);
            end
        end

        function names = getBasisNames(proj, pset, ctaByN)
            names = arrayfun(@(i) sprintf('PC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
