classdef ProjPCA < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
    end

    methods
        function proj = ProjPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
            proj.K = [];
        end
    end

    methods
        
        function coeff = computeProjectionCoefficients(proj, pset, varargin)
            CTAbyN = pset.buildCTAbyN('conditionIdx', proj.buildFromConditionIdx);
            
            if exist('pca', 'file') == 2
               coeff = pca(CTAbyN, 'Rows', 'complete');
            else
                coeff = princomp(CTAbyN);
            end
            
            if ~isempty(proj.K)
                coeff = coeff(:, 1:proj.K);
            end
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('PC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
