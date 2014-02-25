classdef ProjDPCA < StateSpaceProjection

    properties
        P = 10;
    end

    methods
        function proj = ProjDPCA(varargin)
            proj = proj@StateSpaceProjection();
        end
        
        function setBuildFromConditionIdx(proj, varargin) %#ok<INUSD>
            error('Not supported');
        end
        
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

        function coeff = calculateProjectionCoefficients(proj, pset, varargin) 
            NbyTAbyAttr = pset.buildNbyTAbyConditionAttributes();
            coeff = dpca_nanSafe(NbyTAbyAttr, proj.P, [], []);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
