classdef ProjPCA < StateSpaceProjection

    properties
        K % if empty, keep all components. Otherwise, keep only first K components.
    end

    methods
        function proj = ProjPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
            proj.K = [];
        end
        
        function pset = preparePsetForInference(proj, pset) 
            pset = pset.meanSubtractBases('conditionIdx', proj.buildFromConditionIdx);
        end
    end

    methods
        
        function coeff = computeProjectionCoefficients(proj, pset, varargin)
            % run pca on valid bases
            CTAbyNvalid = pset.buildCTAbyN('conditionIdx', proj.buildFromConditionIdx, ...
                'validBasesOnly', true);
            
            idx = find(all(isnan(CTAbyNvalid), 1));
            if ~isempty(idx)
                error('No valid trial average timepoints found for %d bases', numel(idx));
            end
            
            if exist('pca', 'file') == 2
                coeffValid = pca(CTAbyNvalid, 'Rows', 'complete');
            else
                coeffValid = princomp(CTAbyNvalid);
            end
            
            % filter down to K output bases
            if ~isempty(proj.K)
                coeffValid = coeffValid(:, 1:proj.K);
            end
            
            % make coefficients for invalid bases 0 so that multiply
            % suppresses invalid bases automatically
            % if NaN is used, they will need to be masked out when
            % multiplying or you'll get a NaN result
            coeff = zeros(pset.nBases, size(coeffValid, 2));
            coeff(pset.basisValid, :) = coeffValid;
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('PC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
