classdef ProjPCA < StateSpaceProjection

    methods
        function proj = ProjPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            pset = pset.meanSubtractBases();
        end
    end
    
    methods(Static)
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjPCA();
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjPCA();
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, varargin{:});
        end
        
        function [denoised, proj, stats] = denoiseViaLowRankApproximation(pset, K, varargin)
            if nargin < 2
                K = 10;
            end
            [proj, stats] = ProjPCA.createFrom(pset, 'nBasesProj', K, 'computeStatistics', nargout >= 3);
            denoised = proj.projectInAndOut(pset);
        end
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
            p = inputParser;
            p.addParamValue('nBasesProj', NaN, @isscalar);
            p.parse(varargin{:});
            K = p.Results.nBasesProj;
            
            % run pca on valid bases
            CTAbyNvalid = pset.buildCTAbyN('validBasesOnly', true);
            
            idx = find(all(isnan(CTAbyNvalid), 1));
            if ~isempty(idx)
                error('No valid trial average timepoints found for %d bases', numel(idx));
            end
            
            ctaKeep = ~any(isnan(CTAbyNvalid), 2);
            CTAbyNvalid = CTAbyNvalid(ctaKeep, :);
            CTAbyNvalid = bsxfun(@minus, CTAbyNvalid, mean(CTAbyNvalid, 1));
            
            if exist('pca', 'file') == 2
                [coeffValid] = pca(CTAbyNvalid, 'Rows', 'complete');
            else
                [coeffValid] = princomp(CTAbyNvalid);
            end
            
            % filter down to K output bases, unless K is too small
            if ~isnan(K) && size(coeffValid, 2) > K
                coeffValid = coeffValid(:, 1:K);
            end
            
            % make coefficients for invalid bases 0 so that multiply
            % suppresses invalid bases automatically
            % if NaN is used, they will need to be masked out when
            % multiplying or you'll get a NaN result
            coeff = zeros(pset.nBases, size(coeffValid, 2));
            coeff(pset.basisValid, :) = coeffValid;
            
            % now coeff is N by K, to make the decoderKbyN, we take the transpose
            % the encoder is simply the transpose of the decoder for PCA
            decoderKbyN = coeff';
            encoderNbyK = coeff;
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('PC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end
