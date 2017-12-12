classdef ProjPCA < StateSpaceProjection

    properties
        meanSubtract = true;
    end
    
    methods
        function proj = ProjPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            if proj.meanSubtract
                if sum(pset.nTimeDataMean) == 1
                    warning('Doing mean subtraction on single timepoint data. Set .meanSubtract = false');
                end
                    
                pset = pset.meanSubtractBases();
            end
        end
    end
    
    methods(Static)
        function [proj, unmatched] = parseCreateParams(proj, varargin)
            p = inputParser;
            p.addParameter('meanSubtract', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unmatched = p.Unmatched;
            
            proj.meanSubtract = p.Results.meanSubtract;
        end
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjPCA();
            [proj, unmatched] = ProjPCA.parseCreateParams(proj, varargin{:});
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjPCA();
            [proj, unmatched] = ProjPCA.parseCreateParams(proj, varargin{:});
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, 'computeStatistics', nargout >= 3, unmatched);
        end
        
        function [denoised, proj, stats] = denoiseViaLowRankApproximation(pset, K, varargin)
            if nargin < 2
                K = 10;
            end
            [proj, stats] = ProjPCA.createFrom(pset, 'nBasesProj', K, 'computeStatistics', nargout >= 3);
            denoised = proj.projectInAndOut(pset);
        end
        
%         function [proj, psetProjected, stats] = createFromAndProjectCaptureThresholdVariance(pset, threshFractionVar, varargin)
%             % truncate the output bases so as to capture a certain threshFractionVar of the
%             % variance. Look at proj.nBasesProj to get this number
%             p = inputParser;
%             p.addParamValue('threshForSignalVariance', false, @isscalar); % set threshhold for signal variance. if false, for overall variance
%             p.addParamValue('threshForMarginalizationVariance', [], @(x) isempty(x) || StateSpaceProjectionStatistics.isMarginalizationSpec(x)); 
%             p.parse(varargin{:});
%             
%             [proj, psetProjected, stats] = ProjPCA.createFromAndProject(pset, varargin{:});
%             proj.truncate
%         end
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
            p = inputParser;
            p.addParameter('nBasesProj', NaN, @isscalar);
            p.parse(varargin{:});
            K = p.Results.nBasesProj;
            
            % run pca on valid bases
            CTAbyNvalid = pset.arrangeCTAbyN('validBasesOnly', true);
            
            idx = find(all(isnan(CTAbyNvalid), 1));
            if ~isempty(idx)
                error('No valid trial average timepoints found for %d bases', numel(idx));
            end
            
            ctaKeep = ~any(isnan(CTAbyNvalid), 2);
            CTAbyNvalid = CTAbyNvalid(ctaKeep, :);
            if proj.meanSubtract
                CTAbyNvalid = bsxfun(@minus, CTAbyNvalid, mean(CTAbyNvalid, 1));
            end
            
            if size(CTAbyNvalid, 1) == 1
                % single timepoint, just normalize as the direction
                coeffValid = CTAbyNvalid' / sqrt(sum(CTAbyNvalid(:).^2));
            elseif exist('pca', 'file') == 2
                [coeffValid] = pca(CTAbyNvalid, 'Rows', 'complete');
            else
                [coeffValid] = princomp(CTAbyNvalid); %#ok<PRINCOMP>
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

        function [nameStem, names] = generateBasisNameProjStem(proj, pset) %#ok<INUSD>
            names = {};
            nameStem = 'PC';
        end
    end

end
