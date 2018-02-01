classdef ProjJPCA < StateSpaceProjection

    properties
        subtractCrossConditionMean
        normalize
        softNormOffset

        suppressBWrosettes    % if present and true, the black & white rosettes are not plotted
        suppressHistograms    % if present and true, the blue histograms are not plotted
        suppressText          % if present and true, no text is output to the command window 
    end
    
    methods
        function proj = ProjJPCA(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            if sum(pset.nTimeDataMean) == 1
                error('jPCA not supported for short time windows');
            end
                
            pset = pset.meanSubtractBases(); % this is subtracting the global mean off of each basis, not the cross condition mean
            if proj.normalize
                pset = pset.normalizeBasesSoftRange(proj.softNormOffset);
            end
        end
    end
    
    methods(Static)
        function [proj, unmatched] = parseCreateParams(proj, varargin)
            p = inputParser;
            p.addParameter('subtractCrossConditionMean', false, @islogical);
            p.addParameter('normalize', true, @islogical);
            p.addParameter('softNormOffset', 10, @isscalar); % 10 is a good value here
            
            p.addParameter('suppressBWrosettes', true, @islogical);
            p.addParameter('suppressHistograms', true, @islogical);
            p.addParameter('suppressText', true, @islogical);
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unmatched = p.Unmatched;
            
            assert(~p.Results.subtractCrossConditionMean, 'Not currently supported');

            proj.subtractCrossConditionMean = p.Results.subtractCrossConditionMean;
            proj.normalize = p.Results.normalize;
            proj.softNormOffset = p.Results.softNormOffset;
            proj.suppressBWrosettes = p.Results.suppressBWrosettes;
            proj.suppressHistograms = p.Results.suppressHistograms;
            proj.suppressText = p.Results.suppressText;
        end
        
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjJPCA();
            [proj, unmatched] = ProjJPCA.parseCreateParams(proj, varargin{:});
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjJPCA();
            [proj, unmatched] = ProjJPCA.parseCreateParams(proj, varargin{:});
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, 'computeStatistics', nargout >= 3, unmatched);
        end
        
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
            p = inputParser;
            p.addParameter('nBasesProj', NaN, @isscalar);
            p.parse(varargin{:});
            K = p.Results.nBasesProj;
            if isnan(K)
                K = 6;
            end
            
            % run pca on valid bases
            NvalidbyCbyTA = pset.arrangeNbyCbyTA('validBasesOnly', true);
            
            idx = find(all(isnan(NvalidbyCbyTA), 1));
            if ~isempty(idx)
                error('No valid trial average timepoints found for %d bases', numel(idx));
            end
            
            taKeep = TensorUtils.allMultiDim(~isnan(NvalidbyCbyTA), [1 2]);
            NvalidbyCbyTA = NvalidbyCbyTA(:, :, taKeep);
            NvalidbyCbyTA = bsxfun(@minus, NvalidbyCbyTA, mean(NvalidbyCbyTA, 3));
            
            decoderKbyNvalid = TrialDataUtilities.jPCA.jPCA(NvalidbyCbyTA, K, ...
                'suppressBWrosettes', proj.suppressBWrosettes, ...
                'suppressHistograms', proj.suppressHistograms, ...
                'suppressText', proj.suppressText);

            % make coefficients for invalid bases 0 so that multiply
            % suppresses invalid bases automatically
            % if NaN is used, they will need to be masked out when
            % multiplying or you'll get a NaN result
            
            % now coeff is N by K, to make the decoderKbyN, we take the transpose
            % the encoder is simply the transpose of the decoder for PCA
            decoderKbyN = TensorUtils.inflateMaskedTensor(decoderKbyNvalid, 2, pset.basisValid);
            encoderNbyK = decoderKbyN';
        end

        function [nameStem, names] = generateBasisNameProjStem(proj, pset) %#ok<INUSD>
            names = {};
            nameStem = 'jPC';
        end
    end

end
