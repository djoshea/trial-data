classdef ProjRegressionPTStoPTS < StateSpaceProjection

    properties
        regressMeanSubtractedData
        target
        basisNameProjStem
    end
    
    methods
        function proj = ProjRegressionPTStoPTS(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            if proj.regressMeanSubtractedData
                pset = pset.meanSubtractBases();
            end
        end
    end
    
    methods(Static)
        function [proj, unmatched] = parseCreateParams(proj, varargin)
            p = inputParser;
            p.addParameter('regressMeanSubtractedData', true, @islogical);
            p.addParameter('target', [], @(x) isa(x, 'PopulationTrajectorySet'));
            p.addParameter('basisNameProjStem', 'pre dict_', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unmatched = p.Unmatched;
            
            proj.regressMeanSubtractedData = p.Results.regressMeanSubtractedData;
            proj.target = p.Results.target;
            proj.basisNameProjStem = p.Results.basisNameProjStem;
        end
        
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjRegressionPTStoPTS();
            [proj, unmatched] = ProjRegressionPTStoPTS.parseCreateParams(proj, varargin{:});
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjRegressionPTStoPTS();
            [proj, unmatched] = ProjRegressionPTStoPTS.parseCreateParams(proj, varargin{:});
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, 'computeStatistics', nargout >= 3, unmatched);
        end
        
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, source, varargin)
            p = inputParser;
            p.addParameter('nBasesProj', NaN, @isscalar); % this is optional for both the caller to provide and the projection to obey
            p.parse(varargin{:});
            
            target = proj.target; %#ok<*PROPLC>
            nBasesProj = p.Results.nBasesProj;
            if isnan(nBasesProj)
                nBasesProj = target.nBases;
            else
                assert(nBasesProj <= target.nBases, 'nBasesProj must be <= target.nBasesProj');
            end
            
            if proj.regressMeanSubtractedData
                targetMeans = target.computeMeanByBasis();
                target = target.meanSubtractBases();
                
            end
            
            assert(isequal(source.nTimeDataMean, target.nTimeDataMean), 'Time point counts do not match between source and target pts');
            assert(source.nConditions == target.nConditions, 'Condition counts do not match between source and target pts');
               
            % extract source data
            source_CTAbyNvalid = source.arrangeCTAbyN('validBasesOnly', true);
            % we can only use timepoints and conditions where all valid bases have data
            source_ctaKeep = all(~isnan(source_CTAbyNvalid), 2);
            
            % extract target data
            target_CTAbyNvalid = target.arrangeCTAbyN('validBasesOnly', true);
            target_CTAbyNvalid = target_CTAbyNvalid(:, 1:nBasesProj);
            
            % we can only use timepoints and conditions where all valid bases have data
            target_ctaKeep = all(~isnan(target_CTAbyNvalid), 2);
            
            ctaKeep = source_ctaKeep & target_ctaKeep;
            if ~any(ctaKeep)
                error('No valid timepoints simultaneously in source and target pts');
            end
            source_CTAbyNvalid = source_CTAbyNvalid(ctaKeep, :);
            target_CTAbyNvalid = target_CTAbyNvalid(ctaKeep, :);
            
            % want decoder to best reconstruct target_CTAbyN = source_CTAbyN * decoder_KbyN'
            decoderNSvalidbyNTvalid = source_CTAbyNvalid \ target_CTAbyNvalid;
            
            % make coefficients for invalid bases 0 so that multiply
            % suppresses invalid bases automatically
            % if NaN is used, they will need to be masked out when
            % multiplying or you'll get a NaN result
            decoderNTbyNS = TensorUtils.inflateMaskedTensor(decoderNSvalidbyNTvalid', [1 2], {target.basisValid, source.basisValid}, 0);
            
            decoderKbyN = decoderNTbyNS;
            encoderNbyK = decoderKbyN';
            
            if ~isempty(target.translationNormalization)
                proj.translationNormalizationPostProject = StateSpaceTranslationNormalization.buildManual(targetMeans, ones(nBasesProj, 1));
            end
        end

        function [nameStem, names] = generateBasisNameProjStem(proj, pset) %#ok<INUSD>
            names = {};
            nameStem = 'regression';
        end
    end

end
