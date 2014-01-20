classdef StateSpaceProjection 
% Abstract base class representing a set of coefficients per basis with which
% to project PopulationTrajectorySet instances. Typically these are built using 
% fromPopulationTrajectorySet to compute the coefficients, in a manner determined
% by subclasses

    properties
        translationNormalization
        coeff % N x K matrix of component weights; coeff(i,j) is component j from neuron i 
    end

    properties(Dependent)
        nBasesSource
        nBasesProj
    end
    
    methods(Abstract)
        names = getBasisNames(proj, pset, data)
    end

    methods
        function psetProjected = projectPopulationTrajectorySet(proj, pset, varargin)
            assert(pset.nBases == proj.nBasesSource, 'Number of bases must match in order to project');

            psetProjected = PopulationTrajectorySet('alignDescriptorSet', pset.alignDescriptorSet, 'conditionDescriptor', pset.conditionDescriptor);
            
            % WARNING: THIS IS A TEMPORARY HACK
            psetProjected.alignTimeInfoData = pset.alignTimeInfoData(1:proj.nBasesProj, :, :);
            
            psetProjected.dataValid = true(proj.nBasesProj, pset.nConditions, pset.nAlign);
            psetProjected.nTrialsData = nan(proj.nBasesProj, pset.nConditions, pset.nAlign);
            psetProjected.tMinDataManual = nan(proj.nBasesProj, pset.nConditions, pset.nAlign);
            psetProjected.tMaxDataManual = nan(proj.nBasesProj, pset.nConditions, pset.nAlign);
            
            [tByNEachCA, tvecByConditionAlign] = proj.extractDataForProjection(pset);
           
            [projData projTimeData] = deal(cell(proj.nBasesProj, pset.nConditions, pset.nAlign));
            for iAlign = 1:pset.nAlign
                for iCondition = 1:pset.nConditions
                    tByN = tByNEachCA{iCondition, iAlign};
                    if isempty(tByN)
                        continue;
                    end
                    tByN = proj.prepareBases(tByN, 2);
                    timeVec = tvecByConditionAlign{iCondition, iAlign};

                    projCond = tByN * proj.coeff;  
                    projData(:, iCondition, iAlign) = mat2cell(projCond, ...
                        size(tByN, 1), ones(proj.nBasesProj, 1));
                    projTimeData(:, iCondition, iAlign) = repmat({timeVec}, proj.nBasesProj, 1);
                end
                
                psetProjected.data = projData;
                psetProjected.timeData = projTimeData;
                psetProjected.basisMeta = proj.basisMeta;
                psetProjected.basisNames = proj.basisNames;
            end
        end
        
        function [tByNEachCA, tvecByConditionAlign] = extractDataForProjection(proj, pset)
            [tByNEachCA, tvecByConditionAlign] = pset.buildTByNEachCA('timeValidAcrossConditions', proj.useCommonValidTimeWindow);
        end 
        
        function ctaByN = extractDataForSelection(proj, pset)
            % selectively take only common valid time windows if requested
            ctaByN = pset.buildCTAByN('timeValidAcrossConditions', proj.useCommonValidTimeWindow);
        end
        
        function dataTensor = extractDataForMarginalization(proj, pset)
            proj.assertCommonValidTime();
            dataTensor = pset.buildNByTAByAttributeTensor();
        end
        
        function fromPopulationTrajectorySet(proj, pset, varargin)
            % build this projection matrix based on an existing PopulationTrajectorySet
            % defers to calculateProjectionMatrix for the actual basis computation
            p = inputParser;
            p.addRequired('pset', @(x) isa(x, 'PopulationTrajectorySet'));
            p.parse(pset, varargin{:});

            results = StateSpaceProjectionResults();

            % make any necessary transformations, particularly translation / normalization
            pset = proj.preparePsetForInference(pset);

            % extract the translation normalization that will be applied before projection
            proj.translationNormalization = pset.translationNormalization;

            % build psth matrix nConditions*nTimepoints*nAlignments x nNeurons
            CTAbyN = pset.buildCTAbyN();

            proj.coeff = proj.calculateBasisCoefficients(pset, ctaByN);
            
            [results.covSource results.corrSource] = proj.calculateSourceBasisCovariance(pset, ctaByN);
            results.latent = proj.calculateLatentVariance();
            
            results.explained = results.latent / trace(results.covSource);
            
            [covFull results.covMarginalized results.basisMixtureNames] = ...
                proj.calculateMarginalizedCovariances(pset);
            [results.latentMarginalized results.basisMixtures] = ...
                proj.calculateLatentVarianceMarginalized();
        end

    end

    methods
        function n = get.nBasesSource(proj)
            n = size(proj.coeff, 1);
        end

        function n = get.nBasesProj(proj)
            n = size(proj.coeff, 2);
        end

        function offsets = preparePsetForInference(proj, pset)
            pset = pset.meanSubtractBases();
            nBases = size(ctaByN, 2);
            switch proj.offsetMode
                case ''
                    offsets = zeros(nBases, 1);
                case 'meanSubtract'
                    offsets = nanmean(ctaByN, 1);
                otherwise
                    error('Unknown offsetMode %s', proj.offsetMode);
            end
            
            offsets = makecol(offsets);
        end

        function normalization = calculateBasisNormalization(proj, pset, ctaByN)
            switch proj.normalizationMode
                case 'none'
                    normalization = ones(1, size(ctaByN, 2));
                    
                case 'softMax'
                    normalization = max(abs(ctaByN), [], 1) + proj.softMaxAlpha;
               
                case 'softMaxCrossConditionVariance'
                    ccvByUnit = pset.crossConditionVariance();
                    maxCCV = cellfun(@(ccv) nanmax(ccv), ccvByUnit);
                    normalization = maxCCV' + proj.softMaxAlpha;

                case 'softMaxCrossConditionStd'
                    ccsByUnit = pset.crossConditionStd();
                    maxCCS = cellfun(@(ccs) nanmax(ccs), ccsByUnit);
                    normalization = maxCCS' + proj.softMaxAlpha;
                    
                case 'softMaxVariance'
                    normalization = var(ctaByN, [], 1) + proj.softMaxAlpha;
                    
                case 'softMaxStd'
                    normalization = std(ctaByN, [], 1) + proj.softMaxAlpha;
            end
            
            normalization = makecol(normalization);
        end

        function [covariance corr] = calculateSourceBasisCovariance(proj, pset, ctaByN)
            covariance = nancov(ctaByN, 1); 
            corr = corrcoef(ctaByN, 'rows', 'complete');
        end

        function data = prepareBases(proj, data, basesInDim)
            offsets = shiftdim(makecol(proj.basisOffsetsSource), -(basesInDim-1));
            normalization = shiftdim(makecol(proj.basisNormalizationSource), -(basesInDim-1));
            data = bsxfun(@minus, data, offsets);
            data = bsxfun(@rdivide, data, normalization);
        end

        function latent = calculateLatentVariance(proj)
            % covSource is N*N, coeff is N*K, latent is K*K
            latent = diag(proj.coeff' * proj.covSource * proj.coeff);
        end
        
        function assertCommonValidTime(proj)
            assert(proj.useCommonValidTimeWindow, 'Must build using .useCommonValidTimeWindow = true for marginalization');
        end
            
        function [latentMarginalized basisMixtures] = calculateLatentVarianceMarginalized(proj)
            proj.assertCommonValidTime();
            
            nMarginalized = length(proj.covMarginalized);
            latentMarginalized = nan(proj.nBasesProj, nMarginalized);
            for iCov = 1:length(proj.covMarginalized)
                latentMarginalized(:, iCov) = diag(proj.coeff' * proj.covMarginalized{iCov} * proj.coeff);
            end
            
            basisMixtures = bsxfun(@rdivide, latentMarginalized, sum(latentMarginalized, 2));
        end
        
        function [covFull covMarginalized covMarginalizedNames] = calculateMarginalizedCovariances(proj, pset)
            % Use dpca_covs to compute each covariance matrix for us
            proj.assertCommonValidTime();
            
            dataTensor = proj.extractDataForMarginalization(pset);
            dataTensor = proj.prepareBases(dataTensor, 1);
            
            attrNames = [ {'time'}, pset.conditionDescriptor.groupByList{:} ];
            [covMarginalizedMap covFull attrSets] = dpca_covs_nanSafe(dataTensor);

            nCov = length(covMarginalizedMap);
            covMarginalized = cell(nCov, 1);
            for iCov = 1:nCov
                covMarginalized{iCov} = covMarginalizedMap(mat2str(attrSets{iCov}));
            end
            
            % generate mixture names
            covMarginalizedNames = cellfun(@(attrInds) strjoin(attrNames(attrInds-1), ' x '), ...
                attrSets, 'UniformOutput', false);
            covMarginalizedNames = makecol(covMarginalizedNames);
        end
        
        function plotBasisMixtures(proj, varargin)
            p = inputParser;
            p.addParamValue('basisIdx', [1:min([10 proj.nBasesProj])], @(x) isvector(x) && ...
                all(inRange(x, [1 proj.nBasesSource])));
            p.parse(varargin{:});
            basisIdx = p.Results.basisIdx;
            if islogical(basisIdx)
                basisIdx = find(basisIdx);
            end
            
            proj.assertCommonValidTime();
            cla;
%             p = panel();
%             p.pack(1,1);
%             p(1,1).select();
            
            cumBasisMix = cumsum(proj.basisMixtures(basisIdx,:), 2);
            nCov = size(proj.basisMixtures, 2);
            %nBases = proj.nBasesProj;
            nBases = length(basisIdx);
            rowHeight = 0.8;
            xMin = 0;
            xMax = 1;
            hPatch = nan(nCov, 1);
            for iIdxB = 1:nBases
                iB = basisIdx(iIdxB);
                for iCov = 1:nCov
                    % SW, NW, NE, SE order for rectangle
                    if iCov == 1
                        patchX = [0 0 cumBasisMix(iB, 1) cumBasisMix(iB, 1)];
                    else 
                        patchX = [cumBasisMix(iB, iCov-1) cumBasisMix(iB, iCov-1) ...
                            cumBasisMix(iB, iCov) cumBasisMix(iB, iCov)];
                    end
                    patchY = [nBases-iB nBases-iB+rowHeight nBases-iB+rowHeight nBases-iB];
                    
                    hPatch(iCov) = patch(patchX, patchY, proj.basisMixtureColors(iCov, :)); %, 'EdgeColor', 'none');
                    hold on
                end
                
                h = text(-0.03, nBases-iB+rowHeight/2, proj.basisNames{iB}, ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Right'); 
                extent = get(h, 'Extent');
                xMin = min(extent(1), xMin);
                
                h = text(1.03, nBases-iB+rowHeight/2, sprintf('%.2f%%', proj.explained(iB)*100), ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Left'); 
                extent = get(h, 'Extent');
                xMax = max(extent(1)+extent(3), xMax);
            end
            
            xlim([xMin xMax]);
            ylim([0 nBases]);
            axis off;
            title('Basis Mixtures');
            legend(hPatch, nCov, proj.basisMixtureNames, 'Location', 'NorthEastOutside');
            legend boxoff;
            
        end
        
    end

end
