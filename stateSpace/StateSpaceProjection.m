classdef StateSpaceProjection 
% Abstract base class representing a set of coefficients per basis with which
% to project PopulationTrajectorySet instances. Typically these are built using 
% fromPopulationTrajectorySet to compute the coefficients, in a manner determined
% by subclasses

    % settings which control how the StateSpaceProjection behaves
    % These must be set before building from PopulationTrajectorySet
    properties
        % if true, computes the dataErrorHigh/Low for the projected bases
        % using resampling methods and the interval size from the pset
        propagateDataError = false; 
    end

    properties(SetAccess=protected)
        buildFromConditionIdx % if non-empty, build from selected conditions only
        initialized = false;
        translationNormalization % stored translation / normalization stored when building and used when projecting
        coeff % N x K matrix of component weights; coeff(i,j) is weight for component j from input basis i 
    end

    properties(Dependent)
        nBasesSource
        nBasesProj
    end

    % Simple dependent property getters
    methods
        function n = get.nBasesSource(proj)
            if isempty(proj.coeff)
                n = NaN;
            else
                n = size(proj.coeff, 1);
            end
        end

        function n = get.nBasesProj(proj)
            if isempty(proj.coeff)
                n = NaN;
            else
                n = size(proj.coeff, 2);
            end
        end
        
        % we implement this as a method so that subclasses can disable this
        % functionality if they don't support it
        function proj = setBuildFromConditionIdx(proj, idx)
            proj.warnIfNoArgOut(nargout);
            proj.buildFromConditionIdx = idx;
        end   
    end
    
    methods(Abstract)
        % return a list of basis names for the new basis
        names = getBasisNames(proj, pset, data)

        % compute the N * K matrix of basis coefficients for the projection
        coeff = computeProjectionCoefficients(pset)
    end

    methods
        function [proj, stats] = buildFromPopulationTrajectorySet(proj, pset, varargin)
            % build this projection matrix based on an existing PopulationTrajectorySet
            % defers to calculateProjectionMatrix for the actual basis computation
            
            proj.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('pset', @(x) isa(x, 'PopulationTrajectorySet'));
            p.parse(pset, varargin{:});

            if isempty(proj.buildFromConditionIdx)
                proj.buildFromConditionIdx = truevec(pset.nConditions);
            end
            
            % make any necessary transformations, particularly translation / normalization
            pset = proj.preparePsetForInference(pset);

            % extract the translation normalization that will be applied before projection
            proj.translationNormalization = pset.translationNormalization;

            % compute the coefficients for the projection
            proj.coeff = proj.computeProjectionCoefficients(pset);
            
            % results will store statistics and useful quantities related to the
            % projection
            stats = proj.computeProjectionStatistics(pset, true);

            proj.initialized = true;
        end

        function [psetProjected, stats] = projectPopulationTrajectorySet(proj, pset)
            assert(pset.nBases == proj.nBasesSource, ...
                'Number of bases must match in order to project');

            % replace translation normalization
            pset = pset.clearTranslationNormalization().applyTranslationNormalization(proj.translationNormalization);
            
            % copy basic settings from pset 
            b = PopulationTrajectorySetBuilder.copySettingsDescriptorsFromPopulationTrajectorySet(pset);

            b.basisNames = proj.getBasisNames(pset);
            b.basisUnits = proj.getBasisUnits(pset); 
            
            % copy/compute trial averaged data 
            
            b.tMinForDataMean = pset.tMinForDataMean;
            b.tMaxForDataMean = pset.tMaxForDataMean;

            % sum across bases for projection (dataNTrials is nAlign x nBases x nConditions)
            b.dataNTrials = repmat(sum(pset.dataNTrials, 2), [1, proj.nBasesProj, 1]);
            % all across bases for validity
            b.dataValid = repmat(all(pset.dataValid, 2), [1, proj.nBasesProj, 1]);
        
            % project dataMean, leave dataSem as NaN
            [b.dataMean, b.dataSem] = deal(cell(pset.nAlign, 1));
            for iAlign = 1:pset.nAlign

                % mat is N x CT, coeff is N x K
                mat = reshape(pset.dataMean{iAlign}, pset.nBases, ...
                    pset.nConditions * pset.nTimeDataMean(iAlign));
                projMat = proj.coeff' * mat;
                b.dataMean{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));

                % leave sem as NaN
                b.dataSem{iAlign} = nan(proj.nBasesProj, pset.nConditions, pset.nTimeDataMean(iAlign));
            end
          
            % project randomized data, recompute intervals
            if proj.propagateDataError
                [b.dataMeanRandomized, b.dataIntervalLow, b.dataIntervalHigh] = deal(cell(pset.nAlign, 1));
                for iAlign = 1:pset.nAlign
                    % mat is N x CTS, coeff is N x K, where S is nRandomSamples 
                    mat = reshape(pset.dataMeanRandomized{iAlign}, pset.nBases, ...
                        pset.nConditions*pset.nTimeDataMean(iAlign)*pset.nRandomSamples);
                    projMat = proj.coeff' * mat;
                    b.dataMeanRandomized{iAlign} = reshape(projMat, proj.nBasesProj, ...
                        pset.nConditions, pset.nTimeDataMean(iAlign), pset.nRandomSamples);

                    quantiles = quantile(b.dataMeanRandomized{iAlign}, ...
                        [pset.dataIntervalQuantileLow, pset.dataIntervalQuantileHigh], 4);
                    b.dataIntervalQuantileLow{iAlign} = quantiles(:, :, :, 1);
                    b.dataIntervalQuantileHigh{iAlign} = quantiles(:, :, :, 2);
                end
            end

            % aggregate AlignSummary data. Each projected basis samples trials from all original
            % trials, so we aggregate all AlignSummary instances into one
            b.alignSummaryData = pset.alignSummaryAggregated';
            b.basisAlignSummaryLookup = ones(pset.nBases, 1);
            
            psetProjected = b.buildManualWithTrialAveragedData();
            stats = proj.computeProjectionStatistics(pset, false);
        end
        
        function [psetProjected, proj, statsBuild, statsProject] = buildFromAndProjectPopulationTrajectorySet(proj, pset, varargin)
            [proj, statsBuild] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
            [psetProjected, statsProject] = proj.projectPopulationTrajectorySet(pset);
        end
    end

    methods(Access=protected, Sealed)
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~ishandle(obj)
                message = sprintf('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
                expr = sprintf('debug(''%s'')', message);
                evalin('caller', expr); 
            end
        end
    end

    methods
        function proj = filterBases(proj, idx)
            % after building from a PopulationTrajectorySet, mask the bases
            % within to reduce the total contribution
            proj.warnIfNoArgOut(nargout);
            assert(proj.initialized, 'Call filterBases after building / initializing');
            
            proj.translationNormalization = proj.translationNormalization.filterBases(idx);
            proj.coeff = proj.coeff(idx, :);
        end
        
        function names = getBasisUnits(proj, pset)  %#ok<INUSL>
            names = repmat({''}, pset.nBases, 1);
        end
    
        function pset = preparePsetForInference(proj, pset) 
            % apply any appropriate translations, normalizations, or other adjustments
            % to pset before inferring coefficients for projection. The .translationNormalization
            % found in pset after this function runs will be used to normalize all psets
            % that are projected via this StateSpaceProjection. By default, this will
            % will not do any normalization, but will perform mean-subtraction. 
            % the caller may manually specify the normalization in the pset before 
            % building the StateSpaceProjection. Subclasses may wish to override this method
            % to prevent mean-subtraction or add basis normalization, if necessary

            pset = pset.meanSubtractBases('conditionIdx', proj.buildFromConditionIdx);
        end
        
        function s = computeProjectionStatistics(proj, pset, forBuild)
            % for buildIndicates if we are computing statistics at the time
            % of projection building, false indicates we are simply
            % projecting a pset with already computed coefficients.
            
            s = StateSpaceProjectionStatistics();
            
            % get the names and counts of source -> proj bases
            s.basisNamesSource = pset.basisNames;
            s.basisNamesProj = proj.getBasisNames(pset);
            s.nBasesSource = proj.nBasesSource;
            s.nBasesProj = proj.nBasesProj;

            if forBuild
                CTAbyN = pset.buildCTAbyN('conditionIdx', proj.buildFromConditionIdx);
            else
                CTAbyN = pset.buildCTAbyN();
            end

            s.covSource = nancov(CTAbyN, 1); 
            s.corrSource= corrcoef(CTAbyN, 'rows', 'complete');

            % covSource is N*N, coeff is N*K, latent is K*K
            s.latent = diag(proj.coeff' * s.covSource * proj.coeff);
            s.explained = s.latent / trace(s.covSource);
            
            % Use dpca_covs_nanSafe to compute each covariance matrix for us
            NbyTAbyAttr = pset.buildNbyTAbyConditionAttributes();
            attrNames = [ {'time'}, pset.conditionDescriptor.axisNames{:} ];
            [covMarginalizedMap, ~, attrSets] = dpca_covs_nanSafe(NbyTAbyAttr);

            % generate mixture names
            covMarginalizedNames = cellfun(@(attrInds) ...
                strjoin(attrNames(attrInds-1), ' x '), ...
                attrSets, 'UniformOutput', false);
            s.covMarginalizedNames = makecol(covMarginalizedNames);

            % unpack covariance map into nMixtures x 1 cell array
            nCov = length(covMarginalizedMap);
            s.covMarginalized = cell(nCov, 1);
            for iCov = 1:nCov
                s.covMarginalized{iCov} = covMarginalizedMap(mat2str(attrSets{iCov}));
            end
            
            % compute the marginalized variance explained, i.e. the amount
            % of source variance marginalized according to dpca_cov, in each direction 
            s.latentMarginalized = nan(proj.nBasesProj, nCov);
            for iCov = 1:nCov
                s.latentMarginalized(:, iCov) = ...
                    diag(proj.coeff' * s.covMarginalized{iCov} * proj.coeff);
            end
        end

%         function normalization = calculateBasisNormalization(proj, pset, ctaByN)
%             switch proj.normalizationMode
%                 case 'none'
%                     normalization = ones(1, size(ctaByN, 2));
%                     
%                 case 'softMax'
%                     normalization = max(abs(ctaByN), [], 1) + proj.softMaxAlpha;
%                
%                 case 'softMaxCrossConditionVariance'
%                     ccvByUnit = pset.crossConditionVariance();
%                     maxCCV = cellfun(@(ccv) nanmax(ccv), ccvByUnit);
%                     normalization = maxCCV' + proj.softMaxAlpha;
% 
%                 case 'softMaxCrossConditionStd'
%                     ccsByUnit = pset.crossConditionStd();
%                     maxCCS = cellfun(@(ccs) nanmax(ccs), ccsByUnit);
%                     normalization = maxCCS' + proj.softMaxAlpha;
%                     
%                 case 'softMaxVariance'
%                     normalization = var(ctaByN, [], 1) + proj.softMaxAlpha;
%                     
%                 case 'softMaxStd'
%                     normalization = std(ctaByN, [], 1) + proj.softMaxAlpha;
%             end
%             
%             normalization = makecol(normalization);
%         end
    end

end
