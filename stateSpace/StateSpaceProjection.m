classdef StateSpaceProjection 
% Abstract base class representing a set of coefficients per basis with which
% to project PopulationTrajectorySet instances. Typically these are built using 
% fromPopulationTrajectorySet to compute the coefficients, in a manner determined
% by subclasses

    % settings which control how the StateSpaceProjection behaves
    % These must be set before building from PopulationTrajectorySet
    properties
        axisCombinations % see TrialDataUtilities.DPCA.dpca for dimListsToCombineList argument
    end

    properties(SetAccess=protected)
        buildFromConditionIdx % if non-empty, build from selected conditions only
        initialized = false;
        translationNormalization % stored translation / normalization stored when building and used when projecting
        coeff % N x K matrix of component weights; coeff(i,j) is weight for component j from input basis i 
        
        basisValid % N x 1 logical vector indicating which bases were considered valid in the projection
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
    
    methods
        function proj = set.axisCombinations(proj, v)
            assert(iscell(v), 'Axis Combinations must be a cell of vectors');
            proj.axisCombinations = v;
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
            
            assert(size(proj.coeff, 1) == pset.nBases, 'Coefficient matrix returned by computeProjectionCoefficients must match pset.nBases along dim 1');
            
            % copy the basis valid mask
            proj.basisValid = pset.basisValid;
            
            % results will store statistics and useful quantities related to the
            % projection
            stats = proj.computeProjectionStatistics(pset, true);

            proj.initialized = true;
        end

        function [psetProjected, stats] = projectPopulationTrajectorySet(proj, pset, applyTranslationNormalization)
            if nargin < 3
                applyTranslationNormalization = true;
            end
            
            assert(pset.nBases == proj.nBasesSource, ...
                'Number of bases must match in order to project');

            % replace translation normalization
            if applyTranslationNormalization
                debug('Applying Translation/Normalization to data\n');
                pset = pset.clearTranslationNormalization().applyTranslationNormalization(proj.translationNormalization);
            end
            
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

                % mat is N x CT, coeff is Nvalid x K
                mat = reshape(pset.dataMean{iAlign}, pset.nBases, ...
                    pset.nConditions * pset.nTimeDataMean(iAlign));

                % coeff must have 0 for invalid bases, but we'll zero out
                % mat because it has NaNs for invalid bases that will mess
                % up the matrix multiply
                mat(~proj.basisValid, :) = 0;

                projMat = proj.coeff' * mat;
                b.dataMean{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));

                % use sqrt(sd1^2 / n1 + sd2^2 / n2 + ...) formula
                % which equals sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                mat = reshape(pset.dataSem{iAlign}, pset.nBases, ...
                    pset.nConditions * pset.nTimeDataMean(iAlign));
                mat(~proj.basisValid, :) = 0;
                projMat = sqrt(abs(proj.coeff)' * (mat.^2));
                b.dataSem{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));
            end
          
            % project randomized data, recompute intervals
            if ~isempty(pset.dataMeanRandomized)
                [b.dataMeanRandomized, b.dataIntervalLow, b.dataIntervalHigh] = deal(cell(pset.nAlign, 1));
                for iAlign = 1:pset.nAlign
                    % mat is N x CTS, coeff is N x K, where S is nRandomSamples 
                    mat = reshape(pset.dataMeanRandomized{iAlign}, pset.nBases, ...
                        pset.nConditions*pset.nTimeDataMean(iAlign)*pset.nRandomSamples);
                    mat(~proj.basisValid, :) = 0;
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
        
        function [proj, psetProjected, statsBuild, statsProject] = buildFromAndProjectPopulationTrajectorySet(proj, pset, varargin)
            [proj, statsBuild] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
            [psetProjected, statsProject] = proj.projectPopulationTrajectorySet(pset, false); % false --> don't apply translation normalization since we'll pull that from this pset to begin with
        end
    end

    methods(Access=protected, Sealed)
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~ishandle(obj)
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
            end
        end
    end

    methods
        function proj = filterBases(proj, idx)
            % select on output bases
            proj.warnIfNoArgOut(nargout);
            assert(proj.initialized, 'Call filterBases after building / initializing');
            proj.coeff = proj.coeff(:, idx); % select on output bases
        end
        
        function names = getBasisUnits(proj, pset)  %#ok<INUSL>
            names = repmat({''}, pset.nBases, 1);
        end
    
        function pset = preparePsetForInference(proj, pset) 
            % apply any appropriate translations, normalizations, or other adjustments
            % to pset before inferring coefficients for projection. The .translationNormalization
            % found in pset after this function runs will be used to normalize all psets
            % that are projected via this StateSpaceProjection. By default, this will
            % will not do anything.
            
            % the caller may manually specify the normalization in the pset before 
            % building the StateSpaceProjection. Subclasses may wish to override this method
            % to do mean-subtraction or add basis normalization, if necessary

           % pset = pset.meanSubtractBases('conditionIdx', proj.buildFromConditionIdx);
        end
        
        function s = computeProjectionStatistics(proj, pset, forBuild)
            % for buildIndicates if we are computing statistics at the time
            % of projection building, false indicates we are simply
            % projecting a pset with already computed coefficients.
            %
            % we'll work with valid bases only and then inflate the results
            % to the full size
            
            import TrialDataUtilities.DPCA.*;
            
            s = StateSpaceProjectionStatistics();
            
            % get the names and counts of source -> proj bases
            s.basisNamesSource = pset.basisNames;
            s.basisNamesProj = proj.getBasisNames(pset);
            s.nBasesSource = proj.nBasesSource;
            s.nBasesProj = proj.nBasesProj;

            if forBuild
                CTAbyNv = pset.buildCTAbyN('validBasesOnly', true, 'conditionIdx', proj.buildFromConditionIdx);
            else
                CTAbyNv = pset.buildCTAbyN('validBasesOnly', true);
            end

            s.covSource = TensorUtils.inflateMaskedTensor(nancov(CTAbyNv, 1), [1 2], proj.basisValid); 
            s.corrSource = TensorUtils.inflateMaskedTensor(corrcoef(CTAbyNv, 'rows', 'complete'), [1 2], proj.basisValid);

            % covSource is N*N, coeff is N*K, latent is K*1
            s.latent = diag(proj.coeff(proj.basisValid, :)' * s.covSource(proj.basisValid, proj.basisValid) * proj.coeff(proj.basisValid, :));
            s.explained = s.latent / trace(s.covSource(proj.basisValid, proj.basisValid));
            
            % Use dpca_covs_nanSafe to compute each covariance matrix for us
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            
            % filter for non-singular axes
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            dimIdx = find(dimMask);
            
            % merge each covariate with each covariate + time mixture
            [combinedParams, s.covMarginalizedNames] = dpca_generateTimeCombinedParams(...
                dimIdx, 'dimNames', pset.conditionDescriptor.axisNames(:), 'combineEachWithTime', true);
            covMarginalizedValidCell = dpca_marginalizedCov(NvbyTAbyAttr, 'combinedParams', combinedParams);
            s.covMarginalized = cellfun(@(x) TensorUtils.inflateMaskedTensor(x, [1 2], proj.basisValid), ...
                covMarginalizedValidCell, 'UniformOutput', false);

            % compute the marginalized variance explained, i.e. the amount
            % of source variance marginalized according to dpca_cov, in each direction 
            nCov = numel(s.covMarginalized);
            s.latentMarginalized = nan(proj.nBasesProj, nCov);
            for iCov = 1:nCov
                s.latentMarginalized(:, iCov) = ...
                    diag(proj.coeff(proj.basisValid, :)' * s.covMarginalized{iCov}(proj.basisValid, proj.basisValid) * proj.coeff(proj.basisValid, :));
            end
        end
        
        function proj = orthonormalize(proj)
            proj.warnIfNoArgOut(nargout);
            proj.coeff = orth(proj.coeff);
        end
        
        function tf = testIsOrthogonal(proj)
            assert(proj.initialized, 'Call after building / initializing');
            
            thresh = 1e-10;
            dp = proj.coeff' * proj.coeff;
            dp = abs(dp - diag(diag(dp)));
            tf = max(dp(:)) < thresh;
        end
    
        function proj = set.coeff(proj, v)
            %assert(isequal(size(proj.coeff), size(v)), 'Size of coeff must be [%d %d]', size(proj.coeff, 1), size(proj.coeff, 2));
            proj.coeff = v;
        end
        
        function proj = orthonormalizeOutputBases(proj, basisIdx)
            proj.coeff(:, basisIdx) = orth(proj.coeff(:, basisIdx));
        end
    end
    
end
