classdef StateSpaceProjectionStatistics
% The following properties described the amount of variance
% captured along each basis in the projection. This technique
% borrows heavily from DPCA, in particular some code is taken
% wholesale from https://github.com/wielandbrendel/dPCA and then
% modified to provide some additional statistics and deal with edge
% cases. The math is presented in http://arxiv.org/abs/1410.6031
%
% In computing these variance statistics,
% there are two different ways to make use of the timepoints over
% which data were provided. In 'allTimepoints', we make use of
% every timepoint that exists for all neurons on EACH condition. In
% this case, if one condition is much longer than another, we'll
% take all timepoints from the longer condition, and fewer
% timepoints from the shorter condition. This enables us to make
% use of more of the data, and biases any projections designed to
% explain variance towards the longer condition which contributes
% more data. The other possibility, 'sharedTimepoints', uses only
% timepoints which are captured in ALL conditions. In this case,
% we'll only use timepoints in common between the short and long
% condition (and all other conditions). This uses less of the data,
% and has the effect of each condition contributing equally to the
% variances explained. The advantage of the sharedTimepoints
% approach is that it enables us to marginalize the variance along
% each covariate separately. Thus if conditions are described by
% stimulus and decision, we can separate variance due to time
% alone, due to the stimulus alone, and due to the decision alone,
% and due to combinations of these. To do this requires that all
% conditions have data at each timepoint being considered.
% Consequently, when looking at the marginalized variances, be
% sure to compare this with the corresponding _sharedTimepoints
% total variance measurements.
%
% pca: we compute the cumulative variance explained by pca bases 
%   here as a reference for the best-case scenario. a set of bases
%   whose cumulative variance is close to pca is doing about as
%   well as it could in terms of variance explained. For shared
%   timepoints, note that this is pca done only on those shared
%   timepoints, so PCA bases found on all timepoints may not reach
%   this bound since it was asked to explain timepoints no longer
%   being considered.
%
% signal variance: we use pairs of individual trials and trial counts
%   to estimate the amount of noise in our average data. Following 
%   http://arxiv.org/abs/1410.6031 we derive an upper bound on the amount
%   of variance that could be due to noise. We then can derive signal
%   variance estimates explained in the first k < K bases by subtracting
%   the maximal variance that could be explained by noise in the worst case
%   orientation of that noise (since we don't assume simulataneous
%   recordings) from the amount of total of variance explained by those
%   first k bases.
%
% cumulativeVar: cumulative variance explained by taking the first
%   k bases together. It is not always possible to differentiate this
%   to get back to the variance per basis, since bases may not be
%   orthogonal. In the case of signalVariance, differentiation is
%   not possible because we do not actually know how the noise
%   variance is oriented (when not recorded simultaneously),
%   and so the signal variance computation is sensible only in the
%   cumulative "taking these first k bases together" sense.
%
% sizes: 
%   M = number of marginalizations. this is determined both by the
%     number of condition axes and the axisCombinations setting.
%   K = number of output bases

    properties(SetAccess=protected) % basic metadata
        basisNamesSource
        basisNamesProj
        
        nBasesSource
        nBasesProj
        
        basisValid

        conditionDescriptor
        marginalizationNames
        combinedParams
        marginalizationColorsStored
        
        % stored parameters that were passed in
        axisAttributeSets % what attributes were grouped along each axis
        axesIgnore % manually specified axes ignored for marginalization
        axesCombineSpecificMarginalizations
        axesCombineAllMarginalizations
        combineAxesWithTime % true or false
        
        axisIncludeMask % axes considered for marginalization, typically those with size > 1
        axisNames % names of original axes
    end
    
    properties(Constant)
        propertyListSettings = {'basisNamesSource', 'basisNamesProj', 'nBasesSource', 'nBasesProj', ...
                'conditionDescriptor', 'marginalizationNames', 'combinedParams', 'marginalizationColorsStored'};
    end

    properties(Dependent)
        marginalizationColors
        hasSignalVariance
        hasStatisticsRandomized
        nRandomSamples
        
        nConditions
        
        nBasesValid
        
        % nBasesValidSource * nTimepoints * nConditions
        nDataPoints_allTimepoints 
        nDataPoints_sharedTimepoints
    end
            
    properties(Dependent)
        nMarginalizations
    end
    
    methods
        function v = get.nMarginalizations(s)
            if isempty(s.marginalizationNames)
                v = NaN;
            else
                v = numel(s.marginalizationNames);
            end
        end
        
        function tf = get.hasSignalVariance(s)
            tf = ~isempty(s.totalSignalVar_allTimepoints);
        end
        
        function tf = get.hasStatisticsRandomized(s)
            tf = ~isempty(s.statisticsRandomized);
        end
        
        function n = get.nRandomSamples(s)
            n = numel(s.statisticsRandomized);
        end
        
        function n = get.nConditions(s)
            n = s.conditionDescriptor.nConditions;
        end
        
        function n = get.nBasesValid(s)
            n = nnz(s.basisValid);
        end

        function n = get.nDataPoints_allTimepoints(s)
            n = s.nBasesValid * s.nTimepointsAllConditions_allTimepoints;
        end
        
        function n = get.nDataPoints_sharedTimepoints(s)
            n = s.nBasesValid * s.nTimepointsAllConditions_sharedTimepoints;
        end
    end
    
    methods
        function s = set.marginalizationColors(s, v)
            if ~isa(v, 'function_handle')
                assert(size(v, 1) == s.nMarginalizations && size(v, 2) == 3, 'Colors must be nMarginalizations==%d x 3', s.nMarginalizations);
            end
            s.marginalizationColors = v;
        end
        
        function v = get.marginalizationColors(s)
            if isempty(s.marginalizationColorsStored)
                %v = distinguishable_colors(s.nMarginalizations, [0, 0, 0; 1, 1, 1]);
                v = cbrewer('qual', 'Dark2', s.nMarginalizations);
            else
                v = s.marginalizationColorsStored;
                if isa(v, 'function_handle')
                    v = v(s.nMarginalizations);
                end
            end
        end
        
    end

    methods(Access=protected)
        function s = StateSpaceProjectionStatistics()
        end
    end
    
    methods
        function sCopy = copyExcludingComputedStatistics(s)
            props = StateSpaceProjectionStatistics.propertyListSettings;
            sCopy = StateSpaceProjectionStatistics();
            for iP = 1:numel(props)
                sCopy.(props{iP}) = s.(props{iP});
            end
        end
    end
    
    properties(SetAccess=protected) % computed statistics
        % nBases x nBases covariance matrix
        corrSource
        covSource
        covMarginalized
        
        % variance estimates
        
        %%%%%
        % All timepoints
        %%%%%
        
        % number of total timepoints across all conditions in the "all
        % timepoints" case. Note that different conditions contribute differing
        % numbers of timepoints.
        nTimepointsAllConditions_allTimepoints % scalar
        timeMaskByCondition_allTimepoints % T x ConditionAttr logical
        
        % total variances
        totalVar_allTimepoints % scalar
        totalSignalVar_allTimepoints % scalar (requires single trials)
        totalNoiseVar_allTimepoints % scalar (requires single trials)
        
        % variance per component
        componentVarByBasis_allTimepoints % 1 x K (do not cumsum these to get cumulative variance)
        
        % cumulative raw variance
        cumulativeVarByBasis_allTimepoints % 1x K
        cumulativeNoiseVarByBasis_allTimepoints % 1 x K
        cumulativeSignalVarByBasis_allTimepoints % 1 x K
        
        % cumulative fraction 0 --> 1
        cumulativeFractionVarByBasis_allTimepoints % 1 x K
        cumulativeFractionSignalVarByBasis_allTimepoints % 1 x K
        
        % pca upper bounds
        pca_cumulativeVarByBasis_allTimepoints % 1 x K
        pca_cumulativeSignalVarByBasis_allTimepoints % 1 x K
        pca_cumulativeFractionVarByBasis_allTimepoints % 1 x K
        pca_cumulativeFractionSignalVarByBasis_allTimepoints % 1 x K
        
        %%%%%
        % Shared timepoints
        %%%%%
        
        nTimepointsAllConditions_sharedTimepoints % scalar
        nTimepointsPerCondition_sharedTimepoints % scalar
        timeMask_sharedTimepoints % T x 1 logical
        
        % total variances
        totalVar_sharedTimepoints % scalar
        totalSignalVar_sharedTimepoints % scalar
        totalNoiseVar_sharedTimepoints % scalar
        
        % variance per component
        componentVarByBasis_sharedTimepoints % 1 x K (do not cumsum these to get cumulative variance)
        
        % cumulative raw variance
        cumulativeVarByBasis_sharedTimepoints
        cumulativeSignalVarByBasis_sharedTimepoints
        cumulativeNoiseVarByBasis_sharedTimepoints

        % cumulative fraction 0 -> 1
        cumulativeFractionVarByBasis_sharedTimepoints
        cumulativeFractionSignalVarByBasis_sharedTimepoints
        
        % pca upper bounds
        pca_cumulativeVarByBasis_sharedTimepoints
        pca_cumulativeSignalVarByBasis_sharedTimepoints
        pca_cumulativeFractionVarByBasis_sharedTimepoints
        pca_cumulativeFractionSignalVarByBasis_sharedTimepoints
        
        %%%%%
        % Marginalized variance in shared timepoints
        %%%%%
        
        % total variances
        totalMarginalizedVar_sharedTimepoints
        totalMarginalizedNoiseVar_sharedTimepoints
        totalMarginalizedSignalVar_sharedTimepoints
        
        % variance per component
        componentMarginalizedVarByBasis_sharedTimepoints
            
        % cumulative raw variance
        cumulativeMarginalizedVarByBasis_sharedTimepoints
        cumulativeMarginalizedSignalVarByBasis_sharedTimepoints
        cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints
        
        % cumulative fraction 0 --> 1
        cumulativeFractionMarginalizedSignalVarByBasis_sharedTimepoints
        cumulativeFractionMarginalizedVarByBasis_sharedTimepoints
        
        % StateSpaceProjectionStatistics instances for each randomization
        % in dataMeanRandomized
        statisticsRandomized
    end

    methods % Post-build utilities
        function idxList = lookupMarginalizationIdx(marginalizationSpec)
            if StateSpaceProjectionStatistics.isMarginalizationSpec(marginalizationSpec)
                % convert to list so we can loop over it
                marginalizationSpec = {marginalizationSpec};
            elseif ~StateSpaceProjectionStatistics.isListOfMarginalizationSpec(marginalizationSpec)
                error('Input must be marginalization spec or list of marginalizationSpec');
            end
            
            for iM = 1:numel(marginalizationSpec)
                margSpec = marginalizationSpec{iM};
                
            end
        end
    end
    
    methods % Plotting results 
        function plotCumulativeVariance(s, varargin)
            p = inputParser();
            p.addParameter('fractional', true, @islogical);
            p.addParameter('showPCABound', true, @islogical);
            p.addParameter('timepoints', 'all', @ischar); % or 'shared'
            p.addParameter('varianceType', 'total', @ischar);
            p.addParameter('marginalize', [], @(x) isempty(x) || islogical(x));
            p.addParameter('showRandomQuantiles', 0.95, @isscalar); % show traces individual for individual randomization statistics
            p.addParameter('showTotals', true, @islogical); % show horizontal lines at total for each variance
            p.parse(varargin{:});
            
            varianceType = p.Results.varianceType;
            marginalize = p.Results.marginalize;
            timepoints = p.Results.timepoints;
            
            markerSize = 4;
            
            if isempty(marginalize)
                marginalize = strcmp(timepoints, 'shared');
            elseif marginalize
                assert(strcmp(timepoints, 'shared'), 'Timepoints must be set to shared in order to show marginalized variances');
            end
            
            axh = newplot;
            
            switch varianceType
                case 'total'
                    varStr = 'Var';
                    varStrP = 'Variance';
                 case 'signal'
                    assert(s.hasSignalVariance, 'Signal variance not computed, PopulationTrajectorySet must have single trials for this to be possible');
                    varStr = 'SignalVar';
                    varStrP = 'Signal Variance';
                case 'noise' 
                    assert(s.hasSignalVariance, 'Noise variance not computed, PopulationTrajectorySet must have single trials for this to be possible');
                    varStr = 'NoiseVar';
                    varStrP = 'Noise Variance';
                otherwise
                    error('Unknown varianceType argument %s, supported values are total, signal, noise', p.Results.varianceType);
            end
            
            switch timepoints
                case 'all'
                    timeStr = 'allTimepoints';
                    timeStrP = 'All Timepoints';
                case 'shared';
                    timeStr = 'sharedTimepoints';
                    timeStrP = 'Shared Timepoints';
                otherwise
                    error('Unknown timepoints argument %s, supported values are all and shared', p.Results.timepoints);
            end
            
            if p.Results.fractional
                fracStr = 'Fraction';
                fracStrP = 'Fraction ';
            else
                fracStr = '';
                fracStrP = '';
            end
            
            fieldName = sprintf('cumulative%s%sByBasis_%s', fracStr, varStr, timeStr);
            v = s.(fieldName);
            
            h = plot(axh, 1:s.nBasesProj, v, 'o-', ...
                'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.1 0.1 0.1], ...
                'MarkerEdgeColor', 'none', 'MarkerSize', markerSize, 'LineWidth', 1);
            TrialDataUtilities.Plotting.showInLegend(h, 'Total');

             hold(axh, 'on');
            
            if marginalize
                mFieldName = sprintf('cumulative%sMarginalized%sByBasis_%s', fracStr, varStr, timeStr);
                vm = s.(mFieldName);
                
                for m = 1:s.nMarginalizations
                    color = s.marginalizationColors(m, :);
                    h = plot(axh, 1:s.nBasesProj, vm(m, :), 'o-', ...
                    'Color', color, 'MarkerFaceColor', color, ...
                    'MarkerEdgeColor', 'none', 'MarkerSize', markerSize, 'LineWidth', 1);
                    TrialDataUtilities.Plotting.showInLegend(h, s.marginalizationNames{m});
                end
            else
                mFieldName = ''; %#ok<NASGU>
            end
            
            if p.Results.showPCABound
                pcaFieldName = ['pca_' fieldName];
                try
                    vp = s.(pcaFieldName);
%                     vp = s.(pcaFieldName);
                    h = plot(axh, 1:s.nBasesProj, vp, '-', ...
                        'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                    TrialDataUtilities.Plotting.showInLegend(h, 'PCA Upper Bound');
                catch
                    pcaFieldName = ''; %#ok<NASGU>
                end
            else
                pcaFieldName = ''; %#ok<NASGU>
            end
            
            if ~isempty(p.Results.showRandomQuantiles) && ~isempty(s.statisticsRandomized) && ...
                    strcmp(varianceType, 'total')
                % if asked, show each randomization quantiles, should be
                % nBasesProj x nQuantiles matrix
                randKbyQ = TensorUtils.squeezeDims(s.computeQuantilesForRandomizedStatistics(fieldName, p.Results.showRandomQuantiles), 1);
                plot(axh, 1:s.nBasesProj, randKbyQ, '-', 'Color', [0.5 0.5 0.5]);
            end
            
            if p.Results.showTotals
                fieldName = sprintf('total%s_%s', varStr, timeStr);
                total = s.(fieldName);
                if ~p.Results.fractional
                    h = TrialDataUtilities.Plotting.horizontalLine(total, ...
                        'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'LineStyle', '--');
                    TrialDataUtilities.Plotting.hideInLegend(h);
                end
                
                if marginalize
                    mFieldName = sprintf('totalMarginalized%s_%s', varStr, timeStr);
                    totalMarg = s.(mFieldName);
                    if p.Results.fractional
                        totalMarg = totalMarg / total;
                    end
                    
                    for m = 1:s.nMarginalizations
                        color = s.marginalizationColors(m, :);
                        h = TrialDataUtilities.Plotting.horizontalLine(totalMarg(m), ... 
                            'LineStyle', '--', 'Color', color, 'LineWidth', 0.5);
                        TrialDataUtilities.Plotting.hideInLegend(h);
                    end
                end
            end
            
            xlabel(axh, 'Number of Bases');
            ylabel(axh, sprintf('Cumulative %sVariance Explained', fracStrP))
            titleStr = sprintf('Cumulative %s%s, %s', fracStrP, varStrP, timeStrP);
            title(titleStr);
            
            hold(axh, 'off');
            xlim([0 s.nBasesProj]);
            if p.Results.fractional
                ylim([0 1]);
            end
            
            au = AutoAxis.replace(axh);
            
            % legend
            if marginalize
                au.addColoredLabels(cat(1, {'total'}, s.marginalizationNames), cat(1, [0 0 0], s.marginalizationColors), ...
                    'posY', AutoAxis.PositionType.Top, 'posX', AutoAxis.PositionType.Left);
                au.update();
            end
            
            % fieldsPlotted = {fieldName; mFieldName; pcaFieldName}
        end
        
        function plotCumulativeSignalVariance_SharedTimepoints(s, varargin)
            s.plotCumulativeVariance('timepoints', 'shared', 'marginalize', true, 'fractional', true, 'varianceType', 'signal', varargin{:});
        end
        
        function plotCumulativeVariance_SharedTimepoints(s, varargin)
            s.plotCumulativeVariance('timepoints', 'shared', 'marginalize', true, 'fractional', true, 'varianceType', 'total', varargin{:});
        end
        
        function plotCumulativeSignalVariance_AllTimepoints(s, varargin)
            s.plotCumulativeVariance('timepoints', 'all', 'fractional', true, 'varianceType', 'signal', varargin{:});
        end
        
        function plotCumulativeVariance_AllTimepoints(s, varargin)
            s.plotCumulativeVariance('timepoints', 'all', 'fractional', true, 'varianceType', 'total', varargin{:});
        end
        
        function plotCovSource(proj, varargin)
            clf;
            pmat(proj.covSource);
            box off;
            title('Source Covariance');
            hold off;
        end
        
        function plotBasisMixtures(s, varargin)
            p = inputParser;
            p.addParameter('basisIdx', 1:min([20 s.nBasesProj]), @(x) isvector(x) && ...
                all(inRange(x, [1 s.nBasesSource])));
            p.parse(varargin{:});
            basisIdx = p.Results.basisIdx;
            if islogical(basisIdx)
                basisIdx = find(basisIdx);
            end
            
            axh = newplot;
            hold(axh, 'on');
            
            cumBasisMix = cumsum(s.componentMarginalizedVarByBasis_sharedTimepoints, 1)';
            nCov = s.nMarginalizations;
            nBases = length(basisIdx);
            rowHeight = 0.8;
            xMin = 0;
            xMax = max(cumBasisMix(:));
            hPatch = nan(nCov, 1);
            
            cmap = s.marginalizationColors;
            
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
                    
                    hPatch(iCov) = patch(patchX, patchY, cmap(iCov, :), 'Parent', axh, 'EdgeColor', 'none');
                end
                
%                 h = text(0, nBases-iB+rowHeight/2, s.basisNamesProj{iB}, ...
%                     'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Right');
%                 extent = get(h, 'Extent');
%                 %xMin = min(extent(1), xMin);
                
                componentFractionTotalVariance = s.componentVarByBasis_sharedTimepoints(iB) / s.totalVar_sharedTimepoints;
                text(cumBasisMix(iB, end), nBases-iB+rowHeight/2, sprintf('  %.2f%%', componentFractionTotalVariance*100), ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Left', 'Parent', axh); 
               % extent = get(h, 'Extent');
                %xMax = max(extent(1)+extent(3), xMax);
            end
            
            xlim([xMin xMax*1.1]);
            ylim([0 nBases]);
            axis off;
            title('Basis Mixtures');
            
            au = AutoAxis(gca);
            au.addTitle();
            au.addTicklessLabels('y', 'tick', nBases-0.5:-1:0.5, 'tickLabel', s.basisNamesProj);
            
            au.addColoredLabels(s.marginalizationNames, cmap, ...
                'posY', AutoAxis.PositionType.Bottom, 'posX', AutoAxis.PositionType.Right);

            au.update();
            au.installCallbacks();
            hold off;
        end
    end
    
    methods(Static) % Build from StateSpaceProjection and PopulationTrajectorySet
        function s = build(proj, pset, varargin)
            p = inputParser();
            p.addParameter('meanSubtract', true, @islogical); % setting this to false only makes sense for situations where the data is already normalized relative to some absolute baseline, such as a difference between two conditions
            p.addParameter('computeForRandomized', true, @islogical);
            p.addParameter('axesIgnore', {}, @(x) true);
            p.addParameter('axesCombineSpecificMarginalizations', {}, @(x) true);
            p.addParameter('axesCombineAllMarginalizations', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('combineAxesWithTime', true, @(x) islogical(x) || iscell(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            s = StateSpaceProjectionStatistics();
            
            % get the names and counts of source -> proj bases
            s.basisNamesSource = pset.basisNames;
            s.basisNamesProj = proj.getBasisNames(pset);
            s.nBasesSource = proj.nBasesSource;
            
            s.basisValid = proj.basisValid;
            s.nBasesProj = proj.nBasesProj;
            s.conditionDescriptor = pset.conditionDescriptor;
                      
            % filter for non-singular axes
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            
            % build the list of covariates and covariate interactions to marginalize along
            s.axesIgnore = p.Results.axesIgnore;
            s.axesCombineSpecificMarginalizations = p.Results.axesCombineSpecificMarginalizations;
            s.axesCombineAllMarginalizations = p.Results.axesCombineAllMarginalizations;
            s.combineAxesWithTime = p.Results.combineAxesWithTime;
            s.axisNames = pset.conditionDescriptor.axisNames;
            s.axisIncludeMask = dimMask;
            [s.combinedParams, s.marginalizationNames] = StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization( ...
                pset.conditionDescriptor.axisAttributes, ...
                'axisIncludeMask', dimMask, ...
                'axisNames', pset.conditionDescriptor.axisNames, ...
                'axesIgnore', p.Results.axesIgnore, ...
                'axesCombineSpecificMarginalizations', p.Results.axesCombineSpecificMarginalizations, ...
                'axesCombineAllMarginalizations', p.Results.axesCombineAllMarginalizations, ...
                'combineAxesWithTime', p.Results.combineAxesWithTime);
                   
            % Extract trial averaged data
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            
            % Save covariance matrices
            CTAbyNv = pset.buildCTAbyN('validBasesOnly', true);
            if size(CTAbyNv, 1) == 1
                CTAbyNv = repmat(CTAbyNv, 2, 1); % some methods below transpose automatically if the matrix looks like a row vector
            end
            
            s.covSource = TensorUtils.inflateMaskedTensor(nancov(CTAbyNv, 1), [1 2], proj.basisValid); 
            s.corrSource = TensorUtils.inflateMaskedTensor(corrcoef(CTAbyNv, 'rows', 'complete'), [1 2], proj.basisValid);
            covMarginalizedValidCell = TrialDataUtilities.DPCA.dpca_marginalizedCov(NvbyTAbyAttr, 'combinedParams', s.combinedParams);
            s.covMarginalized = cellfun(@(x) TensorUtils.inflateMaskedTensor(x, [1 2], proj.basisValid), ...
                covMarginalizedValidCell, 'UniformOutput', false);

            % filter encoder / decoder by basis valid
            encoderNvbyK = proj.encoderNbyK(proj.basisValid, :);
            decoderKbyNv = proj.decoderKbyN(:, proj.basisValid);
            
            % if possible, use individual trials to get signal variance as
            % well
            if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimate)
                % grab difference of trials scaled estimate, and shape into
                % tensor (N x TA x C --> N x TA x Attr)
                scaledNoiseEstimate_NbyTAbyAttr = reshape(pset.dataDifferenceOfTrialsScaledNoiseEstimate, ...
                    [pset.nBases, sum(pset.nTimeDataMean), pset.conditionDescriptor.conditionsSize]);
                scaledNoiseEstimate_NvbyTAbyAttr = TensorUtils.selectAlongDimension(scaledNoiseEstimate_NbyTAbyAttr, 1, proj.basisValid);
                
                s = s.computeStatistics(NvbyTAbyAttr, ...
                    decoderKbyNv, encoderNvbyK, 'combinedParams', s.combinedParams, ...
                    'scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr', scaledNoiseEstimate_NvbyTAbyAttr, ...
                    'meanSubtract', p.Results.meanSubtract);
            else
                s = s.computeStatistics(NvbyTAbyAttr, ...
                    decoderKbyNv, encoderNvbyK, 'combinedParams', s.combinedParams, ...
                    'meanSubtract', p.Results.meanSubtract);
            end
            
            if pset.hasDataRandomized && p.Results.computeForRandomized
                % do same statistics computation on 
                prog = ProgressBar(pset.nRandomSamples, 'Computing projection statitics on dataMeanRandomized');
                for iR = 1:pset.nRandomSamples
                    prog.update(iR);
                    NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true, ...
                        'type', 'meanRandom', 'dataRandomIndex', iR);
                    
                    sRand(iR) = s.copyExcludingComputedStatistics(); %#ok<AGROW>
                    sRand(iR) = sRand(iR).computeStatistics(NvbyTAbyAttr, ...
                        decoderKbyNv, encoderNvbyK, 'combinedParams', s.combinedParams, ...
                        'meanSubtract', p.Results.meanSubtract, ...
                        'verbose', false, ...
                        'scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr', pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized(:, :, :, iR)); %#ok<AGROW>
                end
                prog.finish();
                
                s.statisticsRandomized = makecol(sRand);
            end
        end
    end
    
    methods(Access=protected) % internal statistics computation function
        function s = computeStatistics(s, NbyTAbyAttr, decoderKbyN, encoderNbyK, varargin)
            % based heavily on math inside dpca_explainedVariance, though modified heavily
            %
            % stats = computeProjectionStatistics(NbyTAbyAttr, decoderKbyN, encoderNbyK, varargin) 
            %   computes various measures and stores them in
            %   dataNbyT is the data
            %   matrix, decoder is the decoder matrix (K x N), encoder is the encoder matrix (N x K).
            %
            %  'combinedParams' - cell array of cell arrays specifying 
            %                     which marginalizations should be added up together,
            %                     e.g. for the three-parameter case with parameters
            %                           1: stimulus
            %                           2: decision
            %                           3: time
            %                     one could use the following value:
            %                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.
            %
            %  'singleTrials_NbyTAbyRbyAttr'        - array of single trials. Has one extra dimension as 
            %                     compared with X and stores individual single trial
            %                     firing rates, as opposed to the trial average. If
            %                     provided, "signal variance" will be computed. Must
            %                     have **exactly two trials** per neuron, per
            %                     condition. These two trials should be
            %                     chosen randomly.
            %
            %  'trialCountsTotal_NbyAttr' - must be provided. total number of trials
            %                     taken into average
            %

            % restore original input shape
            N = size(decoderKbyN,2);
            K = size(decoderKbyN,1);
            maxNK = max(K, N);

            assert(size(NbyTAbyAttr, 1) == N, 'Shape of decoder must be K by N');
            assert(size(encoderNbyK, 1) == N, 'Shape of encoder must be N by K');
            assert(size(decoderKbyN, 1) == size(encoderNbyK, 2), 'Encoder and decoder differ on size K');
            
            szVec = TensorUtils.sizeNDims(NbyTAbyAttr, 3);
            nConditions = prod(szVec(3:end));
            
            p = inputParser();
            p.addParameter('meanSubtract', true, @islogical); % setting this to false only makes sense for situations where the data is already normalized relative to some absolute baseline, such as a difference between two conditions
            p.addParameter('combinedParams', {}, @(x) true);
            p.addParameter('scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr', [], @(x) isempty(x) || isnumeric(x));
            p.addParameter('verbose', true, @islogical);
            p.addParameter('marginalize', true, @islogical);
            p.addParameter('pcaBound', true, @islogical);
            p.parse(varargin{:});

            dataTensor = NbyTAbyAttr;
            
            %%%%%%%%
            % All timepoints variances
            %%%%%%%%
            
            dataNbyT = dataTensor(:,:);

            % dataNbyT is now N x T
            % remove nan timepoints from the flattened traces, this enables all
            % timepoints from all conditions to be used as long as all neurons have a
            % sample at that time
            allTimepoints_NxT_keepMaskT = makecol(~any(isnan(dataNbyT), 1));
            dataNbyT_allTimepoints = dataNbyT(:, allTimepoints_NxT_keepMaskT);
            
            s.nTimepointsAllConditions_allTimepoints = nnz(allTimepoints_NxT_keepMaskT);
            
            % go from TC x 1 to T x Cattr
            s.timeMaskByCondition_allTimepoints = reshape(allTimepoints_NxT_keepMaskT, szVec(2:end));
            
            % subtract means for "all timepoints" over time for each basis,
            % condition
            if p.Results.meanSubtract
                dataNbyT_allTimepoints = bsxfun(@minus, dataNbyT_allTimepoints, mean(dataNbyT_allTimepoints, 2));
            end

            % total variance
            s.totalVar_allTimepoints = sum(dataNbyT_allTimepoints(:).^2);

            if p.Results.pcaBound
                % PCA explained variance
                if p.Results.verbose, debug('Computing trial-average SVD\n'); end
                Spca_all = svd(dataNbyT_allTimepoints', 0);
                if numel(Spca_all) < K
                    Spca_all = cat(1, Spca_all, zerosvec(K-numel(Spca_all)));
                end
                s.pca_cumulativeVarByBasis_allTimepoints = cumsum(Spca_all.^2');
                s.pca_cumulativeFractionVarByBasis_allTimepoints = s.pca_cumulativeVarByBasis_allTimepoints / s.totalVar_allTimepoints;
            end
            
            % explained variance along cumulative sets of bases
            Z = decoderKbyN*dataNbyT_allTimepoints;
            s.componentVarByBasis_allTimepoints = sum(Z.^2, 2)';
            
            prog = ProgressBar(size(decoderKbyN,1), 'Computing cumulative variance explained all timepoints');
            for i=1:size(decoderKbyN,1)
                prog.update(i);
                s.cumulativeVarByBasis_allTimepoints(i) = s.totalVar_allTimepoints - sum(sum((dataNbyT_allTimepoints - encoderNbyK(:,1:i)*Z(1:i,:)).^2));    
            end
            prog.finish();
            s.cumulativeFractionVarByBasis_allTimepoints = s.cumulativeVarByBasis_allTimepoints / s.totalVar_allTimepoints;

            %%%%%%%
            % Shared timepoints variance
            %%%%%%%
            
            % find timepoints where all conditions that have any samples and all neurons have samples, any
            % along non-time dimensions
            dataTensor_NxTxC = dataTensor(:, :, :);
            maskConditionsConsider = all(any(~isnan(dataTensor_NxTxC), 2), 1);
            sharedTimepoints_tensor_keepMaskT = makecol(squeeze(~TensorUtils.anyMultiDim(isnan(dataTensor_NxTxC(:, :, maskConditionsConsider)), [1 3])));
            dataTensor_shared = TensorUtils.selectAlongDimension(dataTensor, 2, sharedTimepoints_tensor_keepMaskT);

            s.nTimepointsPerCondition_sharedTimepoints = nnz(sharedTimepoints_tensor_keepMaskT);
            s.nTimepointsAllConditions_sharedTimepoints = s.nTimepointsPerCondition_sharedTimepoints * nConditions;
            s.timeMask_sharedTimepoints = sharedTimepoints_tensor_keepMaskT;
            
            if any(sharedTimepoints_tensor_keepMaskT)
                % subtract means for "all timepoints" over time for each basis,
                % condition. do this once for the tensor
                if p.Results.meanSubtract
                    dataTensor_shared = TensorUtils.centerSlicesOrthogonalToDimension(dataTensor_shared, 1);
                end

                % dataNxT_shared has only participating conditions
                dataNxTxC_shared = dataTensor_shared(:, :, maskConditionsConsider);
                dataNxT_shared = dataNxTxC_shared(:, :);
                
                s.totalVar_sharedTimepoints = sum(dataNxT_shared(:).^2);
                    
                if p.Results.pcaBound
                    % PCA explained variance upper bound
                    if p.Results.verbose, debug('Computing trial-average SVD shared timepoints\n'); end

                    % the size of this will be min(N, T) x 1, need it to be
                    % at least length K
                    Sshared = svd(dataNxT_shared', 0);
                    if numel(Sshared) <= K
                        Sshared = cat(1, Sshared, zerosvec(K-numel(Sshared)));
                    end
                    s.pca_cumulativeVarByBasis_sharedTimepoints = cumsum(Sshared(1:K).^2');
                    s.pca_cumulativeFractionVarByBasis_sharedTimepoints = s.pca_cumulativeVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;
                end
                
                % explained variance along bases
                prog = ProgressBar(size(decoderKbyN, 1), 'Computing cumulative variance explained shared timepoints');
                Z = decoderKbyN*dataNxT_shared;
                s.componentVarByBasis_sharedTimepoints = sum(Z.^2, 2)'; % along each basis INDIVIDUALLY
                for i=1:K
                    prog.update(i);
                    s.cumulativeVarByBasis_sharedTimepoints(i) = s.totalVar_sharedTimepoints - sum(sum((dataNxT_shared - encoderNbyK(:,1:i)*Z(1:i,:)).^2));    
                end
                prog.finish();
                s.cumulativeFractionVarByBasis_sharedTimepoints = s.cumulativeVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;
            
                %%%%%%
                % Marginalized Variances
                %%%%%%

                if p.Results.marginalize
                    % for marginalization though, we need to only use timepoints where all conditions
                    % are present (non-nan). This means the total variance over all
                    % marginalizations will be different than for the non marginalized
                    % variance, since we throw away timepoints not shared across all conditions
                    % hence the _sharedTimepoints vs allTimepoints distinction

                    % marginalizing the trimmed tensor, should already be
                    % mean-subtracted
                    dataMarginalizedTensors_shared = TrialDataUtilities.DPCA.dpca_marginalize(dataTensor_shared, ...
                        'meanSubtract', false, ... % should already be mean subtracted since it's a difference
                        'combinedParams', p.Results.combinedParams, 'ifFlat', 'no');
                    nMarginalizations = length(dataMarginalizedTensors_shared); %#ok<*PROPLC>

                    Zshared = decoderKbyN*dataNxT_shared;
                    s.componentVarByBasis_allTimepoints = sum(Zshared.^2, 2)';

                    % total marginalized variance
                    s.totalMarginalizedVar_sharedTimepoints = nan(nMarginalizations, 1);
                    dataMarg_NxTxC = cellvec(nMarginalizations);
                    for i=1:nMarginalizations
                        dataMarg_NxTxC{i} = dataMarginalizedTensors_shared{i}(:, :, maskConditionsConsider);
                        s.totalMarginalizedVar_sharedTimepoints(i) = sum(dataMarg_NxTxC{i}(:).^2);
                    end

                    % marginalized variance of each component : shared timepoints
                    for i=1:nMarginalizations
                        s.componentMarginalizedVarByBasis_sharedTimepoints(i,:) = sum((decoderKbyN * dataMarg_NxTxC{i}(:, :)).^2, 2)'; % @ djoshea shouldn't be normalized
                    end

                    % explained variance along cumulative sets of bases
                    prog = ProgressBar(nMarginalizations, 'Computing variance over %d marginalizations', nMarginalizations);
                    for m=1:nMarginalizations
                        prog.update(m);
%                         thisMarg = dataMarginalizedTensors_shared{m}(:, :);
                        thisMarg = dataMarg_NxTxC{m}(:, :);
                        totalVarThisMarg = sum(thisMarg(:).^2);
                        Zmarg = decoderKbyN*thisMarg;
                        progK = ProgressBar(K, 'Computing cumulative marginalized variance');
                        for i=1:K
                            progK.update(i);
                            s.cumulativeMarginalizedVarByBasis_sharedTimepoints(m, i) = totalVarThisMarg - sum(sum((thisMarg - encoderNbyK(:,1:i)*Zmarg(1:i,:)).^2));    
                        end
                        progK.finish();
                    end
                    prog.finish();
                    s.cumulativeFractionMarginalizedVarByBasis_sharedTimepoints = s.cumulativeMarginalizedVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;
                end
            else
                warning('No shared timepoints found in PopulationTrajectorySet, skipping marginalization');
            end
            
            %%%%%%%
            % Noise and signal variances, if single trial data provided
            %%%%%%%
            if ~isempty(p.Results.scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr)
                if p.Results.verbose, debug('Computing signal variance via noise-floor\n'); end
                noiseTensor = p.Results.scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr;
                
                % filter using the same mask we used before for the all timepoint
                noiseNxT = noiseTensor(:, :);
                noiseNxT_all = noiseNxT(:, allTimepoints_NxT_keepMaskT);

                % this may need to be replaced with something better, there's no
                % guarantee that individual trials won't have NaNs at places where the
                % trial-averaged traces do not, since there could still be enough
                % trials at that timepoint to form a non-NaN average. However, if we
                % drop additional timepoints here, then we'll have less total variance
                % to explain and need to correct for this. #todo
                nanMaskT = any(isnan(noiseNxT_all), 1);
                if all(nanMaskT(:))
                    warning('Single trial data has NaN values at all timepoints where trial-averaged data did not (%d / %d time points). Noise variance will not be computed.\n', nnz(nanMaskT), numel(nanMaskT));
                    noiseNxT_all_nonNan = noiseNxT_all(:, ~nanMaskT);
                elseif any(nanMaskT(:))
                    debug('Single trial data has NaN values at timepoints where trial-averaged data did not (%d / %d time points). Noise variance will be normalized accordingly.\n', nnz(nanMaskT), numel(nanMaskT));
                    noiseNxT_all_nonNan = noiseNxT_all(:, ~nanMaskT);
                else
                    noiseNxT_all_nonNan = noiseNxT_all;
                end
                noiseVarMultiplier_allTimepoints = numel(nanMaskT) / nnz(~nanMaskT);

                % we mean subtract here because we will also mean subtract
                % when we look at noise in shared timepoints, and we want
                % the noise variances to be comparable when all timepoints
                % are shared timepoints. we mean subtract there because the
                % marginalization essentially accomplishes the same effect,
                %Removed since this is
    %                 already the difference between two trials so
    %                 distribution mean should be zero
%                 noiseNxT_all_nonNan = bsxfun(@minus, noiseNxT_all_nonNan, mean(noiseNxT_all_nonNan, 2));

                %%%%%%
                % Noise variance, all timepoints
                %%%%%%
                
                if ~all(nanMaskT(:))
                     % total noise variance
                    s.totalNoiseVar_allTimepoints = sum(noiseNxT_all_nonNan(:).^2) * noiseVarMultiplier_allTimepoints;

                    % total signal variance
                    s.totalSignalVar_allTimepoints = max(0, s.totalVar_allTimepoints - s.totalNoiseVar_allTimepoints);

                    % PCA explained signal variance
                    if p.Results.verbose, debug('Computing noise SVD\n'); end
                    Snoise_all = svd(noiseNxT_all_nonNan', 0) * noiseVarMultiplier_allTimepoints;
                    if numel(Snoise_all) < K % account for too many output bases
                        Snoise_all = cat(1, Snoise_all, zerosvec(K-numel(Snoise_all)));
                    end
                else
                    s.totalNoiseVar_allTimepoints = NaN;
                    s.totalSignalVar_allTimepoints = NaN;
                    Snoise_all = nanvec(maxNK);
                end
                
                pcaSignal = Spca_all(1:K).^2 - Snoise_all(1:K).^2;
                s.pca_cumulativeSignalVarByBasis_allTimepoints = min(s.totalSignalVar_allTimepoints, cumsum(pcaSignal'));
                s.pca_cumulativeSignalVarByBasis_allTimepoints = TensorUtils.makeNonDecreasing(s.pca_cumulativeSignalVarByBasis_allTimepoints);
                
                s.pca_cumulativeFractionSignalVarByBasis_allTimepoints = s.pca_cumulativeSignalVarByBasis_allTimepoints / s.totalSignalVar_allTimepoints;
                
                % variance explained cumulatively over basis variance
                for i=1:K
                    s.cumulativeNoiseVarByBasis_allTimepoints(i) = sum(Snoise_all(1:i).^2);
                end
                s.cumulativeSignalVarByBasis_allTimepoints = min(s.totalSignalVar_allTimepoints, ...
                    s.cumulativeVarByBasis_allTimepoints - s.cumulativeNoiseVarByBasis_allTimepoints);
                
                % ensure cumulative signal var is non-decreasing and adjust
                % noise var to ensure signal+noise = total
                s.cumulativeSignalVarByBasis_allTimepoints = TensorUtils.makeNonDecreasing(max(0, s.cumulativeSignalVarByBasis_allTimepoints));
                s.cumulativeNoiseVarByBasis_allTimepoints = s.cumulativeVarByBasis_allTimepoints - s.cumulativeSignalVarByBasis_allTimepoints;
                
                % compute fraction of total signal var
                s.cumulativeFractionSignalVarByBasis_allTimepoints = s.cumulativeSignalVarByBasis_allTimepoints / s.totalSignalVar_allTimepoints;
                if s.totalSignalVar_allTimepoints == 0
                    % otherwise will be nan
                    s.cumulativeFractionSignalVarByBasis_allTimepoints = 0 * s.cumulativeSignalVarByBasis_allTimepoints;
                    s.pca_cumulativeFractionSignalVarByBasis_allTimepoints = 0 * s.cumulativeSignalVarByBasis_allTimepoints;
                end
                %%%%%%%%%%
                % Noise variance, shared timepoints
                %%%%%%%%%%

                if any(sharedTimepoints_tensor_keepMaskT)
                    % noiseTensor is NbyTAbyAttr
                    % filter in time to match dataTensor_shared and hopefully remove all NaNs 
                    % in shared timepoints
                    noiseTensor_shared = TensorUtils.selectAlongDimension(noiseTensor, 2, ...
                        sharedTimepoints_tensor_keepMaskT);
                    
    %                 % subtract mean by basis now - Removed since this is
    %                 already the difference between two trials so
    %                 distribution mean should be zero
    %                 noiseTensor_shared = bsxfun(@minus, noiseTensor_shared, nanmean(noiseTensor_shared,2));

                    % this may need to be replaced with something better, there's no
                    % guarantee that individual trials won't have NaNs at places where the
                    % trial-averaged traces do not, since there could still be enough
                    % trials at that timepoint to form a non-NaN average. However, if we
                    % drop additional timepoints here, then we'll have less total variance
                    % to explain and need to correct for this. Thus the
                    % correction factor noiseVarMultiplier below
                    noiseTensor_shared_NxTxC = noiseTensor_shared(:, :, maskConditionsConsider);
                    nanMaskT = makecol(squeeze(TensorUtils.anyMultiDim(isnan(noiseTensor_shared_NxTxC), [1 3])));
                    %  nanMaskT = TensorUtils.anyMultiDim(isnan(noiseTensor_shared), [1 3:ndims(noiseTensor_shared)]); % old version
                    if all(nanMaskT)
                        warning('Single trial data has NaN values at all shared timepoints where trial-averaged data did not (%d / %d time points). Noise variance will not be computed.\n', nnz(nanMaskT), numel(nanMaskT));
                        noiseTensor_shared_maskT = TensorUtils.selectAlongDimension(noiseTensor_shared, 2, ~nanMaskT);
                    elseif any(nanMaskT)
                        debug('Single trial data has NaN values at shared timepoints where trial-averaged data did not (%d / %d time points). Noise variance will be normalized accordingly.\n', nnz(nanMaskT), numel(nanMaskT));
                        noiseTensor_shared_maskT = TensorUtils.selectAlongDimension(noiseTensor_shared, 2, ~nanMaskT);
                    else
                        noiseTensor_shared_maskT = noiseTensor_shared;
                    end
                    noiseVarMultiplier_sharedTimepoints = numel(nanMaskT) / nnz(~nanMaskT);

                    noiseData_NxTxC_nonNaN = noiseTensor_shared_maskT(:, :, maskConditionsConsider);
                    if ~all(nanMaskT)
                        s.totalNoiseVar_sharedTimepoints = sum(noiseData_NxTxC_nonNaN(:).^2) * noiseVarMultiplier_sharedTimepoints;

                        % ensure signal var is nonnegative using max(0, ...) and
                        % update the noise var to match if necessary
                        s.totalSignalVar_sharedTimepoints = max(0, s.totalVar_sharedTimepoints - s.totalNoiseVar_sharedTimepoints);
                        s.totalNoiseVar_sharedTimepoints =  s.totalVar_sharedTimepoints - s.totalSignalVar_sharedTimepoints;
                    else
                        s.totalSignalVar_sharedTimepoints = NaN;
                        s.totalNoiseVar_sharedTimepoints = NaN;
                    end
                    
                    % PCA explained signal variance
                    if p.Results.verbose, debug('Computing noise SVD shared timepoints\n'); end
                    if ~isempty(noiseData_NxTxC_nonNaN)
                        noiseNxT_shared = noiseData_NxTxC_nonNaN(:, :);
                        [~, Dnoise_shared, Vnoise_shared] = svd(noiseNxT_shared', 0);
                        if isvector(Dnoise_shared);
                            Snoise_shared = makecol(Dnoise_shared) * noiseVarMultiplier_sharedTimepoints;
                        else
                            Snoise_shared = diag(Dnoise_shared) * noiseVarMultiplier_sharedTimepoints;
                        end
                        if numel(Snoise_shared) < K
                            Snoise_shared = cat(1, Snoise_shared, zerosvec(K-numel(Snoise_shared)));
                        end
                        
                        if size(Vnoise_shared, 1) < K
                            Vnoise_shared = cat(2, Vnoise_shared, zeros(size(Vnoise_shared, 1), K-size(Vnoise_shared, 2)));
%                             Vnoise_shared = cat(1, Vnoise_shared, zeros(K-size(Vnoise_shared, 1), size(Vnoise_shared, 2)));
                        end
                    else
                        % not possible to compute noise
                        Snoise_shared = nanvec(maxNK);
                        Vnoise_shared = nan(maxNK, maxNK);
                    end
                        
                    pcaSignal = Sshared(1:K).^2 - Snoise_shared(1:K).^2;
                    s.pca_cumulativeSignalVarByBasis_sharedTimepoints = min(s.totalSignalVar_sharedTimepoints, cumsum(pcaSignal'));
                    s.pca_cumulativeSignalVarByBasis_sharedTimepoints = TensorUtils.makeNonDecreasing(s.pca_cumulativeSignalVarByBasis_sharedTimepoints);

                    s.pca_cumulativeFractionSignalVarByBasis_sharedTimepoints = s.pca_cumulativeSignalVarByBasis_sharedTimepoints / s.totalSignalVar_sharedTimepoints;

                    % variance explained cumulatively over basis variance
                    for i=1:K
                        s.cumulativeNoiseVarByBasis_sharedTimepoints(i) = sum(Snoise_shared(1:i).^2);
                    end
                    % explained signal var cannot exceed total, which can
                    % happen if the rates of signal and noise variance don't
                    % keep pace with each other in the first K bases.
                    s.cumulativeSignalVarByBasis_sharedTimepoints = min(s.totalSignalVar_sharedTimepoints, ...
                        s.cumulativeVarByBasis_sharedTimepoints - s.cumulativeNoiseVarByBasis_sharedTimepoints);

                    % ensure cumulative signal var is non-decreasing and adjust
                    % noise var to ensure signal+noise = total
                    s.cumulativeSignalVarByBasis_sharedTimepoints = TensorUtils.makeNonDecreasing(max(0, s.cumulativeSignalVarByBasis_sharedTimepoints));
                    s.cumulativeNoiseVarByBasis_sharedTimepoints = s.cumulativeVarByBasis_sharedTimepoints - s.cumulativeSignalVarByBasis_sharedTimepoints;

                    % compute fraction of total signal variance
                    s.cumulativeFractionSignalVarByBasis_sharedTimepoints = s.cumulativeSignalVarByBasis_sharedTimepoints / s.totalSignalVar_sharedTimepoints;

                    if s.totalSignalVar_sharedTimepoints == 0 % otherwise will be NaN
                        s.cumulativeFractionSignalVarByBasis_sharedTimepoints = 0 * s.cumulativeSignalVarByBasis_sharedTimepoints;
                        s.pca_cumulativeFractionSignalVarByBasis_sharedTimepoints = 0 * s.cumulativeSignalVarByBasis_sharedTimepoints;
                    end
                    
                    %%%%%%%%%%
                    % marginalized signal variance
                    %%%%%%%%%%

                    if p.Results.marginalize
                        % marginalize the noise tensor
                        noiseMarginalizedTensors_shared = TrialDataUtilities.DPCA.dpca_marginalize(noiseTensor_shared_maskT, ...
                            'meanSubtract', false, 'ifFlat', 'no', ... % should already be mean subtracted since it's a difference
                            'combinedParams', p.Results.combinedParams); % Theta_phi in paper

                        % total signal variance, marginalized
                        s.totalMarginalizedNoiseVar_sharedTimepoints = nan(nMarginalizations, 1);
                        noiseMarg_NxTxC = cellvec(nMarginalizations);
                        for m=1:nMarginalizations
                            noiseMarg_NxTxC{m} = noiseMarginalizedTensors_shared{m}(:, :, maskConditionsConsider);
                            s.totalMarginalizedNoiseVar_sharedTimepoints(m) = sum(noiseMarg_NxTxC{m}(:).^2) * noiseVarMultiplier_sharedTimepoints;
                        end

                        % put floor at zero and adjust noise bound to ensure they
                        % sum to total
                        s.totalMarginalizedSignalVar_sharedTimepoints = max(0, s.totalMarginalizedVar_sharedTimepoints - s.totalMarginalizedNoiseVar_sharedTimepoints);
                        s.totalMarginalizedNoiseVar_sharedTimepoints = s.totalMarginalizedVar_sharedTimepoints - s.totalMarginalizedSignalVar_sharedTimepoints; 

                        % cumulative marginalized signal variance
                        prog = ProgressBar(nMarginalizations, 'Computing signal variance over %d marginalizations', nMarginalizations);
                        %Snoise_marg = cell(nMarginalizations, 1);
                        s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints = nan(nMarginalizations, K);
                        for m=1:nMarginalizations
                            prog.update(m);

        %                     Snoise_marg{m} = svd(noiseMarginalizedTensors_shared{m}(:, :)', 0);
        %                     Snoise_marg{m} = Snoise_marg{m}(1:K);
        %                     Z = decoderKbyN*dataMarginalizedTensors_shared{m};
                            progK = ProgressBar(K, 'Computing cumulative marginalized signal variance');
                            for d = 1:K
                                progK.update(d);
                                % project marginalized noise into SVD bases found
                                % on the full noise matrix and compute the noise
                                % variance in those bases
                                s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d) = sum(sum((Vnoise_shared(:, 1:d)' * noiseMarg_NxTxC{m}(:,:)).^2)) * noiseVarMultiplier_sharedTimepoints;

        %                        s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints(m, d) = max(0, ...
        %                            sum(dataMarginalizedTensors_shared{m}(:).^2) - sum(sum((dataMarginalizedTensors_shared{m} - encoderNbyK(:,1:d)*Z(1:d,:)).^2)) - ...
        %                            s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d));

                                s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints(m, d) = ...
                                    max(0, min(s.cumulativeSignalVarByBasis_sharedTimepoints(d), ...
                                    s.cumulativeMarginalizedVarByBasis_sharedTimepoints(m, d) - ...
                                    s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d)));
                            end
                            progK.finish();
                        end
                        prog.finish();

                        % ensure cumulative signal var is non-decreasing and adjust
                        % noise var to ensure signal+noise = total
                        s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints = ...
                            TensorUtils.makeNonDecreasing(s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints, 2);
                        s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints = ...
                            s.cumulativeMarginalizedVarByBasis_sharedTimepoints - ...
                            s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints;

                        % compute fraction of total signal var
                        s.cumulativeFractionMarginalizedSignalVarByBasis_sharedTimepoints = ...
                            s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints / s.totalSignalVar_sharedTimepoints;
                    end
                end
                % end of signal, noise variance section
            end
                
            %%%%%%%%
            % Sanity checks
            %%%%%%%%
            if p.Results.verbose
                debug('Running sanity checks on explained variance\n'); 
            end
            smallMult = 1.001;
            small = 1e-6;
            checkFields(s, searchFields(s, 'total.*'), @(x) x >= -small);
            checkFields(s, searchFields(s, '.*fraction*'), @(x) x >= -small & x <= 1+small);
            checkFields(s, searchFields(s, '.*cumulativeSignal*'), @(x) isNonDecreasing(x, 2));
            checkFields(s, searchFields(s, '.*cumulativeVar*'), @(x) isNonDecreasing(x, 2));

            % cannot exceed total var
            checkFields(s, searchFields(s, '.*VarByBasis_allTimepoints'), ...
                @(x) x <= s.totalVar_allTimepoints * smallMult);
            checkFields(s, searchFields(s, '.*VarByBasis_sharedTimepoints'), ...
                @(x) x <= s.totalVar_sharedTimepoints * smallMult);

            if p.Results.marginalize
                % marginalized var total less than var total
                checkFields(s, 'totalMarginalizedVar_sharedTimepoints', @(x) x <= s.totalVar_sharedTimepoints * smallMult);

                % marginalized vars must be less than total in that marginalization
                checkFields(s, 'componentMarginalizedVarByBasis_sharedTimepoints', @(x) bsxfun(@le, x, s.totalMarginalizedVar_sharedTimepoints * smallMult));
            end
            
            if p.Results.pcaBound
                % check that appropriate fields are less than their pca counterparts
                cFields = searchFields(s, 'cumulative.*');
                pcaFields = cellfun(@(x) strcat('pca_', x), cFields, 'UniformOutput', false);
                mask = isfield(s, pcaFields);
                checkFields2(s, cFields(mask), pcaFields(mask), @(c, pc) c <= pc * smallMult); % allow small fudge for numerical error
            end
            
            % Utility functions for sanity checks
            function tf = isNonDecreasing(vec, dim)
                if nargin < 2, dim = find(size(vec) > 1, 1, 'first'); end
                delta = diff(vec, 1, dim);
                tf = all(delta(:) >= -small);
            end

            function names = searchFields(in, rex)
                fields = fieldnames(in);
                matches = ~cellfun(@isempty, regexpi(fields, rex, 'once'));
                names = fields(matches);
            end

            function checkFields2(in, fieldNames1, fieldNames2, testFn)
                if ~iscell(fieldNames1), fieldNames1 = {fieldNames1}; end
                if ~iscell(fieldNames2), fieldNames2 = {fieldNames2}; end
                for iF = 1:numel(fieldNames1)
                    %debug('Running test %s on fields %s and %s\n', func2str(testFn), fieldNames1{iF}, fieldNames2{iF});
                    tf = testFn(in.(fieldNames1{iF}), in.(fieldNames2{iF}));
                    assertWarning(all(tf(:)), 'Test %s failed on fields %s and %s', func2str(testFn), fieldNames1{iF}, fieldNames2{iF});
                end
            end

            function checkFields(in, fieldNames, testFn)
                if ~iscell(fieldNames), fieldNames = {fieldNames}; end
                for iF = 1:numel(fieldNames)
                    %debug('Running test %s on field %s\n', func2str(testFn), fieldNames{iF});
                    tf = testFn(in.(fieldNames{iF}));
                    assertWarning(all(tf(:)), 'Test %s failed on field %s', func2str(testFn), fieldNames{iF});
                end
            end
            
            function assertWarning(condition, varargin)
                if ~condition
                    warning(varargin{:});
                end
            end
        end
    end
    
    methods % for computing statistics on randomized data
        function varargout = aggregateRandomizedStatistics(s, varargin)
            % [agg1, agg2, ...] = s.aggregateRandomizedStatistics(s, prop1Name, prop2Name, ...)
            % for s.(propName), looks into s.statisticsRandomized and
            % aggregates these values along dimension 3
            
            assert(s.hasStatisticsRandomized, 'Randomized statistics were not computed in this StateSpaceProjectionStatistics');
            varargout = cell(numel(varargin), 1);
            for i = 1:numel(varargin)
                varargout{i} = cat(3, s.statisticsRandomized.(varargin{i}));
            end
        end
        
        function vals = computeQuantilesForRandomizedStatistics(s, name, quantileValues)
            data = s.aggregateRandomizedStatistics(name);
            vals = quantile(data, quantileValues, 3);
        end
    end
        
    methods(Static) % Internal static utilities 
        % below are validation functions for marginalization specifications
        function tf = isAxisSpec(x) 
            % does x describe a specific axis? scalar index or cell of strings
            tf = isempty(x) || isscalar(x) || iscellstr(x);
        end
        
        function tf = isMarginalizationSpec(x)
            % marginalizationSpec is a cell array of axesSpecs that are
            % considered one type of marginalized variance
            tf = isempty(x) || (iscell(x) && all(cellfun(@StateSpaceProjectionStatistics.isAxisSpec, x)));
        end
        
        function tf = isListOfAxisSpecs(x)
            % list of axisSpecs is a cell array of axesSpecs
            tf = StateSpaceProjectionStatistics.isMarginalizationSpec(x);
        end
        
        function tf = isMarginalizationCombination(x)
            % a marginalizationCombination is a cell array of
            % marginalization specs
            tf = isempty(x) || (iscell(x) && all(cellfun(@StateSpaceProjectionStatistics.isMarginalizationSpec, x)));
        end
        
        function tf = isListOfMarginalizationSpec(x)
            % cell array of marginalizationSpec
            tf = StateSpaceProjectionStatistics.isMarginalizationCombination(x);
        end
        
        function tf = isListOfMarginalizationCombinations(x)
            % list of marginalizationCombinations is a cell array of
            % marginalizationCombinations
            tf =  isempty(x) || (iscell(x) && all(cellfun(@StateSpaceProjectionStatistics.isMarginalizationCombination, x)));
        end
        
        % translate axisSpec cellstrs into numeric axis numbers (one of the numbers in dims)
        function axisIndex = staticAxisLookup(axisSpec, axisAttributeSets, axisInclude)
            axisInclude = TensorUtils.vectorMaskToIndices(axisInclude);
            if isempty(axisSpec)
                axisIndex = [];
                return
            end

            if isnumeric(axisSpec)
                assert(ismember(axisSpec, axisInclude), 'axisSpec %s not found in dims [%s]', num2str(axisSpec), vec2str(axisInclude));
                axisIndex = axisSpec;
            else
                axisSpec = sort(axisSpec);
                if numel(axisSpec) == 1 && strcmp(axisSpec{1}, 'time')
                    axisIndex = 0;
                    return;
                end
                for iA = 1:numel(axisAttributeSets)
                    if isequal(sort(axisAttributeSets{iA}), axisSpec)
                        axisIndex = iA;
                        return;
                    end
                end
                error('axisSpec {%s} not found in axisAttributeSets', strjoin(axisSpec, ','));
            end
        end
        
        function axisIndices = staticLookupMarginalizationSpec(axisSpecList, axisAttributeSets, axisInclude)
            assert(StateSpaceProjectionStatistics.isMarginalizationSpec(axisSpecList), 'AxisSpecList must be marginalizationSpec');
            if iscell(axisSpecList)
                axisIndices = makecol(cellfun(@(x) StateSpaceProjectionStatistics.staticLookupAxis(x, axisAttributeSets, axisInclude), axisSpecList));
            else
                axisIndices = makecol(arrayfun(@(x) StateSpaceProjectionStatistics.staticLookupAxis(x, axisAttributeSets, axisInclude), axisSpecList));
            end
        end
        
        function listOfAxisIndexVectors = staticLookupMarginalizationCombination(margComb, axisAttributeSets, axisInclude)
            assert(StateSpaceProjectionStatistics.isMarginalizationCombination(margComb), 'Input must be marginalizationCombination');
            listOfAxisIndexVectors = makecol(cellfun(@(x) StateSpaceProjectionStatistics.staticLookupMarginalizationSpec(x, axisAttributeSets, axisInclude), margComb, 'UniformOutput', false));
        end
        
         function listOfListsOfAxisIndexVectors = staticLookupListOfMarginalizationCombinations(listMargComb, axisAttributeSets, axisInclude)
             assert(StateSpaceProjectionStatistics.isListOfMarginalizationCombinations(listMargComb), 'Input must be list of marginalizationCombinations');
             listOfListsOfAxisIndexVectors = makecol(cellfun(@(x) StateSpaceProjectionStatistics.staticLookupListOfMarginalizationCombinations(x, axisAttributeSets, axisInclude), listMargComb, 'UniformOutput', false));
         end

        function [fullList, names, axisIncludeList] = generateCombinedParamsForMarginalization(axisAttributeSets, varargin)
            % function [fullList, names] = generateCombinedParamsForMarginalization(dims, varargin)
            % builds combinedParams argument for dpca, which specifies the covariate and
            % covariate interactions to consider collectively.
            %
            % In the parameters below, an axisSpec indicates a specific axis. This is
            % done either via scalar numeric index (0 = time, 1 = first axis, etc.) or
            % by specifying the set of axisAttributes that comprise that axis as a
            % cellstr. Lone strings is not accepted, it must be wrapped in { } to 
            % resolve ambiguities. E.g. {'axis1Attr'} or {'axis1AttrA', 'axis1AttrB'}
            % or [1].
            %
            % An marginalizationSpec is a cell list of axisSpecs, which therefore indicates a
            % subset of the list of axes. This indicates a subset of axes, referrring
            % to covariation that arises from the interactions between the individual
            % axes included, e.g. { {'time'}, {'stimulus'} } or { 0, 1 }, would refer to
            % variance resulting from the interaction of time and stimulus (axis 1).
            %
            % A marginalizationCombination is a list of marginalizationSpecs. This is
            % used to indicate a set of marginalizationSpecs that we wish to consider
            % collectively in further analysis, that is, we will consider collectively the
            % variance arising due each marginalizationSpec alone. For example, 
            % { {{'time'}}, {{'time'}, {'stimulus'}} }, or { {0}, {0, 1} } would
            % indicate that we wish to collectively consider pure-time variance with
            % variance arising from interactions of time and stimulus. Note that this
            % is different from { {{'time'}}, {{'stimulus'}} } or { {0}, {1} }, which
            % would consider pure-time and pure-stimulus variance collectively, but
            % would not reserve time-stimulus interaction as a different
            % marginalization.
            % 
            % INPUTS:
            %
            % axisAttributeSets: cellvec of axisSpec, each indicating the set of axis
            %   attributes in each axis, e.g. { {'stimulus'}, {'targetDirection', 'targetDistance'} }
            %
            % axisMask: logical vector masking axes to marginalize over,
            %   Must match length of axisAttributeSets, e.g. [true; true]
            %
            % 'axesCombineSpecificMarginalizations': a list of marginalization
            %   combinations to apply.
            %
            % 'axesCombineAllMarginalizations': a cellvec of cellvec of axisSpec. For each
            %   cellvec of axes in the list, the appropriate marginalizationCombinations
            %   will be applied to prevent the marginalization from ever distinguishing
            %   variance due to one axis in the list from the others.
            %
            % 'combineAxesWithTime': logical scalar or logical vector or
            %   cell list of axis specs, indicating which axes should be
            %   combined with pure time (e.g. axis and axis+time will be
            %   combined)
            %
            % 'axisNames': nCov cellstr. optional, specifies the names for each
            %   covariate. should match length of axisAttributeSets, i.e.
            %   subselected for dims
            %

            import(getPackageImportString);

            p = inputParser();
            p.addParameter('axisIncludeMask', [], @(x) isempty(x) || islogical(x));
            p.addParameter('axesIgnore', {}, @(x) true);
            p.addParameter('axesCombineSpecificMarginalizations', {}, @(x) true);
            p.addParameter('axesCombineAllMarginalizations', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('combineAxesWithTime', true, @(x) islogical(x) || iscell(x));
            p.addParameter('axisNames', {}, @iscellstr);
            p.parse(varargin{:});

            axesCombineSpecificMarginalizations = p.Results.axesCombineSpecificMarginalizations;
            axesCombineAllMarginalizations = p.Results.axesCombineAllMarginalizations;

            assert(StateSpaceProjectionStatistics.isListOfAxisSpecs(axisAttributeSets), 'axisAttributeSets must be list of axisSpecs, see help.');
            axesIgnore = p.Results.axesIgnore;
            assert(StateSpaceProjectionStatistics.isListOfAxisSpecs(axesIgnore), 'axisIgnore must be list of axisSpecs, see help.');
            
            assert(StateSpaceProjectionStatistics.isListOfMarginalizationCombinations(axesCombineSpecificMarginalizations), ...
                'axesCombineSpecificMarginalizations must be a list of marginalizationCombinations, see help');
            % this isnt the right nomenclature since we're specifying list of lists
            % of axisSpec, rather than precise marginalizations, but the test is
            % equivalent
            assert(StateSpaceProjectionStatistics.isMarginalizationCombination(axesCombineAllMarginalizations), ...
                'axesCombineAllMarginalizations must be a list of list of axisSpecs, see help');

            combineAxesWithTime = p.Results.combineAxesWithTime;
            assert(islogical(combineAxesWithTime) || StateSpaceProjectionStatistics.isListOfAxisSpecs(combineAxesWithTime), 'combineAxesWithTime must be logical or list of axisSpecs, see help');

            % utility function for printing lists
            function pr(axisSpecList) %#ok<DEFNU>
                for iA = 1:numel(axisSpecList)
                    fprintf('%s\n', strjoin(cellfun(@(x) strjoin(x, ','), axisSpecList{iA}, 'UniformOutput', false), ' / '));
                end
            end

            % build dims referencing axesMask
            dims = 1:numel(axisAttributeSets);
            if ~isempty(p.Results.axisIncludeMask)
                dims = dims(p.Results.axisIncludeMask);
            end
            
            % then drop axes listed in axesIgnore
            axesIgnore = StateSpaceProjectionStatistics.staticLookupMarginalizationSpec(axesIgnore, axisAttributeSets, dims);
            dims = setdiff(dims, axesIgnore);
            axisIncludeList = dims;
            
            combineSpecific = StateSpaceProjectionStatistics.staticLookupMarginalizationCombination(axesCombineSpecificMarginalizations, axisAttributeSets, axisIncludeList);
%             combineSpecific = map(@(x) map(@(x) doAxisLookupMultiple(x, axisAttributeSets, dims), axesCombineSpecificMarginalizations));
            
            combineAll = StateSpaceProjectionStatistics.staticLookupListOfMarginalizationCombinations(axesCombineAllMarginalizations, axisAttributeSets, axisIncludeList);
%            combineAll = map(@doAxisLookupMultiple, axesCombineAllMarginalizations);

            if isempty(combineSpecific)
                combineSpecific = cell(0, 1); % to allow for vertcat later
            end

            % generate names for each dim if not provided
            nDims = numel(dims);
            if isempty(p.Results.axisNames)
                covariateNames = arrayfun(@(i) sprintf('Axis %d', i), ...
                    1:numel(axisAttributeSets), 'UniformOutput', false);
            else
                covariateNames = p.Results.axisNames;
                assert(numel(covariateNames) >= nDims, 'Provided covariateNames needs at least %d entries', nDims);
            end
            covariateNames = makecol(covariateNames);

            % add dim, time+dim combinations to the list 
            if iscell(combineAxesWithTime)
                combineAxesWithTime = StateSpaceProjectionStatistics.staticLookupMarginalizationSpec(combineAxesWithTime, axisAttributeSets, axisIncludeList);
%                 combineAxesWithTime = doAxisLookupMultiple(combineAxesWithTime);
            end
            if islogical(combineAxesWithTime) && isscalar(combineAxesWithTime)
                combineAxesWithTime = repmat(combineAxesWithTime, numel(dims), 1);
            end
            list = supersets(num2cell(dims(combineAxesWithTime)));
            % add all the {0, each-of-supersets(d)} combinations
            timeList = cellfun(@(dim) {dim; [0; dim]}, list, 'UniformOutput', false);
            combineSpecific = [combineSpecific; timeList];

            % generate full list
            fullList = num2cell(subsets(union(0, dims)));
            
            % merge axis all combinations to the list
            for iSet = 1:numel(combineAll)
                set = combineAll{iSet}; % numberic vector of axes   
                theseAxesCombos = subsets(set);
                combineSpecific = [combineSpecific; {theseAxesCombos}]; %#ok<AGROW>
                
                otherAxisCombos = subsets(setdiff([0, dims], set));
                for iOther = 1:numel(otherAxisCombos)
                    for iThis = 1:numel(theseAxesCombos)
                        withSubset = union(theseAxesCombos{iThis}, otherAxisCombos{iOther});
                        withAll = union(set, otherAxisCombos{iOther});
                        combineSpecific = [combineSpecific; {{withSubset; withAll}}]; %#ok<AGROW>
                    end
                end
            end

            % remove singular lists (nothing being combined)
            nToCombine = cellfun(@numel, combineSpecific);
            combineSpecific = combineSpecific(nToCombine > 1);

           
            for iC = 1:numel(combineSpecific)
                % search for any superset of the elements of the combine list
                combine = combineSpecific{iC};
                searchFor = supersets(combine);

                % search for rows of full list that match
                %isSubset = @(x, of) all(ismember(x, of));
                matches = falsevec(numel(fullList));
                for iV = 1:numel(searchFor)
                    value = searchFor{iV};
                    rowMatchesFn = @(list) any(cellfun(@(v) isequal(v, value), list));
                    matches = matches | cellfun(rowMatchesFn, fullList);
                end

                % combine all the matches
                combinedMatches = uniqueCell(cat(1, fullList {matches}));

                % remove the matching rows and add the combined row
                fullList = [fullList(~matches); {combinedMatches}];
            end
        %     pr(fullList)

            % generate marginalization names
            if nargout > 1
                covariateNames = [{'time'}; covariateNames];
                names = cellvec(numel(fullList));
                for iF = 1:numel(fullList)
                    pieceStr = cellvec(numel(fullList{iF}));
                    for iPiece = 1:numel(fullList{iF})
                        index = fullList{iF}{iPiece} + 1;
                        % if we're combining with time automatically, don't
                        % include time in the marginalization names
                        if any(combineAxesWithTime) && ~isequal(index, 1) && ismember(1, index)
                            pieceStr{iPiece} = '';
                        else
                            pieceStr{iPiece} = strjoin(covariateNames(index), ' x ');
                        end
                    end
                    pieceStr = pieceStr(~cellfun(@isempty, pieceStr));
                    names{iF} = strjoin(pieceStr, ', ');
                end
            end

            % change to being time = 1, dim 1 = 2 indexed
            for iF = 1:numel(fullList)
                for iJ = 1:numel(fullList{iF})
                    fullList{iF}{iJ} = fullList{iF}{iJ} + 1;
                end
            end

            function S = subsets(X)

                % S = subsets(X) returns a cell array of all subsets of vector X apart
                % from the empty set. Subsets are ordered by the number of elements in
                % ascending order.
                %
                % subset([1 2 3]) = {[1], [2], [3], [1 2], [1 3], [2 3], [1 2 3]}

                X = makecol(X);    
                d = length(X);
                pc = dec2bin(1:2^d-1) - '0';
                [~, ind] = sort(sum(pc, 2));
                pc = fliplr(pc(ind,:));
                S = cellvec(length(pc));
                for iP=1:length(pc)
                    S{iP} = makecol(X(logical(pc(iP,:))));
                end
            end

            function S = supersets(X)
                % S = supersets(X) returns a cell array of all supersets of the cell contents of X apart
                % from the empty set.
                %
                % superset({1, [1 2], [1 3]}) = {[1], [1 2], [1 3], [1 2 3]}

                X = cellfun(@makecol, X, 'UniformOutput', false);
                idxSubsets = subsets(1:numel(X));

                S = cellvec(numel(idxSubsets));
                for iS = 1:numel(idxSubsets)
                    S{iS} = unique(cat(1, X{idxSubsets{iS}}));
                end

                S = uniqueCell(S);
            end

            function [B, I, J] = uniqueCell(A)
                B = cell(0, 1);
                I = [];
                J = zeros(numel(A),1);
                for iA = 1:numel(A)
                    idx = find(cellfun(@(b) isequal(A{iA}, b), B), 1, 'first');
                    if isempty(idx)
                        B{end+1} = A{iA}; %#ok<AGROW>
                        I(end+1) = iA; %#ok<AGROW>
                        J(iA) = numel(B); 
                    else
                        J(iA) = idx;
                    end
                end
                B = makecol(B);
            end 
        end
    end
end
