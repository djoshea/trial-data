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

        conditionDescriptor
        marginalizationNames
        combinedParams
        axisCombinations
        combineCovariatesWithTime = true;
        
        marginalizationColorsStored
    end
    
    properties(Dependent)
        marginalizationColors
        hasSignalVariance
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
    end
    
    methods(Access=protected)
        function s = StateSpaceProjectionStatistics()
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
    end

    methods
        function plotCumulativeVariance(s, varargin)
            p = inputParser();
            p.addParameter('fractional', true, @islogical);
            p.addParameter('showPCABound', true, @islogical);
            p.addParameter('timepoints', 'all', @ischar); % or 'shared'
            p.addParameter('varianceType', 'total', @ischar);
            p.addParameter('marginalize', [], @(x) isempty(x) || islogical(x));
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
                    pcaFieldName = '';
                end
            else
                pcaFieldName = '';
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
            au.addColoredLabels(s.marginalizationNames, s.marginalizationColors, ...
                'posY', AutoAxis.PositionType.Top, 'posX', AutoAxis.PositionType.Left);
            au.update();
            
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
            
%             legend(hPatch, nCov, s.covMarginalizedNames, 'Location', 'NorthEastOutside');
%             legend boxoff;
            hold off;
        end
    end
    
    methods(Static) % compute statistics!
        function s = build(proj, pset, varargin)
            p = inputParser();
            p.addParameter('axisCombinations', {}, @iscell);
            p.addParameter('combineCovariatesWithTime', true, @islogical);
            p.parse(varargin{:});
            
            s = StateSpaceProjectionStatistics();
            
            % get the names and counts of source -> proj bases
            s.basisNamesSource = pset.basisNames;
            s.basisNamesProj = proj.getBasisNames(pset);
            s.nBasesSource = proj.nBasesSource;
            s.nBasesProj = proj.nBasesProj;
            s.conditionDescriptor = pset.conditionDescriptor;
                      
            % filter for non-singular axes
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            dimIdx = find(dimMask);
            
            % merge each covariate with each covariate + time mixture
            s.axisCombinations = p.Results.axisCombinations;
            axisListsToCombine = cellvec(numel(s.axisCombinations));
            for i = 1:numel(s.axisCombinations)
                list = s.axisCombinations{i};
                axisListsToCombine{i} = cellvec(numel(list));
                for j = 1:numel(list)
                    if strcmp(list{j}, 'time') || isequal(list{j}, 0)
                        axisListsToCombine{i}{j} = 0;
                    else
                        axisListsToCombine{i}{j} = pset.conditionDescriptor.axisLookupByAttributes(list{j});
                    end
                end
            end
            
            % build the list of covariates and covariate interactions to marginalize along
            [s.combinedParams, s.marginalizationNames] = TrialDataUtilities.DPCA.dpca_generateTimeCombinedParams(...
                dimIdx, 'covariateNames', pset.conditionDescriptor.axisNames(:), ...
                'combineEachWithTime', p.Results.combineCovariatesWithTime, ...
                'combine', axisListsToCombine); %#ok<FNDSB>
                   
            % Extract trial averaged data
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            
            % Save covariance matrices
            CTAbyNv = pset.buildCTAbyN('validBasesOnly', true);
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
                    'scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr', scaledNoiseEstimate_NvbyTAbyAttr);
            else
                s = s.computeStatistics(NvbyTAbyAttr, ...
                    decoderKbyNv, encoderNvbyK, 'combinedParams', s.combinedParams);
            end
        end
    end
    
    methods(Access=protected)
        function s = computeStatistics(s, NbyTAbyAttr, decoderKbyN, encoderNbyK, varargin)
            % based heavily on dpca_explainedVariance, though modified
            % heavily
            %
            % stats = computeProjectionStatistics(NbyTAbyAttr, decoderKbyN, encoderNbyK, varargin) 
            % computes various measures and stores them in
            % dataNbyT is the data
            % matrix, decoder is the decoder matrix (K x N), encoder is the encoder matrix (N x K).
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
            assert(size(NbyTAbyAttr, 1) == N, 'Shape of decoder must be K by N');
            assert(size(encoderNbyK, 1) == N, 'Shape of encoder must be N by K');
            assert(size(decoderKbyN, 1) == size(encoderNbyK, 2), 'Encoder and decoder differ on size K');
            
            p = inputParser();
            p.addParameter('combinedParams', {}, @(x) true);
            p.addParameter('scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr', [], @(x) isempty(x) || isnumeric(x));
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
            allTimepoints_NxT_keepMaskT = ~any(isnan(dataNbyT), 1);
            dataNbyT_allTimepoints = dataNbyT(:, allTimepoints_NxT_keepMaskT);
            
            % subtract means for "all timepoints" over time for each basis,
            % condition
            dataNbyT_allTimepoints = bsxfun(@minus, dataNbyT_allTimepoints, mean(dataNbyT_allTimepoints, 2));

            % total variance
            s.totalVar_allTimepoints = sum(dataNbyT_allTimepoints(:).^2);

            % PCA explained variance
            debug('Computing trial-average SVD\n');
            Spca_all = svd(dataNbyT_allTimepoints', 0);
            s.pca_cumulativeVarByBasis_allTimepoints = cumsum(Spca_all.^2');
            s.pca_cumulativeFractionVarByBasis_allTimepoints = s.pca_cumulativeVarByBasis_allTimepoints / s.totalVar_allTimepoints;
            
            % explained variance along cumulative sets of bases
            Z = decoderKbyN*dataNbyT_allTimepoints;
            s.componentVarByBasis_allTimepoints = sum(Z.^2, 2)';
            
            for i=1:size(decoderKbyN,1)
                s.cumulativeVarByBasis_allTimepoints(i) = s.totalVar_allTimepoints - sum(sum((dataNbyT_allTimepoints - encoderNbyK(:,1:i)*Z(1:i,:)).^2));    
            end
            s.cumulativeFractionVarByBasis_allTimepoints = s.cumulativeVarByBasis_allTimepoints / s.totalVar_allTimepoints;

            %%%%%%%
            % Shared timepoints variance
            %%%%%%%
            
            % find timepoints where all conditions and all neurons have samples, any
            % along non-time dimensions
            sharedTimepoints_tensor_keepMaskT = squeeze(~TensorUtils.anyMultiDim(isnan(dataTensor), [1 3:ndims(dataTensor)]));
            dataTensor_shared = TensorUtils.selectAlongDimension(dataTensor, 2, sharedTimepoints_tensor_keepMaskT);

            % subtract means for "all timepoints" over time for each basis,
            % condition. do this once for the tensor
            dataTensor_shared = TensorUtils.centerSlicesOrthogonalToDimension(dataTensor_shared, 1);

            % PCA explained variance upper bound
            debug('Computing trial-average SVD shared timepoints\n');
            dataNxT_shared = dataTensor_shared(:, :);
            Sshared = svd(dataNxT_shared', 0);
            s.pca_cumulativeVarByBasis_sharedTimepoints = cumsum(Sshared(1:K).^2');

            s.totalVar_sharedTimepoints = sum(dataNxT_shared(:).^2);
            s.pca_cumulativeFractionVarByBasis_sharedTimepoints = s.pca_cumulativeVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;

            % compute total variance separately for the shared timepoints only
            s.totalVar_sharedTimepoints = sum(dataTensor_shared(:).^2);

            % explained variance along bases
            Z = decoderKbyN*dataNxT_shared;
            s.componentVarByBasis_sharedTimepoints = sum(Z.^2, 2)'; % along each basis INDIVIDUALLY
            for i=1:K
                s.cumulativeVarByBasis_sharedTimepoints(i) = s.totalVar_sharedTimepoints - sum(sum((dataNxT_shared - encoderNbyK(:,1:i)*Z(1:i,:)).^2));    
            end
            s.cumulativeFractionVarByBasis_sharedTimepoints = s.cumulativeVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;

            %%%%%%
            % Marginalized Variances
            %%%%%%
            
            % for marginalization though, we need to only use timepoints where all conditions
            % are present (non-nan). This means the total variance over all
            % marginalizations will be different than for the non marginalized
            % variance, since we throw away timepoints not shared across all conditions
            % hence the _sharedTimepoints vs allTimepoints distinction

            % marginalizing the trimmed tensor, should already be
            % mean-subtracted
            dataMarginalizedTensors_shared = TrialDataUtilities.DPCA.dpca_marginalize(dataTensor_shared, ...
                'meanSubtract', false, ... % should already be mean subtracted since it's a difference
                'combinedParams', p.Results.combinedParams, 'ifFlat', 'yes');
            nMarginalizations = length(dataMarginalizedTensors_shared);
            
            Zshared = decoderKbyN*dataNxT_shared;
            s.componentVarByBasis_allTimepoints = sum(Zshared.^2, 2)';
            
            % total marginalized variance
            s.totalMarginalizedVar_sharedTimepoints = nan(nMarginalizations, 1);
            for i=1:nMarginalizations
                s.totalMarginalizedVar_sharedTimepoints(i) = sum(dataMarginalizedTensors_shared{i}(:).^2);
            end

            % marginalized variance of each component : shared timepoints
            for i=1:nMarginalizations
                s.componentMarginalizedVarByBasis_sharedTimepoints(i,:) = sum((decoderKbyN * dataMarginalizedTensors_shared{i}).^2, 2)'; % @ djoshea shouldn't be normalized
            end
            
            % explained variance along cumulative sets of bases
            for m=1:nMarginalizations
                thisMarg = dataMarginalizedTensors_shared{m}(:, :);
                totalVarThisMarg = sum(thisMarg(:).^2);
                Zmarg = decoderKbyN*thisMarg;
                for i=1:K
                    s.cumulativeMarginalizedVarByBasis_sharedTimepoints(m, i) = totalVarThisMarg - sum(sum((thisMarg - encoderNbyK(:,1:i)*Zmarg(1:i,:)).^2));    
                end
            end
            s.cumulativeFractionMarginalizedVarByBasis_sharedTimepoints = s.cumulativeMarginalizedVarByBasis_sharedTimepoints / s.totalVar_sharedTimepoints;

            %%%%%%%
            % Noise and signal variances, if single trial data provided
            %%%%%%%
            if ~isempty(p.Results.scaledDifferenceOfTrialsNoiseEstimate_NbyTAbyAttr)
                debug('Computing signal variance via noise-floor\n');
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
                if any(isnan(noiseNxT_all(:)))
                    error('Single trial data has NaN values at timepoints where trial-averaged data did not');
                end

%                 % we mean subtract here because we will also mean subtract
%                 % when we look at noise in shared timepoints, and we want
%                 % the noise variances to be comparable when all timepoints
%                 % are shared timepoints. we mean subtract there because the
%                 % marginalization essentially accomplishes the same effect,
%                 % and then we'll 
%                 noiseNxT_all = bsxfun(@minus, noiseNxT_all, mean(noiseNxT_all, 2));

                %%%%%%
                % Noise variance, all timepoints
                %%%%%%
                
                 % total noise variance
                s.totalNoiseVar_allTimepoints = sum(noiseNxT_all(:).^2);

                % total signal variance
                s.totalSignalVar_allTimepoints = s.totalVar_allTimepoints - s.totalNoiseVar_allTimepoints;

                % PCA explained signal variance
                debug('Computing noise SVD\n');
                Snoise_all = svd(noiseNxT_all', 0);
                pcaSignal = Spca_all(1:K).^2 - Snoise_all(1:K).^2;
                s.pca_cumulativeSignalVarByBasis_allTimepoints = cumsum(pcaSignal');
                s.pca_cumulativeSignalVarByBasis_allTimepoints = TensorUtils.makeNonDecreasing(s.pca_cumulativeSignalVarByBasis_allTimepoints);
                
                s.pca_cumulativeFractionSignalVarByBasis_allTimepoints = s.pca_cumulativeSignalVarByBasis_allTimepoints / s.totalSignalVar_allTimepoints;
                
                
                % variance explained cumulatively over basis variance
                for i=1:K
                    s.cumulativeNoiseVarByBasis_allTimepoints(i) = sum(Snoise_all(1:i).^2);
                end
                s.cumulativeSignalVarByBasis_allTimepoints = s.cumulativeVarByBasis_allTimepoints - s.cumulativeNoiseVarByBasis_allTimepoints;
                
                % ensure cumulative signal var is non-decreasing and adjust
                % noise var to ensure signal+noise = total
                s.cumulativeSignalVarByBasis_allTimepoints = TensorUtils.makeNonDecreasing(s.cumulativeSignalVarByBasis_allTimepoints);
                s.cumulativeNoiseVarByBasis_allTimepoints = s.cumulativeVarByBasis_allTimepoints - s.cumulativeSignalVarByBasis_allTimepoints;
                
                % compute fraction of total signal var
                s.cumulativeFractionSignalVarByBasis_allTimepoints = s.cumulativeSignalVarByBasis_allTimepoints / s.totalSignalVar_allTimepoints;

                %%%%%%%%%%
                % Noise variance, shared timepoints
                %%%%%%%%%%

                % filter in time to match dataTensor_shared and hopefully remove all NaNs 
                % in shared timepoints
                noiseTensor_shared = TensorUtils.selectAlongDimension(noiseTensor, 2, ...
                    sharedTimepoints_tensor_keepMaskT);
%                 % subtract mean by basis now
%                 noiseTensor_shared = bsxfun(@minus, noiseTensor_shared, nanmean(noiseTensor_shared,2));

                % this may need to be replaced with something better, there's no
                % guarantee that individual trials won't have NaNs at places where the
                % trial-averaged traces do not, since there could still be enough
                % trials at that timepoint to form a non-NaN average. However, if we
                % drop additional timepoints here, then we'll have less total variance
                % to explain and need to correct for this. #todo
                if any(isnan(noiseTensor_shared(:)))
                    error('Single trial data has NaN values at timepoints where trial-averaged data did not');
                end

                s.totalNoiseVar_sharedTimepoints = sum(noiseTensor_shared(:).^2);
                
                % ensure signal var is nonnegative using max(0, ...) and
                % update the noise var to match if necessary
                s.totalSignalVar_sharedTimepoints = max(0, s.totalVar_sharedTimepoints - s.totalNoiseVar_sharedTimepoints);
                s.totalNoiseVar_sharedTimepoints =  s.totalVar_sharedTimepoints - s.totalSignalVar_sharedTimepoints;

                % PCA explained signal variance
                debug('Computing noise SVD shared timepoints\n');
                noiseNxT_shared = noiseTensor_shared(:, :);
                [~, Dnoise_shared, Vnoise_shared] = svd(noiseNxT_shared', 0);
                Snoise_shared = diag(Dnoise_shared);
                pcaSignal = Sshared(1:K).^2 - Snoise_shared(1:K).^2;
                s.pca_cumulativeSignalVarByBasis_sharedTimepoints = cumsum(pcaSignal');
                s.pca_cumulativeSignalVarByBasis_sharedTimepoints = TensorUtils.makeNonDecreasing(s.pca_cumulativeSignalVarByBasis_sharedTimepoints);
                
                s.pca_cumulativeFractionSignalVarByBasis_sharedTimepoints = s.pca_cumulativeSignalVarByBasis_sharedTimepoints / s.totalSignalVar_sharedTimepoints;

                % variance explained cumulatively over basis variance
                for i=1:K
                    s.cumulativeNoiseVarByBasis_sharedTimepoints(i) = sum(Snoise_shared(1:i).^2);
                end
                s.cumulativeSignalVarByBasis_sharedTimepoints = s.cumulativeVarByBasis_sharedTimepoints - s.cumulativeNoiseVarByBasis_sharedTimepoints;
                
                % ensure cumulative signal var is non-decreasing and adjust
                % noise var to ensure signal+noise = total
                s.cumulativeSignalVarByBasis_sharedTimepoints = TensorUtils.makeNonDecreasing(s.cumulativeSignalVarByBasis_sharedTimepoints);
                s.cumulativeNoiseVarByBasis_sharedTimepoints = s.cumulativeVarByBasis_sharedTimepoints - s.cumulativeSignalVarByBasis_sharedTimepoints;
                                
                % compute fraction of total signal variance
                s.cumulativeFractionSignalVarByBasis_sharedTimepoints = s.cumulativeSignalVarByBasis_sharedTimepoints / s.totalSignalVar_sharedTimepoints;

                %%%%%%%%%%
                % marginalized signal variance
                %%%%%%%%%%

                % marginalize the noise tensor
                noiseMarginalizedTensors_shared = TrialDataUtilities.DPCA.dpca_marginalize(noiseTensor_shared, ...
                    'meanSubtract', false, ... % should already be mean subtracted since it's a difference
                    'combinedParams', p.Results.combinedParams); % Theta_phi in paper

                % total signal variance, marginalized
                s.totalMarginalizedNoiseVar_sharedTimepoints = nan(nMarginalizations, 1);
                for m=1:nMarginalizations
                    s.totalMarginalizedNoiseVar_sharedTimepoints(m) = sum(noiseMarginalizedTensors_shared{m}(:).^2);
                end
                
                % put floor at zero and adjust noise bound to ensure they
                % sum to total
                s.totalMarginalizedSignalVar_sharedTimepoints = max(0, s.totalMarginalizedVar_sharedTimepoints - s.totalMarginalizedNoiseVar_sharedTimepoints);
                s.totalMarginalizedNoiseVar_sharedTimepoints = s.totalMarginalizedVar_sharedTimepoints - s.totalMarginalizedSignalVar_sharedTimepoints; 
                
                % cumulative marginalized signal variance
                %prog = ProgressBar(nMarginalizations, 'Computing svg for noise marginalizations');
                %Snoise_marg = cell(nMarginalizations, 1);
                s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints = nan(nMarginalizations, K);
                for m=1:nMarginalizations
                    %prog.update(m);

%                     Snoise_marg{m} = svd(noiseMarginalizedTensors_shared{m}(:, :)', 0);
%                     Snoise_marg{m} = Snoise_marg{m}(1:K);
                    Z = decoderKbyN*dataMarginalizedTensors_shared{m};
                    for d = 1:K
                        % project marginalized noise into SVD bases found
                        % on the full noise matrix and compute the noise
                        % variance in those bases
                        s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d) = sum(sum((Vnoise_shared(:, 1:d)' * noiseMarginalizedTensors_shared{m}(:,:)).^2));
                       
%                        s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints(m, d) = max(0, ...
%                            sum(dataMarginalizedTensors_shared{m}(:).^2) - sum(sum((dataMarginalizedTensors_shared{m} - encoderNbyK(:,1:d)*Z(1:d,:)).^2)) - ...
%                            s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d));

                        s.cumulativeMarginalizedSignalVarByBasis_sharedTimepoints(m, d) = ...
                            max(0, min(s.cumulativeSignalVarByBasis_sharedTimepoints(d), ...
                            s.cumulativeMarginalizedVarByBasis_sharedTimepoints(m, d) - ...
                            s.cumulativeMarginalizedNoiseVarByBasis_sharedTimepoints(m, d)));
                    end
                end
                %prog.finish();
                
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
            % end of signal, noise variance section
            
            %%%%%%%%
            % Sanity checks
            %%%%%%%%
            debug('Running sanity checks on explained variance\n');
            smallMult = 1.001;
            small = 1e-08;
            checkFields(s, searchFields(s, 'total.*'), @(x) x >= -small);
            checkFields(s, searchFields(s, '.*fraction*'), @(x) x >= -small & x <= 1+small);
            checkFields(s, searchFields(s, '.*cumulativeSignal*'), @(x) isNonDecreasing(x, 2));
            checkFields(s, searchFields(s, '.*cumulativeVar*'), @(x) isNonDecreasing(x, 2));

            % cannot exceed total var
            checkFields(s, searchFields(s, '.*VarByBasis_allTimepoints'), ...
                @(x) x <= s.totalVar_allTimepoints * smallMult);
            checkFields(s, searchFields(s, '.*VarByBasis_sharedTimepoints'), ...
                @(x) x <= s.totalVar_sharedTimepoints * smallMult);

            % marginalized var total less than var total
            checkFields(s, 'totalMarginalizedVar_sharedTimepoints', @(x) x <= s.totalVar_sharedTimepoints * smallMult);

            % marginalized vars must be less than total in that marginalization
            checkFields(s, 'componentMarginalizedVarByBasis_sharedTimepoints', @(x) bsxfun(@le, x, s.totalMarginalizedVar_sharedTimepoints * smallMult));

            % check that appropriate fields are less than their pca counterparts
            cFields = searchFields(s, 'cumulative.*');
            pcaFields = cellfun(@(x) strcat('pca_', x), cFields, 'UniformOutput', false);
            mask = isfield(s, pcaFields);
            checkFields2(s, cFields(mask), pcaFields(mask), @(c, pc) c <= pc * smallMult); % allow small fudge for numerical error

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
                    assert(all(tf(:)), 'Test %s failed on fields %s and %s', func2str(testFn), fieldNames1{iF}, fieldNames2{iF});
                end
            end

            function checkFields(in, fieldNames, testFn)
                if ~iscell(fieldNames), fieldNames = {fieldNames}; end
                for iF = 1:numel(fieldNames)
                    %debug('Running test %s on field %s\n', func2str(testFn), fieldNames{iF});
                    tf = testFn(in.(fieldNames{iF}));
                    assert(all(tf(:)), 'Test %s failed on field %s', func2str(testFn), fieldNames{iF});
                end
            end
        end


    end
end
