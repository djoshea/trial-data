classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        lambda = [];
        maxTrialsForOptimizeLambda = 50;
        
        nBasesProjKeep = NaN;
        nBasesProjPerMarginalization = NaN; % the number of bases per marginalization to find, before selecting the top nBasesProj. this can be a vector, one for each marginalization type too
        
        nIterationsOptimizeLambda = 5; % higher is more accurate, 5-10 is reasonable
        
        % the properties below control how different types of variance
        % arising along the Condition axes will be treated by DPCA. These
        % have the same meaning as in StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization
        % and you can refer to the documentation there for more details
        
        % list of specific marginalization combinations to combine
        axesCombineSpecificMarginalizations = {};
            
        % cellvec of cellvec of axisSpec. For each
        % cellvec of axes in the list, the appropriate marginalizationCombinations
        % will be applied to prevent the marginalization from ever distinguishing
        % variance due to one axis in the list from the others.
        axesCombineAllMarginalizations = {};
            
        % logical scalar or logical vector or
        % cell list of axis specs, indicating which axes should be
        % combined with pure time (e.g. axis and axis+time will be
        % combined)
        combineAxesWithTime = true; 
    end
    
    properties(Constant)
        defaultLambda = 0;
    end
    
    properties(SetAccess=protected)
        %  'combinedParams' - cell array of cell arrays specifying 
        %                     which marginalizations should be added up together,
        %                     e.g. for the three-parameter case with parameters
        %                           1: stimulus
        %                           2: decision
        %                           3: time
        %                     one could use the following value:
        %                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.
        combinedParams
        
        marginalizationNames
        
        marginalizationIdxByBasis
        
        dataMarginalized
        dataMarginalizedIdx
    end
    
    methods(Static)
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            p = inputParser();
            p.addParamValue('nBasesProjKeep', NaN, @isscalar);
            p.addParamValue('nBasesProjPerMarginalization', NaN, @(x) isscalar(x) || isvector(x));
            p.addParamValue('lambda', [], @(x) isempty(x) || isscalar(x));
            p.addParamValue('axesCombineSpecificMarginalizations', {}, @iscell);
            p.addParamValue('axesCombineAllMarginalizations', {}, @iscell);
            p.addParamValue('combineAxesWithTime', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            proj = ProjDPCA_NonOrthogonal();
            proj.nBasesProjPerMarginalization = p.Results.nBasesProjPerMarginalization;
            proj.nBasesProjKeep = p.Results.nBasesProjKeep;
            proj.lambda = p.Results.lambda;
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, p.Unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            p = inputParser();
            p.addParamValue('nBasesProjKeep', NaN, @isscalar);
            p.addParamValue('nBasesProjPerMarginalization', NaN, @(x) isscalar(x) || isvector(x));
            p.addParamValue('lambda', [], @(x) isempty(x) || isscalar(x));
            p.addParamValue('axesCombineSpecificMarginalizations', {}, @iscell);
            p.addParamValue('axesCombineAllMarginalizations', {}, @iscell);
            p.addParamValue('combineAxesWithTime', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            proj = ProjDPCA_NonOrthogonal();
            proj.nBasesProjPerMarginalization = p.Results.nBasesProjPerMarginalization;
            proj.nBasesProjKeep = p.Results.nBasesProjKeep;
            proj.lambda = p.Results.lambda;
            proj.axesCombineSpecificMarginalizations = p.Results.axesCombineSpecificMarginalizations;
            proj.axesCombineAllMarginalizations = p.Results.axesCombineAllMarginalizations;
            proj.combineAxesWithTime = p.Results.combineAxesWithTime;
            
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, p.Unmatched);
        end
    end

    methods
        function proj = ProjDPCA_NonOrthogonal(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            pset = pset.meanSubtractBases();
        end
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
            
            % filter for non-singular axes
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1;
            
            % build the list of covariates and covariate interactions to marginalize along
            [proj.combinedParams, proj.marginalizationNames] = StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization( ...
                pset.conditionDescriptor.axisAttributes, ...
                'axisIncludeMask', dimMask, ...
                'axisNames', pset.conditionDescriptor.axisNames, ...
                'combineAxesWithTime', proj.combineAxesWithTime, ...
                'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
                'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations);
                   
            M = numel(proj.marginalizationNames);
            nBasesProjKeep = proj.nBasesProjKeep; %#ok<*PROP>
            nBasesProjPerMarginalization = proj.nBasesProjPerMarginalization;
           
            if isnan(nBasesProjPerMarginalization)
                if isnan(nBasesProjKeep)
                    nBasesProjPerMarginalization = repmat(3, M, 1);
                    nBasesProjKeep = sum(nBasesProjPerMarginalization);
                else
                    % specified nBasesProjKeep only, generate the
                    % same number of bases per marginalization to give them
                    % each a shot to capture all of the "kept" set of
                    % nBasesProjKeep bases at the end.
                    nBasesProjPerMarginalization = repmat(nBasesProjKeep, M, 1);
                end
            else
                if numel(nBasesProjPerMarginalization) ~= 1 && numel(nBasesProjPerMarginalization) ~= M
                    error('nBasesProjPerMarginalization must be scalar or a vector with length nMarg = %d (%s)', M, strjoin(proj.marginalizationNames, ', '));
                end
                if isscalar(nBasesProjPerMarginalization)
                    nBasesProjPerMarginalization = repmat(nBasesProjPerMarginalization, M, 1);
                end
                if isnan(nBasesProjKeep)
                     % keep all of the bases generated
                    nBasesProjKeep = sum(nBasesProjPerMarginalization);
                end
            end
            
            proj.nBasesProjKeep = nBasesProjKeep;
            proj.nBasesProjPerMarginalization = nBasesProjPerMarginalization;
            
            % Extract trial averaged data
            debug('Extracting trial-averaged data for DPCA\n');
            NvbyTAbyAttr = pset.buildNbyTAbyConditionAttributes('validBasesOnly', true);
            Nv = size(NvbyTAbyAttr, 1);
            
            % compute optimal lambda
            if isempty(proj.lambda) || isnan(proj.lambda)
                if ~pset.hasDataByTrial
                    debug('Data by trial needed to compute optimal lambda, using default lambda %g\n', proj.defaultLambda);
                    proj.lambda = proj.defaultLambda;
                else
                    [NvbyTAbyAttrbyTrials] = pset.buildNbyTAbyConditionAttributesbyTrials(...
                        'maxTrials', proj.nIterationsOptimizeLambda, ...
                        'minimizeMissingSamples', true, ...
                        'validBasesOnly', true, ...
                        'message', 'Building individual trials tensor for optimizing lambda');
                    
                    nTrials_NvbyC = TensorUtils.squeezeDims(max(pset.dataNTrials(:, pset.basisValid, :), [], 1), 1);
                    nTrials_NvbyAttr = reshape(nTrials_NvbyC, [Nv pset.conditionsSizeNoExpand]);
                    
                    % remove timepoints with missing data for DPCA
                    nanMaskT_avg = TensorUtils.anyMultiDim(isnan(NvbyTAbyAttr), TensorUtils.otherDims(NvbyTAbyAttr, 2));
                    nanMaskT_single = TensorUtils.anyMultiDim(isnan(NvbyTAbyAttrbyTrials), TensorUtils.otherDims(NvbyTAbyAttrbyTrials, 2));
                    nanMaskT = squeeze(nanMaskT_avg | nanMaskT_single);
                    
                    NvbyTAbyAttr = TensorUtils.selectAlongDimension(NvbyTAbyAttr, 2, ~nanMaskT);
                    NvbyTAbyAttrbyTrials = TensorUtils.selectAlongDimension(NvbyTAbyAttrbyTrials, 2, ~nanMaskT);
                    
                    debug('Computing optimal DPCA regularization lambda\n');
                    proj.lambda = TrialDataUtilities.DPCA.dpca_optimizeLambda(...
                        NvbyTAbyAttr, NvbyTAbyAttrbyTrials, nTrials_NvbyAttr, ...
                        'combinedParams', proj.combinedParams, ...
                        'nBasesKeep', nBasesProjKeep, ...
                        'nBasesPerMarginalization', nBasesProjPerMarginalization, ...
                        'numRep', proj.nIterationsOptimizeLambda);
                end
            else
                debug('Using set lambda %g\n', proj.lambda);
            end
            
            % run dpca
            [decoderNvbyK, encoderNvbyK, proj.marginalizationIdxByBasis, proj.dataMarginalized, proj.dataMarginalizedIdx] = TrialDataUtilities.DPCA.dpca(NvbyTAbyAttr, ...
                'nBasesKeep', nBasesProjKeep, ...
                'nBasesPerMarginalization', nBasesProjPerMarginalization, ...
                'combinedParams', proj.combinedParams, ...
                'lambda', proj.lambda);
            
            % inflate from only valid bases to full nBases by filling with NaNs
            decoderKbyN = TensorUtils.inflateMaskedTensor(decoderNvbyK, 1, pset.basisValid)';
            encoderNbyK = TensorUtils.inflateMaskedTensor(encoderNvbyK, 1, pset.basisValid);
        end

        function names = getBasisNames(proj, pset) %#ok<INUSD>
            names = arrayfun(@(i) sprintf('DPC %d', i), ...
                    (1:proj.nBasesProj)', 'UniformOutput', false);
        end
    end

end

