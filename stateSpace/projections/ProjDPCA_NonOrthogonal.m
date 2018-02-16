classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        lambda = [];
        maxTrialsForOptimizeLambda = 50;
        
        % specify this to select the cross-validation optimized lambda for
        % a specific marginalizationSpec
        useOptimizedLambdaForMarginalization = {};
        
        % values found through cross-validation if no labmda is specified
        optimizedLambdaTotal
        optimizedLambdaPerMarginalization
        optimizedLambdaStats % struct describing the actual cross-validation errors
        
        nBasesProjKeep = NaN;
        nBasesProjPerMarginalization = NaN; % the number of bases per marginalization to find, before selecting the top nBasesProj. this can be a vector, one for each marginalization type too
        
        nIterationsOptimizeLambda = 5; % higher is more accurate, 5-10 is reasonable
        
        orderBasesByVariance = true;
    end
    
    properties(Constant)
        defaultLambda = 0;
    end
    
    properties(SetAccess=protected)
        marginalizationIdxByBasis
        
        dataMarginalized
        dataMarginalizedIdx
    end
    
    methods(Static)
        function [proj, unmatched] = parseCreateParams(proj, varargin)
            p = inputParser;
            p.addParameter('nBasesProjKeep', NaN, @isscalar);
            p.addParameter('nBasesProjPerMarginalization', NaN, @(x) isscalar(x) || isvector(x));
            p.addParameter('lambda', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('useOptimizedLambdaForMarginalization', {}, @(x) isscalar(x) || ischar(x) || iscell(x));
            p.addParameter('axesCombineSpecificMarginalizations', {}, @iscell);
            p.addParameter('axesCombineAllMarginalizations', {}, @iscell);
            p.addParameter('combineAxesWithTime', true, @islogical);
            p.addParameter('orderBasesByVariance', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unmatched = p.Unmatched;
            
            proj.nBasesProjPerMarginalization = p.Results.nBasesProjPerMarginalization;
            proj.nBasesProjKeep = p.Results.nBasesProjKeep;
            proj.lambda = p.Results.lambda;
            proj.useOptimizedLambdaForMarginalization = p.Results.useOptimizedLambdaForMarginalization;
            proj.axesCombineSpecificMarginalizations = p.Results.axesCombineSpecificMarginalizations;
            proj.axesCombineAllMarginalizations = p.Results.axesCombineAllMarginalizations;
            proj.combineAxesWithTime = p.Results.combineAxesWithTime;
            proj.orderBasesByVariance = p.Results.orderBasesByVariance;
        end
        
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjDPCA_NonOrthogonal();
            [proj, unmatched] = ProjDPCA_NonOrthogonal.parseCreateParams(proj, varargin{:});
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjDPCA_NonOrthogonal();
            [proj, unmatched] = ProjDPCA_NonOrthogonal.parseCreateParams(proj, varargin{:});
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, unmatched);
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
            M = numel(proj.marginalizationNames);
            nBasesProjKeep = proj.nBasesProjKeep; %#ok<*PROPLC,*PROP>
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
                if isscalar(nBasesProjPerMarginalization)
                    nBasesProjPerMarginalization = repmat(nBasesProjPerMarginalization, M, 1);
                end
                if numel(nBasesProjPerMarginalization) < M
                    warning('Expecting %d marginalizations in nBasesProjPerMarginalization, expanding with zeros', M);
                    nBasesProjPerMarginalization = [makecol(nBasesProjPerMarginalization); zeros(M - numel(nBasesProjPerMarginalization), 1)];
                end
                
                if numel(nBasesProjPerMarginalization) ~= 1 && numel(nBasesProjPerMarginalization) ~= M
                    error('nBasesProjPerMarginalization must be scalar or a vector with length nMarg = %d (%s)', M, strjoin(proj.marginalizationNames, ', '));
                end
                
                if isnan(nBasesProjKeep)
                     % keep all of the bases generated
                    nBasesProjKeep = sum(nBasesProjPerMarginalization);
                end
            end
            
            if any(nBasesProjPerMarginalization > pset.nBasesValid)
                error('nBasesProjPerMarginalization should not exceed nBasesValid');
            end
            
            proj.nBasesProjKeep = nBasesProjKeep;
            proj.nBasesProjPerMarginalization = nBasesProjPerMarginalization;
            
            % Extract trial averaged data
            debug('Extracting trial-averaged data for DPCA\n');
            NvbyTAbyAttr = pset.arrangeNbyTAbyConditionAttributes('validBasesOnly', true);
            Nv = size(NvbyTAbyAttr, 1);
            
            % compute optimal lambda
            if isempty(proj.lambda) || isnan(proj.lambda)
                if ~pset.hasDataByTrial && ~pset.hasCachedDataByTrial
                    debug('Data by trial needed to compute optimal lambda, using default lambda %g\n', proj.defaultLambda);
                    proj.lambda = proj.defaultLambda;
                else
                    % lookup the desired lambda first, in case there's an
                    % issue
                    if ~isempty(proj.useOptimizedLambdaForMarginalization)
                        idxMargLambda = StateSpaceProjectionStatistics.staticFindFlatMarginalizationSpecInList(...
                            proj.useOptimizedLambdaForMarginalization, proj.marginalizationList, 'ignoreMissingTime', true);
                        lambdaMarg = proj.useOptimizedLambdaForMarginalization;
                        if ischar(lambdaMarg)
                            lambdaMarg = {lambdaMarg};
                        end
                        if isempty(idxMargLambda)
                            error('Could not find marginalization %s', strjoin(lambdaMarg, ' x '));
                        end
                        assert(isscalar(idxMargLambda));
                    end  
                    
                    [meansExcludingTrials_A_NbyCbyTbyR, sampledTrials_A_NbyCbyTbyR, nTrials_NbyC_sampled] = ...
                        pset.computeDataMeanExcludingSampledTrials('numTrialsDefault', proj.nIterationsOptimizeLambda, 'validBasesOnly', true);
                    
                    meansExcluding_NbyTAbyCbyR = permute(TensorUtils.cat(3, meansExcludingTrials_A_NbyCbyTbyR{:}), [1 3 2 4]);
                    trials_NbyTAbyCbyR = permute(TensorUtils.cat(3, sampledTrials_A_NbyCbyTbyR{:}), [1 3 2 4]);
                    
%                     [meansExcluding_NbyTAbyCbyR, trials_NbyTAbyCbyR, nTrials_NbyC_sampled] = ...
%                         pset.computeDataMeansExcludingSampledTrials(...
%                         'maxTrials', proj.nIterationsOptimizeLambda, ...
%                         'validBasesOnly', true, ...
%                         'message', 'Building individual trials tensor for optimizing lambda');
                    
                    if any(nTrials_NbyC_sampled < proj.nIterationsOptimizeLambda)
                        error('Not enough trials to optimize lambda');
                    end
                    
                    % remove timepoints with missing data for DPCA
                    nanMaskT_avg = TensorUtils.anyMultiDim(isnan(NvbyTAbyAttr), TensorUtils.otherDims(NvbyTAbyAttr, 2));
                    nanMaskT_single = TensorUtils.anyMultiDim(isnan(trials_NbyTAbyCbyR), TensorUtils.otherDims(trials_NbyTAbyCbyR, 2));
                    nanMaskT = squeeze(nanMaskT_avg | nanMaskT_single);
                    
                    NvbyTAbyAttr = TensorUtils.selectAlongDimension(NvbyTAbyAttr, 2, ~nanMaskT);
                    trials_NbyTAbyCbyR = TensorUtils.selectAlongDimension(trials_NbyTAbyCbyR, 2, ~nanMaskT);
                    meansExcluding_NbyTAbyCbyR = TensorUtils.selectAlongDimension(meansExcluding_NbyTAbyCbyR, 2, ~nanMaskT);
                    
                    sz = size(trials_NbyTAbyCbyR);
                    trials_NbyTAbyAttrbyR = reshape(trials_NbyTAbyCbyR, [sz(1), sz(2), makerow(pset.conditionsSize), sz(4)]); 
                    meansExcluding_NbyTAbyAttrbyR = reshape(meansExcluding_NbyTAbyCbyR, [sz(1), sz(2), makerow(pset.conditionsSize), sz(4)]); 
                    debug('Computing optimal DPCA regularization lambda\n');
                    [proj.optimizedLambdaTotal, proj.optimizedLambdaPerMarginalization, proj.optimizedLambdaStats] = TrialDataUtilities.DPCA.dpca_optimizeLambda(...
                        NvbyTAbyAttr, trials_NbyTAbyAttrbyR, meansExcluding_NbyTAbyAttrbyR, ...
                        'combinedParams', proj.combinedParams, ...
                        'nBasesKeep', nBasesProjKeep, ...
                        'nBasesPerMarginalization', nBasesProjPerMarginalization, ...
                        'numRep', proj.nIterationsOptimizeLambda);
                    
                    if ~isempty(proj.useOptimizedLambdaForMarginalization)
                        % looked up which lambda above
                        debug('Using lambda optimized for marginalization %s\n', strjoin(lambdaMarg, ' x '));
                        proj.lambda = proj.optimizedLambdaPerMarginalization(idxMargLambda);
                    else
                        debug('Using lambda optimized for all marginalizations\n');
                        proj.lambda = proj.optimizedLambdaTotal;
                    end               
                end
            else
                debug('Using set lambda %g\n', proj.lambda);
            end
            
            % run dpca
            [decoderNvbyK, encoderNvbyK, proj.marginalizationIdxByBasis, proj.dataMarginalized, proj.dataMarginalizedIdx] = TrialDataUtilities.DPCA.dpca(NvbyTAbyAttr, ...
                'nBasesKeep', nBasesProjKeep, ...
                'nBasesPerMarginalization', nBasesProjPerMarginalization, ...
                'combinedParams', proj.combinedParams, ...
                'lambda', proj.lambda, ...
                'order', proj.orderBasesByVariance);
            
            % inflate from only valid bases to full nBases by filling with NaNs
            decoderKbyN = TensorUtils.inflateMaskedTensor(decoderNvbyK, 1, pset.basisValid)';
            encoderNbyK = TensorUtils.inflateMaskedTensor(encoderNvbyK, 1, pset.basisValid);
        end
        
        function [nameStem, names] = generateBasisNameProjStem(proj, pset) %#ok<INUSD>
            names = cellvec(proj.nBasesProj);
            for i = 1:proj.nBasesProj
                names{i} = sprintf('%s %d', proj.marginalizationNames{proj.marginalizationIdxByBasis(i)}, i);
            end
            nameStem = 'DPC';
        end
    end

end

