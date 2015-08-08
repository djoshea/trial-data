classdef ProjDPCA_NonOrthogonal < StateSpaceProjection

    properties
        lambda = [];
        optimizeLambda = true;
        maxTrialsForOptimizeLambda = 50;
        
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
    
    properties
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
    end
    
    methods(Static)
        function [proj, stats, psetPrepared] = createFrom(pset, varargin)
            proj = ProjDPCA_NonOrthogonal();
            [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            proj = ProjDPCA_NonOrthogonal();
            [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, varargin{:});
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
        function [decoderKbyN, encoderNbyK] = computeProjectionCoefficients(proj, pset, varargin)
            p = inputParser;
            p.addParamValue('nBasesProj', min(20, pset.nBasesValid), @isscalar);
            p.parse(varargin{:});
            K = p.Results.nBasesProj;
            
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
                   
            % Extract trial averaged data
            debug('Extracting trial-averaged data for DPCA\n');
            NvbyCbyTA = pset.buildNbyCbyTA('validBasesOnly', true);
            Nv = size(NvbyCbyTA, 1);
            TA = size(NvbyCbyTA, 3);
            condSize = makerow(pset.conditionDescriptor.conditionsSize);
            NvbyAttrbyTA = reshape(NvbyCbyTA, [Nv condSize TA]);
            nAxes = numel(condSize);
            
            % compute optimal lambda
            if proj.optimizeLambda && (isempty(proj.lambda) || isnan(proj.lambda))
                if ~pset.hasDataByTrial
                    debug('Data by trial needed to compute optimal lambda, using default lambda %g\n', proj.defaultLambda);
                    proj.lambda = proj.defaultLambda;
                else
                    [NvbyTAbyCbyTrials] = pset.buildNbyTAbyCbyTrials(...
                        'maxTrials', proj.nIterationsOptimizeLambda, ...
                        'minimizeMissingSamples', true, ...
                        'validBasesOnly', true, ...
                        'message', 'Building individual trials tensor for optimizing lambda');
                    
                    maxTrials = size(NvbyTAbyCbyTrials, 4);
                    newSize = [Nv condSize TA maxTrials];
                    NvbyAttrbyTAbyTrials = reshape(permute(NvbyTAbyCbyTrials, [1 3 2 4]), newSize);
                    
                    nTrials_NvbyC = TensorUtils.squeezeDims(max(pset.dataNTrials(:, pset.basisValid, :), [], 1), 1);
                    nTrials_NvbyAttr = reshape(nTrials_NvbyC, [Nv condSize]);
                    
                    % remove timepoints with missing data for DPCA
                    nanMaskT_avg = TensorUtils.anyMultiDim(isnan(NvbyAttrbyTA), 1:(1+nAxes));
                    nanMaskT_single = TensorUtils.anyMultiDim(isnan(NvbyAttrbyTAbyTrials), [1:(1+nAxes), nAxes+3]);
                    nanMaskT = squeeze(nanMaskT_avg | nanMaskT_single);
                    
                    NvbyAttrbyTA = TensorUtils.selectAlongDimension(NvbyAttrbyTA, 2+nAxes, ~nanMaskT);
                    NvbyAttrbyTAbyTrials = TensorUtils.selectAlongDimension(NvbyAttrbyTAbyTrials, 2+nAxes, ~nanMaskT);
                    
                    debug('Computing optimal DPCA regularization lambda\n');
                    proj.lambda = TrialDataUtilities.DPCA.dpca_optimizeLambda(...
                        NvbyAttrbyTA, NvbyAttrbyTAbyTrials, nTrials_NvbyAttr, ...
                        'combinedParams', proj.combinedParams, ...
                        'numRep', proj.nIterationsOptimizeLambda);
                end
            else
                debug('Using set lambda %g\n', proj.lambda);
            end
            
            % run dpca
            [decoderNvbyK, encoderNvbyK, whichMarg] = TrialDataUtilities.DPCA.dpca(NvbyAttrbyTA, K, ...
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

