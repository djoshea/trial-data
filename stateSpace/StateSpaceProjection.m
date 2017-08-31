classdef StateSpaceProjection 
% Abstract base class representing a set of coefficients per basis with which
% to project PopulationTrajectorySet instances. Typically these are built using 
% fromPopulationTrajectorySet to compute the coefficients, in a manner determined
% by subclasses

    properties(SetAccess=protected)
        initialized = false;
        translationNormalization % stored translation / normalization stored when building and used BEFORE projecting
        
        % used AFTER projecting, to transform the output bases, mainly used
        % when inverting a projection
        translationNormalizationPostProject
        
        % given an nBasesSource x T matrix of data in the original basis, we project
        % into the new basis using decoderKbyN * neural_NbyT which yields
        % decoded_KbyT. We can then reconstruct the data in the original
        % basis using reconstruction_NbyT = encoderNbyK * decodedKbyT, or
        % reconstruction_NbyT = encoderNbyK * decoderKbyN * neural_NbyT 
        %
        % These matrices will contain zeros where either 
        % ~basisValidProj (K dimension) or
        % ~basisValid (N dimension)
        decoderKbyN
        encoderNbyK
        
        buildStats % StateSpaceProjectionStatistics instance created at build time
        
        basisValid % nBasesSource x 1 logical vector indicating which bases were considered valid in the projection (typically copied from the pset from which I am built)
        basisInvalidCause % nBasesSource x 1 cell str (typically copied from the pset from which I am built)
        
        % nBasesProj x 1 logical vector indicating which OUTPUT bases will be considered valid
        % this is mostly useful when building inverses, so that you get the
        % original basis valid mask back when inverting a projection
        basisValidProj 
        basisInvalidCauseProj
        
        basisNamesSource % nBasesSource x 1 cellstr of basis names from the pset I was built off
        dataUnitsSource = '';
        
        % stores meta data about the types of marginalizations performed
        % when buildingn statistics during construction
        
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
        
        axisIncludeList
        
        marginalizationList
    end
    
    properties
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
    
    properties(Dependent)
        % nBasesProj x 1 cell of basis names, optional, that I will pick for the basis names of the projected 
        % output. If set, this will override the method "getBasisNames"
        % which typically picks these names based on the input pset.
        % this is primarily used for inverse projections, when we want to
        % restore the basis names of the pset this projection was built off 
        basisNamesProj
    end
    
    properties(SetAccess=protected)
        basisNamesProjManual % stores manual override for basisNamesProj
    end
    
    properties
        basisNamesProjStem % string which defines the stem to use for each basis, ' %d' will be added to each basis if this is set and basisNamesProjManual isempty
    end
    
    properties(Dependent)
        nBasesSource
        nBasesProj
        
        decoderKbyNValid
        encoderNbyKValid
        
        decoderKbyNZeroInvalid
        encoderNbyKZeroInvalid
    end
    
    methods(Abstract, Static)
        % this abstract static method is basically a way of getting around matlab's lack
        % of support for the factory builder pattern. Define this method
        % like this so that buildFromPopulationTrajectorySet() below can do
        % the heavy lifting.
        % methods(Static)
        %     function [proj, stats, psetPrepared] = createFrom(pset, varargin)
        %         proj = YOURSUBCLASSNAME() you can send in varargin{:} if
        %           there are any options you wish to receive
        %         [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, varargin{:});
        %     end
        
        %     function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
        %         proj = YOURSUBCLASSNAME() you can send in varargin{:} if there are options > );
        %         [proj, psetProjected, stats] = proj.buildFromAndProjectPopulationTrajectorySet(pset, varargin{:})
        %     end
        % end
        %
        [proj, stats, psetPrepared] = createFrom(pset, varargin);
        
        [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
    end

    methods(Abstract)
        % compute the N * K matrix of basis coefficients for the projection
        % given an N x T matrix of data in the original basis, we project
        % into the new basis using decoderKbyN * neural_NbyT which yields
        % decoded_KbyT. We can then reconstruct the data in the original
        % basis using reconstruction_NbyT = encoderNbyK * decodedKbyT, or
        % reconstruction_NbyT = encoderNbyK * decoderKbyN * neural_NbyT 
        % 
        % optional params:
        % 'nBasesProj': number of bases to project into, if requested by
        %     user
        [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
    end  
    
    % Methods most classes may wish to override
    methods
        % return either a char nBasesProj x 1 cell str: list of basis names for the new basis  
        function [nameStem, names] = generateBasisNameProjStem(proj, pset) %#ok<INUSD>
            nameStem = 'Basis';
            names = {};
        end
        
        % this function can optionally use pset to create custom names
        % based on pset.basisNames
        function names = generateBasisNamesProj(proj, pset) %#ok<INUSD>
            if ~isempty(proj.basisNamesProjManual) && numel(proj.basisNamesProjManual) == proj.nBasesProj
                names = proj.basisNamesProjManual;
            else
                if ~isempty(proj.basisNamesProjStem)
                    stem = proj.basisNamesProjStem;
                else
                    stem = 'Basis';
                end
                
                names = makecol(arrayfun(@(idx) sprintf('%s %d', stem, idx), 1:proj.nBasesProj, 'UniformOutput', false));
            end
        end
        
        function names = get.basisNamesProj(proj)
            names = makecol(proj.generateBasisNamesProj());
        end
        
        function proj = set.basisNamesProj(proj, names)
            assert(numel(names) == proj.nBasesProj);
            proj.basisNamesProjManual = names;
%             if ~isempty(names)
%                 % only have one or the other set
%                 proj.basisNamesProjStem = '';
%             end
        end
        
        function stem = get.basisNamesProjStem(proj)
            if isempty(proj.basisNamesProjStem)
                stem = '';
            else
                stem = proj.basisNamesProjStem;
            end
        end
                
        function proj = set.basisNamesProjStem(proj, stem)
            proj.basisNamesProjStem = stem;
%             if ~isempty(stem)
%                 proj.basisNamesProjManual = {}; %#ok<MCSUP>
%             end
        end
            
        function proj = markBasesInvalid(proj, mask, cause)
            proj.warnIfNoArgOut(nargout);
  
            mask = makecol(TensorUtils.vectorIndicesToMask(mask, proj.nBasesSource));
            maskOrig = mask; % cache for later
            
            mask = mask & proj.basisValid;
            proj.basisValid(mask) = false;
            
            if ~any(mask)
                return;
            end
            
            if nargin < 3
                cause = 'markBasesInvalid cause unspecified';
            end
            assert(iscellstr(cause) || ischar(cause), 'Cause must be a cellstr or a string');
            if ischar(cause)
                cause = repmat({cause}, nnz(mask), 1);
            end
            cause = makecol(cause);
          
            if numel(cause) == nnz(maskOrig)
                cause = TensorUtils.inflateMaskedTensor(cause, 1, maskOrig, {''});
            end
            assert(numel(cause) == numel(mask), 'Length of cellstr cause must match nnz(mask) or numel(mask)');
           
            proj.basisInvalidCause(mask) = cause(mask);
        end
    end
    
    % Simple dependent property getters
    methods
        function n = get.nBasesSource(proj)
            if isempty(proj.decoderKbyN)
                n = NaN;
            else
                n = size(proj.decoderKbyN, 2);
            end
        end

        function n = get.nBasesProj(proj)
            if isempty(proj.decoderKbyN)
                n = NaN;
            else
                n = size(proj.decoderKbyN, 1);
            end
        end
        
        % Don't remove these, other code may rely on
        % these matrices being zeroed
        function d = get.decoderKbyNZeroInvalid(proj)
            d = proj.decoderKbyN;
            d(~proj.basisValidProj, ~proj.basisValid) = 0;
        end
        
        function e = get.encoderNbyKZeroInvalid(proj)
            e = proj.encoderNbyK;
            e(~proj.basisValid, ~proj.basisValidProj) = 0;
        end
        
        function d = get.decoderKbyNValid(proj)
            d = proj.decoderKbyN(proj.basisValidProj, proj.basisValid);
        end
        
        function e = get.encoderNbyKValid(proj)
            e = proj.encoderNbyK(proj.basisValid, proj.basisValidProj);
        end
        
        function v = get.basisValid(proj)
            if isempty(proj.basisValid)
                if ~isnan(proj.nBasesSource)
                    v = truevec(proj.nBasesSource);
                else
                    v = [];
                end 
            else
                v = proj.basisValid;
            end
        end
        
        function v = get.basisValidProj(proj)
            if isempty(proj.basisValidProj)
                if ~isnan(proj.nBasesProj)
                    v = truevec(proj.nBasesProj);
                else
                    v = [];
                end
            else
                v = proj.basisValidProj;
            end
        end
        
        function v = get.basisInvalidCauseProj(proj)
            if isempty(proj.basisInvalidCauseProj)
                if ~isnan(proj.nBasesProj)
                    v = cellstrvec(proj.nBasesProj);
                else
                    v = [];
                end
            else
                v = proj.basisInvalidCauseProj;
            end
        end
        
        function v = get.basisInvalidCause(proj)
            if isempty(proj.basisInvalidCause)
                if ~isnan(proj.nBasesSource)
                    v = cellstrvec(proj.nBasesSource);
                else
                    v = [];
                end
            else
                v = proj.basisInvalidCause;
            end
        end
        
        function norms = calculateDecoderNorms(proj)
            norms = rownorms(proj.decoderKbyN);
        end
        
        function tf = isDecoderOrthogonal(proj)
            dd = proj.decoderKbyN * proj.decoderKbyN';
            dd = dd - diag(diag(dd));
            tf = max(abs(dd(:))) < 1e-09;
        end
        
        function tf = isDecoderOrthonormal(proj)
            dd = proj.decoderKbyN * proj.decoderKbyN';
            dd = dd - eye(size(dd));
            tf = max(abs(dd(:))) < 1e-09;
        end
    end

    methods
        function [marginalizationNames, combinedParams, axisIncludeList, marginalizationList] = getMarginalizationNames(proj, pset, varargin)
            % other args include
            % 'combineAxesWithTime'
            % 'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
            % 'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations);
            
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1; % filter for non-singular axes
            
            [combinedParams, marginalizationNames, axisIncludeList, marginalizationList] = StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization( ...
                pset.conditionDescriptor.axisAttributes, ...
                'axisIncludeMask', dimMask, ...
                'axisNames', pset.conditionDescriptor.axisNames, ...
                'combineAxesWithTime', proj.combineAxesWithTime, ...
                'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
                'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations);
        end
        
        function [proj, stats, psetPrepared] = buildFromPopulationTrajectorySet(proj, pset, varargin)
            % build this projection matrix based on an existing PopulationTrajectorySet
            % defers to calculateProjectionMatrix for the actual basis computation
            proj.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('pset', @(x) isa(x, 'PopulationTrajectorySet'));
            p.addParameter('computeStatistics', true, @islogical);
            p.addParameter('meanSubtractForStatistics', true, @islogical); % mean subtract data when computing statistics, this makes sense to turn off if the data is already measured relative to some meaningful baseline
            p.addParameter('nBasesProj', NaN, @isscalar); % this is optional for both the caller to provide and the projection to obey
%             p.addParameter('computeStatisticsMarginalized', true, @islogical);
            p.addParameter('computeStatisticsForRandomized', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(pset, varargin{:});
            
            % we do this here simply to force calculation of all compute on
            % demand properties, so that they get computed before we
            % initiate a copy below
            PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            pset.dataDifferenceOfTrialsScaledNoiseEstimate;
           
            % build the list of covariates and covariate interactions to marginalize along
            % in case the projection needs it for its own purposes
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1; % filter for non-singular axes
            [proj.combinedParams, proj.marginalizationNames, proj.axisIncludeList, proj.marginalizationList] = StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization( ...
                pset.conditionDescriptor.axisAttributes, ...
                'axisIncludeMask', dimMask, ...
                'axisNames', pset.conditionDescriptor.axisNames, ...
                'combineAxesWithTime', proj.combineAxesWithTime, ...
                'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
                'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations);
                        
            % make any necessary transformations, particularly translation / normalization
            pset = proj.preparePsetForInference(pset);
            psetPrepared = pset;

            % extract the translation normalization that will be applied before projection
            proj.translationNormalization = pset.translationNormalization;

            % compute the coefficients for the projection
            debug('Computing projection encoder and decoder coefficients\n');
            [decoderKbyN, encoderNbyK, proj] = proj.computeProjectionCoefficients(pset, 'nBasesProj', p.Results.nBasesProj); %#ok<*PROPLC>
            % don't put these in the line above b/c assignment to proj will
            % override it
            proj.decoderKbyN = decoderKbyN;
            proj.encoderNbyK = encoderNbyK;
            
            assert(size(proj.decoderKbyN, 2) == pset.nBases, 'Decoder matrix returned by computeProjectionCoefficients must match pset.nBases along dim 2');
            assert(size(proj.encoderNbyK, 1) == pset.nBases, 'Encoder matrix returned by computeProjectionCoefficients must match pset.nBases along dim 1');
            
            % copy the basis valid mask
            proj.basisValid = pset.basisValid;
            proj.basisInvalidCause = pset.basisInvalidCause;
            proj.basisNamesSource = pset.basisNames;
            proj.dataUnitsSource = pset.dataUnits;
            [proj.basisNamesProjStem, proj.basisNamesProjManual] = proj.generateBasisNameProjStem(pset);
            
            % set coeff to zero on invalid bases
            proj.decoderKbyN(:, ~proj.basisValid) = 0;
            proj.encoderNbyK(~proj.basisValid, :) = 0;
            
            % results will store statistics and useful quantities related to the
            % projection
            if p.Results.computeStatistics
                stats = StateSpaceProjectionStatistics.build(proj, pset, 'meanSubtract', ...
                    p.Results.meanSubtractForStatistics, ... % 'marginalize', p.Results.computeStatisticsMarginalized, ...
                    'computeForRandomized', p.Results.computeStatisticsForRandomized', ...
                    'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations, ...
                    'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
                    'combineAxesWithTime', proj.combineAxesWithTime, ...
                    p.Unmatched);
                proj.buildStats = stats;
            else
                stats = [];
                proj.buildStats = [];
            end

            proj.initialized = true;
        end

        function [psetProjected, stats] = projectPopulationTrajectorySet(proj, pset, varargin)
            p = inputParser();
            p.addParameter('clearBeforeApplyingTranslationNormalization', true, @islogical); % clear existing pset trnorm first
            p.addParameter('applyTranslationNormalization', true, @islogical); % apply proj.trNorm to pset before projecting
            p.addParameter('applyTranslationNormalizationPostProject', true, @islogical); % apply the post projection trNorm, mainly used for inverse projections
            p.addParameter('meanSubtractForStatistics', true, @islogical); % mean subtract data when computing statistics, this makes sense to turn off if the data is already measured relative to some meaningful baseline
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            assert(pset.nBases == proj.nBasesSource, ...
                'Number of bases must match in order to project');

            if any(pset.basisValid & ~proj.basisValid)
                % ensure proj-invalid bases are marked invalid
                % this isn't strictly necessary but just in case
                pset = pset.markBasesPermanentlyInvalid(~proj.basisValid, 'invalidated before state space projection');
            end
            
            if any(proj.basisValid & ~pset.basisValid)
                error('PopulationTrajectorySet has invalid bases not marked invalid in StateSpaceProjection. You should equalize the bases invalid to get consistent results');
            end
            
            % we do this here simply to force calculation of all compute on
            % demand properties, so that they get computed before we
            % initiate a copy below
            PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            pset.dataDifferenceOfTrialsScaledNoiseEstimate;
            
            % replace translation normalization
            if p.Results.applyTranslationNormalization && ~isempty(proj.translationNormalization)
                debug('Applying translation/normalization associated with projection to data\n');
                if p.Results.clearBeforeApplyingTranslationNormalization
                    pset = pset.clearTranslationNormalization();
                end
                pset = pset.applyTranslationNormalization(proj.translationNormalization);
            end
            
            % copy basic settings from pset 
            b = PopulationTrajectorySetBuilder.copySettingsDescriptorsFromPopulationTrajectorySet(pset);
%             b.dataUnits = pset.dataUnits; % user may may want to clear
%             this later, but keep them by default. No just clear them,
%             they are misleading

            b.basisNames = proj.generateBasisNamesProj(pset);
            
            b.basisUnits = proj.getBasisUnits(pset); 
            
            % copy/compute trial averaged data 
            b.tMinForDataMean = pset.tMinForDataMean;
            b.tMaxForDataMean = pset.tMaxForDataMean;

            % sum across bases for projection (dataNTrials is nAlign x nBases x nConditions)
            b.dataNTrials = repmat(sum(pset.dataNTrials, 2), [1, proj.nBasesProj, 1]);
            % all across bases for validity
            b.dataValid = repmat(all(pset.dataValid, 2), [1, proj.nBasesProj, 1]);
        
            % take max / min over all bases involved in each linear combination
            b.tMinValidByAlignBasisCondition = ...
                TensorUtils.linearCombinationApplyScalarFnAlongDimension(...
                pset.tMinValidByAlignBasisCondition, 2, proj.decoderKbyN, @max);
            b.tMaxValidByAlignBasisCondition = ...
                TensorUtils.linearCombinationApplyScalarFnAlongDimension(...
                pset.tMaxValidByAlignBasisCondition, 2, proj.decoderKbyN, @min);
            
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

                projMat = proj.decoderKbyN * mat;
                b.dataMean{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));

                % use sqrt(sd1^2 / n1 + sd2^2 / n2 + ...) formula
                % which here means semNew = sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                mat = reshape(pset.dataSem{iAlign}, pset.nBases, ...
                    pset.nConditions * pset.nTimeDataMean(iAlign));
                mat(~proj.basisValid, :) = 0;
                projMat = sqrt(abs(proj.decoderKbyN) * (mat.^2));
                b.dataSem{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));
            end
            
            % project single trial data if any
            if pset.simultaneous && ~isempty(pset.dataByTrial)
                [b.dataByTrial, b.tMinByTrial, b.tMaxByTrial] = deal(cell(proj.nBasesProj, pset.nAlign));
                
                % dataByTrial is nBases x nAlign cell of nTrials (R) x T_a
                for iAlign = 1:pset.nAlign
                    % will be R x T x N
                    tensor = cat(3, pset.dataByTrial{:, iAlign});
                    tensor(:, :, ~proj.basisValid) = 0;
                    tensor = TensorUtils.linearCombinationAlongDimension(...
                        tensor, 3, proj.decoderKbyN, 'replaceNaNWithZero', true);
                    
                    b.dataByTrial(:, iAlign) = TensorUtils.selectEachAlongDimension(tensor, 3);
                    
                    b.tMinByTrial(:, iAlign) = TensorUtils.selectEachAlongDimension(...
                        TensorUtils.linearCombinationApplyScalarFnAlongDimension(...
                            cat(2, pset.tMinByTrial{:, iAlign}), 2, proj.decoderKbyN, @nanmax), 2);
                    b.tMaxByTrial(:, iAlign) = TensorUtils.selectEachAlongDimension(...
                        TensorUtils.linearCombinationApplyScalarFnAlongDimension(...
                            cat(2, pset.tMaxByTrial{:, iAlign}), 2, proj.decoderKbyN, @nanmin), 2);
                end
                
                b.tMinForDataByTrial = pset.tMinForDataByTrial;
                b.tMaxForDataByTrial = pset.tMaxForDataByTrial;
                b.alignValidByTrial = pset.alignValidByTrial;
                
                % preserve the trial lists
                % nBases x nConditions --> nBasesProj x nConditions
                b.trialLists = repmat(pset.trialLists(1, :), proj.nBasesProj, 1);
            end
            
            % project cached single trial data if any
            if ~isempty(pset.dataCachedSampledTrialsTensor)
                tensor = pset.dataCachedSampledTrialsTensor;
                tensor(~proj.basisValid, :) = 0;
                b.dataCachedSampledTrialsTensor = TensorUtils.linearCombinationAlongDimension(...
                    tensor, 1, proj.decoderKbyN, 'replaceNaNWithZero', true);
                
                tensor = pset.dataCachedMeanExcludingSampledTrialsTensor;
                tensor(~proj.basisValid, :) = 0;
                b.dataCachedMeanExcludingSampledTrialsTensor = TensorUtils.linearCombinationAlongDimension(...
                    tensor, 1, proj.decoderKbyN, 'replaceNaNWithZero', true);
                
                % compute min over all trial counts included in each basis
                b.dataCachedSampledTrialCounts = TensorUtils.linearCombinationApplyScalarFnAlongDimension(...
                    pset.dataCachedSampledTrialCounts, 1, proj.decoderKbyN, @min);
            end
            
            % update difference of trials scaled noise estimates so that we
            % can compute noise variance floors when projecting. since the
            % noise estimates are already scaled by 1/sqrt(2*nTrials), we
            % simply add them together to get the new scaled estimate
            if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimate)
                b.dataDifferenceOfTrialsScaledNoiseEstimate = TensorUtils.linearCombinationAlongDimension(...
                    pset.dataDifferenceOfTrialsScaledNoiseEstimate, 1, abs(proj.decoderKbyN), ...
                    'replaceNaNWithZero', true, ...
                    'keepNaNIfAllNaNs', true, ...
                    'normalizeCoefficientsByNumNonNaN', false);
            end
          
            % project randomized data, recompute intervals
            if ~isempty(pset.dataMeanRandomized)
                [b.dataMeanRandomized] = deal(cell(pset.nAlign, 1));
                for iAlign = 1:pset.nAlign
                    % mat is N x CTS, coeff is N x K, where S is nRandomSamples 
                    mat = reshape(pset.dataMeanRandomized{iAlign}, pset.nBases, ...
                        pset.nConditions*pset.nTimeDataMean(iAlign)*pset.nRandomSamples);
                    mat(~proj.basisValid, :) = 0;
                    projMat = proj.decoderKbyN * mat;
                    b.dataMeanRandomized{iAlign} = reshape(projMat, proj.nBasesProj, ...
                        pset.nConditions, pset.nTimeDataMean(iAlign), pset.nRandomSamples);
                    
                    % use sqrt(sd1^2 / n1 + sd2^2 / n2 + ...) formula
                    % which here means semNew = sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                    mat = reshape(pset.dataSemRandomized{iAlign}, pset.nBases, ...
                        pset.nConditions * pset.nTimeDataMean(iAlign)*pset.nRandomSamples);
                    mat(~proj.basisValid, :) = 0;
                    projMat = sqrt(abs(proj.decoderKbyN) * (mat.^2));
                    b.dataSemRandomized{iAlign} = reshape(projMat, proj.nBasesProj, ...
                        pset.nConditions, pset.nTimeDataMean(iAlign), pset.nRandomSamples);

%                     quantiles = quantile(b.dataMeanRandomized{iAlign}, ...
%                         [pset.dataIntervalQuantileLow, pset.dataIntervalQuantileHigh], 4);
%                     b.dataIntervalLow{iAlign} = quantiles(:, :, :, 1);
%                     b.dataIntervalHigh{iAlign} = quantiles(:, :, :, 2);
                end
                
                if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized)
                    b.dataDifferenceOfTrialsScaledNoiseEstimateRandomized = TensorUtils.linearCombinationAlongDimension(...
                        pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, 1, abs(proj.decoderKbyN), ...
                        'replaceNaNWithZero', true, ...
                        'keepNaNIfAllNaNs', true, ...
                        'normalizeCoefficientsByNumNonNaN', false);
                end
            end

            % aggregate AlignSummary data. Each projected basis samples trials from all original
            % trials, so we aggregate all AlignSummary instances into one
            b.alignSummaryData = pset.alignSummaryAggregated';
            b.basisAlignSummaryLookup = ones(proj.nBasesProj, 1); % needs to match number of projected bases
            
            % ensure there is no translation normalization by default
            b.translationNormalization = [];
            
            b.basisValid = proj.basisValidProj;
            b.basisInvalidCause = proj.basisInvalidCauseProj;
            
            if pset.simultaneous && ~isempty(pset.dataByTrial)
                psetProjected = b.buildManualWithSingleTrialData();
            else
                psetProjected = b.buildManualWithTrialAveragedData();
            end
            
            if p.Results.applyTranslationNormalizationPostProject && ~isempty(proj.translationNormalizationPostProject)
                psetProjected = psetProjected.applyTranslationNormalization(proj.translationNormalizationPostProject);
            end

            % mark output bases invalid if requested
            if ~isempty(proj.basisValidProj)
                psetProjected = psetProjected.markBasesPermanentlyInvalid(~proj.basisValidProj, proj.basisInvalidCauseProj(~proj.basisValidProj));
            end
            
            if nargout > 1
                stats = StateSpaceProjectionStatistics.build(proj, pset, 'meanSubtract', p.Results.meanSubtractForStatistics, ...
                    'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations, ...
                    'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
                    'combineAxesWithTime', proj.combineAxesWithTime, ...
                    p.Unmatched);
            end
        end
        
        function projManual = getAsManual(proj)
            proj.warnIfNoArgOut(nargout);
            projManual = ProjManual.copyFromProjection(proj);
        end
        
        function iproj = getInverse(proj, varargin)
            p = inputParser();
%             p.addParameter('clearBeforeApplyingTranslationNormalization', true, @islogical);
%             p.addParameter('applyTranslationNormalization', true, @islogical);
%             p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            iproj = proj.getAsManual();
            
            % this is the pseudo inverse of the decoder matrix
%             iproj.decoderKbyN = proj.decoderKbyN' * (proj.decoderKbyN * proj.decoderKbyN')^(-1);
%             iproj.encoderNbyK = iproj.decoderKbyN' * (iproj.decoderKbyN * iproj.decoderKbyN')^(-1);

            % just swap the encoder and decoder
            iproj.decoderKbyN = proj.encoderNbyK;
            iproj.encoderNbyK = iproj.decoderKbyN;
            iproj.basisNamesProj = proj.basisNamesSource; % restore the basis names of the source pset when projecting back into the source bases
            iproj.basisNamesSource = proj.basisNamesProj;
            % the inverse projection will invert the post project
            % trans/norm before inverting the projection itself
            if ~isempty(proj.translationNormalizationPostProject)
                iproj.translationNormalization = proj.translationNormalizationPostProject.getInverse();
            else
                iproj.translationNormalization = [];
            end
            
            % and then after projecting, it will invert the original
            % translation normalization, so that the original data is
            % returned
            if ~isempty(proj.translationNormalization)
                iproj.translationNormalizationPostProject = proj.translationNormalization.getInverse();
            else
                iproj.translationNormalizationPostProject = [];
            end
            iproj.basisValidProj = proj.basisValid;
            iproj.basisValid = proj.basisValidProj;
            iproj.basisInvalidCauseProj = proj.basisInvalidCause;
            iproj.basisInvalidCause = proj.basisInvalidCauseProj;
        end
        
        function nproj = getProjectionToNullSpace(proj, varargin)
            nproj = proj.getAsManual();
            
            nproj.decoderKbyN = null(proj.decoderKbyN)';
            if isempty(nproj.decoderKbyN)
                error('Projection decoder matrix is full rank. Cannot project into null space.');
            end
            nproj.encoderNbyK = nproj.decoderKbyN';
            
            nproj.basisNamesProjStem = ['Null ' nproj.basisNamesProjStem];
            nproj.basisValidProj = truevec(nproj.nBasesProj);
            nproj.basisInvalidCauseProj = cellstrvec(nproj.nBasesProj);
            nproj.translationNormalizationPostProject = [];
        end
        
        function [proj, psetProjected, stats] = buildFromAndProjectPopulationTrajectorySet(proj, pset, varargin)
            p = inputParser;
            p.KeepUnmatched = true; 
            p.parse(varargin{:});

            if nargout < 3
                % runs much faster when stats not requested, though note
                % that this will leave .buildStats empty
                [proj, ~, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, 'computeStatistics', false, p.Unmatched);
            else
                [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, p.Unmatched);
            end
            psetProjected = proj.projectPopulationTrajectorySet(psetPrepared, 'applyTranslationNormalization', false, p.Unmatched);
        end
        
        function [psetInOut, statsIn] = projectInAndOut(proj, pset, varargin)
            p = inputParser();
            p.addParameter('clearBeforeApplyingTranslationNormalization', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
           
            % project into PC space, then back out to denoise
            if nargout > 1 % only compute stats if they're requested
                [psetIn, statsIn] = proj.projectPopulationTrajectorySet(pset, p.Results);
            else
                psetIn = proj.projectPopulationTrajectorySet(pset, p.Results);
            end
            iproj = proj.getInverse();
            psetInOut = iproj.projectPopulationTrajectorySet(psetIn);
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

           % pset = pset.meanSubtractBases();
        end
        
        function tf = testIsOrthogonal(proj)
            assert(proj.initialized, 'Call after building / initializing');
            
            thresh = 1e-10;
            dp = proj.decoderKbyN * proj.decoderKbyN';
            dp = abs(dp - diag(diag(dp)));
            tf = max(dp(:)) < thresh;
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
    
    % transformations that return ProjManual copy
    methods
        function proj = filterOutputBases(proj, idx)
            % select on output bases
            proj.warnIfNoArgOut(nargout);
            proj = proj.getAsManual();
            proj.decoderKbyN = proj.decoderKbyN(idx, :); % select on output bases
            proj.encoderNbyK = proj.encoderNbyK(:, idx); % select on output bases
            proj.basisValidProj = proj.basisValidProj(idx);
            proj.basisInvalidCauseProj = proj.basisInvalidCauseProj(idx);
            if ~isempty(proj.basisNamesProjManual)
                proj.basisNamesProj = proj.basisNamesProjManual(idx);
            end
            if ~isempty(proj.translationNormalizationPostProject)
                proj.translationNormalizationPostProject = proj.translationNormalizationPostProject.filterBases(idx);
            end
        end
        
        function proj = truncateOutputBases(proj, K)
            proj.warnIfNoArgOut(nargout);
            proj = proj.getAsManual();
            proj = proj.filterOutputBases(1:K);
        end
        
        function proj = reorderOutputBases(proj, idx)
            proj.warnIfNoArgOut(nargout);
            assert(numel(idx) == proj.nBasesProj, 'Number of output bases must not change, use filterOutputBases instead');
            proj = proj.filterOutputBases(idx);
        end

        function proj = orthonormalize(proj)
            proj.warnIfNoArgOut(nargout);
            proj = proj.getAsManual();
            proj = proj.orthonormalize();
        end
        
        function proj = orthogonalize(proj)
            proj.warnIfNoArgOut(nargout);
            proj = proj.getAsManual();
            proj = proj.orthogonalize();
        end
        
        function proj = normalize(proj)
            proj.warnIfNoArgOut(nargout);
            proj = proj.getAsManual();
            proj = proj.normalize();
        end
    end
    
    
    methods(Static)
        function [marginalizationNames, combinedParams, axisIncludeList, marginalizationList] = determineMarginalizationNames(pset, varargin)
            % other args include
            % 'combineAxesWithTime'
            % 'axesCombineAllMarginalizations', proj.axesCombineAllMarginalizations, ...
            % 'axesCombineSpecificMarginalizations', proj.axesCombineSpecificMarginalizations);
            
            p = inputParser();
            p.addParameter('axesCombineSpecificMarginalizations', {}, @(x) true);
            p.addParameter('axesCombineAllMarginalizations', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('combineAxesWithTime', true, @(x) islogical(x) || iscell(x));
            p.parse(varargin{:});
            args = p.Results;
            
            nConditionsAlongAxis = pset.conditionDescriptor.conditionsSize;
            dimMask = nConditionsAlongAxis > 1; % filter for non-singular axes
            
            [combinedParams, marginalizationNames, axisIncludeList, marginalizationList] = StateSpaceProjectionStatistics.generateCombinedParamsForMarginalization( ...
                pset.conditionDescriptor.axisAttributes, ...
                'axisIncludeMask', dimMask, ...
                'axisNames', pset.conditionDescriptor.axisNames, ...
                'combineAxesWithTime', args.combineAxesWithTime, ...
                'axesCombineAllMarginalizations', args.axesCombineAllMarginalizations, ...
                'axesCombineSpecificMarginalizations', args.axesCombineSpecificMarginalizations);
        end
        
        function proj = concatenate(varargin)
            projCell = varargin;
            
            % proj = concatenate(pset, projCell)
            decoders = cellfun(@(p) p.decoderKbyN, projCell, 'UniformOutput', false);
            decoderKbyN = cat(1, decoders{:});
            encoders = cellfun(@(p) p.encoderNbyK, projCell, 'UniformOutput', false);
            encoderNbyK = cat(2, encoders{:});
            
            basisValidProjCell = cellfun(@(p) p.basisValidProj, projCell, 'UniformOutput', false);
            basisValidProj = cat(1, basisValidProjCell{:});
            basisInvalidCauseProjCell = cellfun(@(p) p.basisInvalidCauseProj, projCell, 'UniformOutput', false);
            basisInvalidCauseProj = cat(1, basisInvalidCauseProjCell{:});
            basisNamesProjCell = cellfun(@(p) p.basisNamesProj, projCell, 'UniformOutput', false);
            basisNamesProj = cat(1, basisNamesProjCell{:});
            
            trNormPreCell = cellfun(@(p) p.translationNormalization, projCell, 'UniformOutput', false);
            trNormPostCell = cellfun(@(p) p.translationNormalizationPostProject, projCell, 'UniformOutput', false);
            
            proj = ProjManual.buildFromEncoderDecoder(projCell{1}, encoderNbyK, decoderKbyN);
            proj.translationNormalization = StateSpaceTranslationNormalization.concatenate(trNormPreCell);
            proj.translationNormalizationPostProject = StateSpaceTranslationNormalization.concatenate(trNormPostCell);
            proj.basisValidProj = basisValidProj;
            proj.basisInvalidCauseProj = basisInvalidCauseProj;
            proj.basisNamesProj = basisNamesProj;
        end

        function proj = compose(varargin)
            % compose several projections into one, i.e. create a single
            % projection that performs the same action as projecting by
            % varargin{1}, then varargin{2}, etc.
            %
            % This can become a bit tricky as the translationNormalization
            % (pre projection) and translationNormalizationPostProject are
            % folded into the decoder matrices as well. To make sense of
            % this code, note that translationNormalizations apply the
            % translation (addition) step first, followed by the
            % normalization.
            projCell = varargin;
            
            % We project using projCell{1} first. We'll asume its
            % translationNormalization as the translationNormalization for the
            % composed projection object, so we can ignore it subsequently.
            trPre = projCell{1}.translationNormalization;
            decoder = projCell{1}.decoderKbyN;
            bias = zerosvec(projCell{1}.nBasesProj); % can ignore bias imposed by first translation, since it will actually be applied to data
            
            % apply translationNormalizationPostProject
            tr = projCell{1}.translationNormalizationPostProject;
            if ~isempty(tr)
                bias = tr.normalizationByBasisNonNaN .* (bias + tr.translationByBasisNonNaN);
                decoder = diag(tr.normalizationByBasis) * decoder;
            end
          
            for iProj = 2:numel(projCell)
                % first handle pre-project translationNormalization
                tr = projCell{iProj}.translationNormalization;
                if ~isempty(tr)
                    bias = tr.normalizationByBasisNonNaN .* (bias + tr.translationByBasisNonNaN);
                    decoder = diag(tr.normalizationByBasisNonNaN) * decoder;
                end
                
                % then handle decoder projection
                decoder = projCell{iProj}.decoderKbyN * decoder;
                bias = projCell{iProj}.decoderKbyN * bias;
                
                if iProj < numel(projCell)
                    % then post-project translationNormalization
                    tr = projCell{1}.translationNormalizationPostProject;
                    if ~isempty(tr)
                        bias = tr.normalizationByBasis .* (bias + tr.translationByBasis);
                        decoder = diag(tr.normalizationByBasis) * decoder;
                    end
                end
            end
            
            % and then combine the last projections
            % translationNormalizationPostProject with the accumulated bias
            trPost = projCell{end}.translationNormalizationPostProject;
            if isempty(trPost)
                trPost = StateSpaceTranslationNormalization.buildIdentityManual(projCell{end}.nBasesProj);
            end
            trPost.translationByBasis = trPost.translationByBasis + bias;
            % and then final normalization will apply after the bias
            
            % now we handle the encoder side
            % first translationNormalizationPostProject is undone, which
            % handles all biases except projCell{1}'s, and projCell{end}'s
            % post project normalization, so the first encoder step should
            % undo projCell{end}'s encoder matrix and pre-project
            % normalization
            
            encoder = projCell{end}.encoderNbyK;
            tr = projCell{end}.translationNormalization;
            if ~isempty(tr)
                encoder = diag(1./tr.normalizationByBasis) * encoder;
            end
            
            for iProj = numel(projCell)-1:-1:1
                % undo the post project normalization
                tr = projCell{iProj}.translationNormalizationPostProject;
                if ~isempty(tr)
                    encoder = diag(1./tr.normalizationByBasis) * encoder;
                end
                
                % then apply the encoder matrix
                encoder = projCell{iProj}.encoderNbyK * encoder;
                
                % then undo the pre-project normalization
                 % projCell{1}'s pre-project translationNormalization is
                 % maintained by the projection and will be undone directly
                if iProj > 1
                    tr = projCell{iProj}.translationNormalization;
                    if ~isempty(tr)
                        encoder = diag(1./tr.normalizationByBasis) * encoder;
                    end
                end
            end
            
            proj = ProjManual.buildFromEncoderDecoder(projCell{1}, encoder, decoder);
            proj.translationNormalization = trPre;
            proj.translationNormalizationPostProject = trPost;
            proj.basisNamesProj = projCell{end}.basisNamesProj;
            proj.basisValidProj = projCell{end}.basisValidProj;
            proj.basisInvalidCauseProj = projCell{end}.basisInvalidCauseProj;
        end
    end
end
