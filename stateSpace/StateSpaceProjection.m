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
        decoderKbyN
        encoderNbyK
        
        buildStats % StateSpaceProjectionStatistics instance from build
        
        basisValid % nBasesSource x 1 logical vector indicating which bases were considered valid in the projection (typically copied from the pset from which I am built)
        basisInvalidCause % nBasesSource x 1 cell str (typically copied from the pset from which I am built)
        
        % nBasesProj x 1 logical vector indicating which OUTPUT bases will be considered valid
        % this is mostly useful when building inverses, so that you get the
        % original basis valid mask back when inverting a projection
        basisValidProj 
        basisInvalidCauseProj
        
        basisNamesSource % nBasesSource x 1 cellstr of basis names from the pset I was built off
        dataUnitsSource = '';
        
        % nBasesProj x 1 cell of basis names, optional, that I will pick for the basis names of the projected 
        % output. If set, this will override the method "getBasisNames"
        % which typically picks these names based on the input pset.
        % this is primarily used for inverse projections, when we want to
        % restore the basis names of the pset this projection was built off 
        basisNamesProj
    end

    properties(Dependent)
        nBasesSource
        nBasesProj
        
        decoderKbyNValid
        encoderNbyKValid
    end
    
    methods(Abstract, Static)
        % this abstract static method is basically a way of getting around matlab's lack
        % of support for the factory builder pattern. Define this method
        % like this so that buildFromPopulationTrajectorySet() below can do
        % the heavy lifting.
        % methods(Static)
        %     function [proj, stats, psetPrepared] = createFrom(pset, varargin)
        %         proj = YOURSUBCLASSNAME() you can send in varargin{:} if there are options > );
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
        % return a nBasesProj x 1 cell str: list of basis names for the new basis  
        function names = getBasisNames(proj, pset) %#ok<INUSD>
            if ~isempty(proj.basisNamesProj)
                names = proj.basisNamesProj;
            else
                names = arrayfun(@(i) sprintf('Basis %d', i), ...
                        (1:proj.nBasesProj)', 'UniformOutput', false);
            end
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
        
        function d = get.decoderKbyNValid(proj)
            d = proj.decoderKbyN(proj.basisValidProj, proj.basisValid);
        end
        
        function d = get.encoderNbyKValid(proj)
            d = proj.encoderNbyK(proj.basisValid, proj.basisValidProj);
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
    end

    methods
        function [proj, stats, psetPrepared] = buildFromPopulationTrajectorySet(proj, pset, varargin)
            % build this projection matrix based on an existing PopulationTrajectorySet
            % defers to calculateProjectionMatrix for the actual basis computation
            
            proj.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('pset', @(x) isa(x, 'PopulationTrajectorySet'));
            p.addParameter('computeStatistics', true, @islogical);
            p.addParameter('meanSubtractForStatistics', true, @islogical); % mean subtract data when computing statistics, this makes sense to turn off if the data is already measured relative to some meaningful baseline
            p.addParameter('nBasesProj', NaN, @isscalar); % this is optional for both the caller to provide and the projection to obey
            p.KeepUnmatched = true;
            p.parse(pset, varargin{:});
            
            % make any necessary transformations, particularly translation / normalization
            pset = proj.preparePsetForInference(pset);
            psetPrepared = pset;

            % extract the translation normalization that will be applied before projection
            proj.translationNormalization = pset.translationNormalization;

            % compute the coefficients for the projection
            debug('Computing projection encoder and decoder coefficients\n');
            [decoderKbyN, encoderNbyK, proj] = proj.computeProjectionCoefficients(pset, 'nBasesProj', p.Results.nBasesProj);
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
            proj.basisNamesProj = proj.getBasisNames(pset);
            
            % set coeff to zero on invalid bases
            proj.decoderKbyN(:, ~proj.basisValid) = 0;
            proj.encoderNbyK(~proj.basisValid, :) = 0;
            
            % results will store statistics and useful quantities related to the
            % projection
            if p.Results.computeStatistics
                stats = StateSpaceProjectionStatistics.build(proj, pset, 'meanSubtract', p.Results.meanSubtractForStatistics, p.Unmatched);
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

            % ensure proj-invalid bases are marked invalid
            % this isn't strictly necessary but just in case
            pset = pset.setBasesInvalid(~proj.basisValid, 'invalidated before state space projection');
            
            if any(proj.basisValid & ~pset.basisValid)
                error('PopulationTrajectorySet has invalid bases not marked invalid in StateSpaceProjection. You should equalize the bases invalid to get consistent results');
            end
            
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
            b.dataUnits = ''; % clear by default, since we're not necessarily willing to call them the same units anymore

            if ~isempty(proj.basisNamesProj)
                b.basisNames = proj.basisNamesProj;
            else
                b.basisNames = proj.getBasisNames(pset);
            end
            
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
            end

            % aggregate AlignSummary data. Each projected basis samples trials from all original
            % trials, so we aggregate all AlignSummary instances into one
            b.alignSummaryData = pset.alignSummaryAggregated';
            b.basisAlignSummaryLookup = ones(pset.nBases, 1);
            
            % ensure there is no translation normalization by default
            b.translationNormalization = [];
            
            psetProjected = b.buildManualWithTrialAveragedData();
            
            if p.Results.applyTranslationNormalizationPostProject && ~isempty(proj.translationNormalizationPostProject)
                psetProjected = psetProjected.applyTranslationNormalization(proj.translationNormalizationPostProject);
            end

            % mark output bases invalid if requested
            if ~isempty(proj.basisValidProj)
                psetProjected = psetProjected.setBasesInvalid(~proj.basisValidProj, proj.basisInvalidCauseProj(~proj.basisValidProj));
            end
            
            if nargout > 1
                stats = StateSpaceProjectionStatistics.build(proj, pset, 'meanSubtract', p.Results.meanSubtractForStatistics, p.Unmatched);
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
            iproj.decoderKbyN = proj.encoderNbyK;
            iproj.encoderNbyK = proj.decoderKbyN;
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
            proj.basisNamesProj = proj.basisNamesProj(idx);
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
    
end
