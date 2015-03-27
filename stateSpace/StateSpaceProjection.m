classdef StateSpaceProjection 
% Abstract base class representing a set of coefficients per basis with which
% to project PopulationTrajectorySet instances. Typically these are built using 
% fromPopulationTrajectorySet to compute the coefficients, in a manner determined
% by subclasses

    properties(SetAccess=protected)
        initialized = false;
        translationNormalization % stored translation / normalization stored when building and used when projecting
        
        % given an N x T matrix of data in the original basis, we project
        % into the new basis using decoderKbyN * neural_NbyT which yields
        % decoded_KbyT. We can then reconstruct the data in the original
        % basis using reconstruction_NbyT = encoderNbyK * decodedKbyT, or
        % reconstruction_NbyT = encoderNbyK * decoderKbyN * neural_NbyT 
        decoderKbyN
        encoderNbyK 
        
        buildStats % StateSpaceProjectionStatistics instance from build
        
        basisValid % N x 1 logical vector indicating which bases were considered valid in the projection (typically copied from the pset from which I am built)
        basisInvalidCause % N x 1 cell str (typically copied from the pset from which I am built)
    end

    properties(Dependent)
        nBasesSource
        nBasesProj
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
    end
    
    methods(Abstract)
        % return a list of basis names for the new basis
        names = getBasisNames(proj, pset, data)

        % compute the N * K matrix of basis coefficients for the projection
        % given an N x T matrix of data in the original basis, we project
        % into the new basis using decoderKbyN * neural_NbyT which yields
        % decoded_KbyT. We can then reconstruct the data in the original
        % basis using reconstruction_NbyT = encoderNbyK * decodedKbyT, or
        % reconstruction_NbyT = encoderNbyK * decoderKbyN * neural_NbyT 
        [decoderKbyN, encoderNbyK] = computeProjectionCoefficients(pset)
    end

    methods
        function [proj, stats] = buildFromPopulationTrajectorySet(proj, pset, varargin)
            % build this projection matrix based on an existing PopulationTrajectorySet
            % defers to calculateProjectionMatrix for the actual basis computation
            
            proj.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('pset', @(x) isa(x, 'PopulationTrajectorySet'));
            p.addParameter('axisCombinations', {}, @iscell);
            p.addParameter('combineCovariatesWithTime',  true, @islogical);
            p.parse(pset, varargin{:});
            
            % make any necessary transformations, particularly translation / normalization
            pset = proj.preparePsetForInference(pset);

            % extract the translation normalization that will be applied before projection
            proj.translationNormalization = pset.translationNormalization;

            % compute the coefficients for the projection
            debug('Computing projection encoder and decoder coefficients\n');
            [proj.decoderKbyN, proj.encoderNbyK] = proj.computeProjectionCoefficients(pset);
            
            assert(size(proj.decoderKbyN, 2) == pset.nBases, 'Decoder matrix returned by computeProjectionCoefficients must match pset.nBases along dim 2');
            assert(size(proj.encoderNbyK, 1) == pset.nBases, 'Encoder matrix returned by computeProjectionCoefficients must match pset.nBases along dim 1');
            
            % copy the basis valid mask
            proj.basisValid = pset.basisValid;
            proj.basisInvalidCause = pset.basisInvalidCause;
            
            % set coeff to zero on invalid bases
            proj.decoderKbyN(:, ~proj.basisValid) = 0;
            proj.encoderNbyK(~proj.basisValid, :) = 0;
            
            % results will store statistics and useful quantities related to the
            % projection
            stats = StateSpaceProjectionStatistics.build(proj, pset, ...
                'combineCovariatesWithTime', p.Results.combineCovariatesWithTime, ...
                'axisCombinations', p.Results.axisCombinations);
            proj.buildStats = stats;

            proj.initialized = true;
        end

        function [psetProjected, stats] = projectPopulationTrajectorySet(proj, pset, varargin)
            p = inputParser();
            p.addParameter('applyTranslationNormalization', true, @islogical);
            p.addParameter('combineCovariatesWithTime',  true, @islogical);
            p.addParameter('axisCombinations', {}, @iscell);
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
            if p.Results.applyTranslationNormalization
                debug('Applying translation/normalization associated with projection to data\n');
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

                projMat = proj.decoderKbyN * mat;
                b.dataMean{iAlign} = reshape(projMat, proj.nBasesProj, ...
                    pset.nConditions, pset.nTimeDataMean(iAlign));

                % use sqrt(sd1^2 / n1 + sd2^2 / n2 + ...) formula
                % which equals sqrt(|coeff1| * sem1^2 + |coeff2| * sem2^2 + ...)
                mat = reshape(pset.dataSem{iAlign}, pset.nBases, ...
                    pset.nConditions * pset.nTimeDataMean(iAlign));
                mat(~proj.basisValid, :) = 0;
                projMat = sqrt(abs(proj.decoderKbyN) * (mat.^2));
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
                    projMat = proj.decoderKbyN * mat;
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
            if nargout > 1
                stats = StateSpaceProjectionStatistics.build(proj, pset, ...
                    'combineCovariatesWithTime', p.Results.combineCovariatesWithTime, ...
                    'axisCombinations', p.Results.axisCombinations);
            end
        end
        
        function [proj, psetProjected, stats] = buildFromAndProjectPopulationTrajectorySet(proj, pset, varargin)
            p = inputParser;
            p.addParameter('axisCombinations', {}, @iscell);
            p.addParameter('combineCovariatesWithTime',  true, @islogical);
            p.parse(varargin{:});

            [proj, stats] = proj.buildFromPopulationTrajectorySet(pset, p.Results);
            psetProjected = proj.projectPopulationTrajectorySet(pset, 'applyTranslationNormalization', false, p.Results); % false --> don't apply translation normalization since we'll pull that from this pset to begin with
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
        function proj = filterOutputBases(proj, idx)
            % select on output bases
            proj.warnIfNoArgOut(nargout);
            assert(proj.initialized, 'Call filterBases after building / initializing');
            proj.decoderKbyN = proj.decoderKbyN(idx, :); % select on output bases
            proj.encoderNbyK = proj.encoderNbyK(:, idx); % select on output bases
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
    
end
