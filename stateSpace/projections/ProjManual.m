classdef ProjManual < StateSpaceProjection
    
    properties
        decoderKbyNManual
        encoderNbyKManual
        meanCenterBases = true;
    end
    
    methods
        function v = get.decoderKbyNManual(proj)
            if ~isempty(proj.decoderKbyN)
                v = proj.decoderKbyN;
            else
                v = proj.decoderKbyNManual;
            end
        end
        
        function v = get.encoderNbyKManual(proj)
            if ~isempty(proj.encoderNbyK)
                v = proj.encoderNbyK;
            else
                v = proj.encoderNbyKManual;
            end
        end
        
        function proj = set.decoderKbyNManual(proj, v)
            % save the old nBasesProj
%             nBasesProjOrig = proj.nBasesProj;
            
            proj.decoderKbyNManual = v;
            proj.decoderKbyN = v;
            
            % c
        end
        
        function proj = set.encoderNbyKManual(proj, v)
            proj.encoderNbyKManual = v;
            proj.encoderNbyK = v;
        end
    end
    
    methods(Static)
        function proj = buildFromEncoderDecoder(psetOrProj, encoderNbyK, decoderKbyN, varargin)
            proj = ProjManual.createFrom(psetOrProj, 'encoderNbyK', encoderNbyK, 'decoderKbyN', decoderKbyN, varargin{:});
        end
        
        function proj = buildSubspaceProjectionFromDecoder(psetOrProj, decoderKbyN, varargin)
            % choose encoder such that encoder * decoder forms the
            % projection matrix onto the row space of the decoder
            [Q, ~] = TrialDataUtilities.Data.qrGramSchmidt(decoderKbyN');
            decoderKbyN = Q';
            encoderNbyK = decoderKbyN'; % * (decoderKbyN * decoderKbyN')^(-1);
            proj = ProjManual.buildFromEncoderDecoder(psetOrProj, encoderNbyK, decoderKbyN, varargin{:});
        end
      
        function proj = buildWithPseudoinverseEncoderForDecoderAndPset(pset, projOrDecoderKbyN, varargin)
            % This builds an encoder matrix based on the equation:
            % (X * D') * E' = X
            % where X and Xhat is CTA x N, D is K x N, E is N x K
            % where E will be estimated by least squares:
            % E' = (X*D') \ X
            % E = ((X*D') \ X)'
            if isa(projOrDecoderKbyN, 'StateSpaceProjection')
                decoderKbyN = projOrDecoderKbyN.decoderKbyN;
            elseif isnumeric(projOrDecoderKbyN)
                decoderKbyN = projOrDecoderKbyN;
            else
                error('Second arg must be StateSpaceProjection or decoderKbyN');
            end
            
            Xv = pset.buildCTAbyN('validBasesOnly', true, 'validTimepointsAllBasesOnly', true);
            
            decoderKbyNv = decoderKbyN(:, pset.basisValid);
            encoderNvbyK = ((Xv*decoderKbyNv') \ (Xv))';

            encoderNbyK = TensorUtils.inflateMaskedTensor(encoderNvbyK, 1, pset.basisValid, NaN);
            
            proj = ProjManual.buildFromEncoderDecoder(pset, encoderNbyK, decoderKbyN, varargin{:});
        end
        
        function proj = buildIdentityFor(psetOrProj)
            if isa(psetOrProj, 'StateSpaceProjection')
                nBases = psetOrProj.nBasesSource;
            else
                nBases = psetOrProj.nBases;
            end
            proj = ProjManual.buildInternal(psetOrProj, 'meanCenterBases', false, ...
                'decoderKbyN', eye(nBases), 'encoderNbyK', eye(nBases));
        end
        
        function [proj] = createFrom(pset, varargin)
            % you can provide the properties listed above as param value
            % pairs
            [proj, ~] = ProjManual.buildInternal(pset, varargin{:});
%             [proj, stats, psetPrepared] = proj.buildFromPopulationTrajectorySet(pset, unmatched);
        end

        function [proj, psetProjected, stats] = createFromAndProject(pset, varargin)
            % you can provide the properties listed above as param value
            % pairs
            [proj, unmatched] = ProjManual.buildInternal(pset, varargin{:});
            [psetProjected, stats] = proj.projectPopulationTrajectorySet(pset, unmatched);
        end
    end
    
    methods(Static, Access=protected)
        function [proj, unmatched] = buildInternal(psetOrProj, varargin)
            if isa(psetOrProj, 'StateSpaceProjection')
                nBases = psetOrProj.nBasesSource;
                proj = ProjManual.copyFromProjection(psetOrProj);
            else
                nBases = psetOrProj.nBases;
                proj = ProjManual.initializeFromPopulationTrajectorySet(psetOrProj);
            end
            
            p = inputParser();
            p.addParamValue('meanCenterBases', true, @islogical);
            p.addParamValue('decoderKbyN', eye(nBases), @(x) ~isempty(x) && ismatrix(x));
            p.addParamValue('encoderNbyK', eye(nBases), @(x) ~isempty(x) && ismatrix(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            unmatched = p.Unmatched;
            
            proj.meanCenterBases = p.Results.meanCenterBases;
            proj.decoderKbyNManual = p.Results.decoderKbyN;
            proj.encoderNbyKManual = p.Results.encoderNbyK;
     
            if isa(psetOrProj, 'PopulationTrajectorySet')
                if ~isempty(psetOrProj.translationNormalization)
                    proj.translationNormalization = psetOrProj.translationNormalization;
                    proj.meanCenterBases = false;
                end
            end
            
            assert(size(proj.decoderKbyN, 2) == nBases, 'Decoder matrix returned by computeProjectionCoefficients must match nBases along dim 2');
            assert(size(proj.encoderNbyK, 1) == nBases, 'Encoder matrix returned by computeProjectionCoefficients must match nBases along dim 1');

            % set coeff to zero on invalid bases
            proj.decoderKbyN(:, ~proj.basisValid) = 0;
            proj.encoderNbyK(~proj.basisValid, :) = 0;
           
            proj.initialized = true;
        end
    end
    
    methods(Static)
        function proj = copyFromProjection(projTemplate, varargin)
            proj = ProjManual();
            
            % copy all properties wholesale!
            meta = ?StateSpaceProjection;
            props = meta.PropertyList;
            for iProp = 1:length(props)
                prop = props(iProp);
                if prop.Dependent || prop.Constant || prop.Transient
                    continue;
                else
                    name = prop.Name;
                    proj.(name) = projTemplate.(name);
                end
            end
            
            proj.decoderKbyNManual = proj.decoderKbyN;
            proj.encoderNbyKManual = proj.encoderNbyK;
        end
        
        function proj = initializeFromPopulationTrajectorySet(pset)
            proj = ProjManual();
            proj.translationNormalization = pset.translationNormalization;
            proj.basisValid = pset.basisValid;
            proj.basisInvalidCause = pset.basisInvalidCause;
            proj.basisNamesSource = pset.basisNames;
            proj.dataUnitsSource = pset.dataUnits;
        end
        
        function proj = concatenateBasesFromProjections(projCell)
            % build from the first projection
            proj = ProjManual.copyFromProjection(projCell{1});
            
            nSource = cellfun(@(proj) proj.nBasesSource, projCell);
            assert(numel(unique(nSource)) == 1, 'Projections differ on .nBasesSource');
            
            encodersNbyK = cellfun(@(proj) proj.encoderNbyK, projCell, 'UniformOutput', false);
            decodersKbyN = cellfun(@(proj) proj.decoderKbyN, projCell, 'UniformOutput', false);
            
            catEncoder = cat(2, encodersNbyK{:});
            catDecoder = cat(1, decodersKbyN{:});
            
            proj.decoderKbyNManual = catDecoder;
            proj.encoderNbyKManual = catEncoder;
            % not necessary due to set. methods
%             proj.decoderKbyN = proj.decoderKbyNManual;
%             proj.encoderNbyK = proj.encoderNbyKManual;
            
            temp = cellfun(@(proj) proj.basisValidProj, projCell, 'UniformOutput', false);
            proj.basisValidProj = cat(1, temp{:});
            temp = cellfun(@(proj) proj.basisInvalidCauseProj, projCell, 'UniformOutput', false);
            proj.basisInvalidCauseProj = cat(1, temp{:});
            temp = cellfun(@(proj) proj.getBasisNames(), projCell, 'UniformOutput', false);
            proj.basisNamesProj = cat(1, temp{:});
        end
    end
    
    methods
        function proj = ProjManual(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
        function pset = preparePsetForInference(proj, pset) 
            if proj.meanCenterBases
                pset = pset.meanSubtractBases();
            end
        end
        
        function proj = filterOutputBases(proj, idx)
            % select on output bases
            proj.warnIfNoArgOut(nargout);
            proj = filterOutputBases@StateSpaceProjection(proj, idx);
            proj.decoderKbyNManual = proj.decoderKbyNManual(idx, :); % select on output bases
            proj.encoderNbyKManual = proj.encoderNbyKManual(:, idx); % select on output bases
        end
        
        function proj = reorderOutputBases(proj, idx)
            proj.warnIfNoArgOut(nargout);
            assert(numel(idx) == proj.nBasesProj, 'Number of output bases must not change, use filterOutputBases instead');
            proj = proj.filterOutputBases(idx);
        end
        
        function proj = orthogonalizeDecoderRows(proj)
            % this will orthonogonalize the decoder matrix row vectors and set the
            % encoder matrix to maintain the same encoder * decoder
            % reconstruction
            proj.warnIfNoArgOut(nargout);
            
            % use modified gram schmidt
            % Q is N x K, R is K x K
            [Q, R] = TrialDataUtilities.Data.qrGramSchmidt(proj.decoderKbyNManual');
            
            % then take the norms of each column from the diagonal terms
            % from R, and move them into Q and out of R
            normsFromR = makerow(diag(R));
            Qmod = bsxfun(@times, Q, normsFromR);
            Rmod = bsxfun(@rdivide, R, normsFromR');
            proj.decoderKbyNManual = Qmod';
            proj.encoderNbyKManual = proj.encoderNbyKManual * Rmod';
        end
        
        function proj = orthonormalize(proj)
            % this will orthonormalize the decoder matrix and set the
            % encoder matrix to maintain the same encoder * decoder
            % reconstruction
            proj.warnIfNoArgOut(nargout);

            % use modified gram schmidt, not sure if necessary
            [Q, R] = TrialDataUtilities.Data.qrGramSchmidt(proj.decoderKbyNValid');
            decoderValid = Q';
            encoderValid = proj.encoderNbyKValid * R';
            
            proj.decoderKbyNManual = TensorUtils.inflateMaskedTensor(decoderValid, [1 2], {proj.basisValidProj, proj.basisValid}, 0);
            proj.encoderNbyKManual = TensorUtils.inflateMaskedTensor(encoderValid, [1 2], {proj.basisValid, proj.basisValidProj}, 0);
        end
        
        function proj = normalize(proj)
            % this will row normalize the decoder matrix and set the
            % encoder matrix to maintain the same encoder * decoder
            % reconstruction
            proj.warnIfNoArgOut(nargout);

            normsByBasis = sqrt(sum(proj.decoderKbyNManual.^2, 2)); % sum over N neurons for each of K bases
            proj.decoderKbyNManual = bsxfun(@rdivide, proj.decoderKbyNManual, normsByBasis);
            proj.encoderNbyKManual = bsxfun(@times, proj.encoderNbyKManual, normsByBasis');
        end
        
        function proj = combineEncoderDecoder(proj)
            % this will set the decoder to be the encoder*decoder product and 
            % make the encoder matrix the identity
            proj.warnIfNoArgOut(nargout);

            proj.decoderKbyNManual = proj.encoderNbyK * proj.decoderKbyN;
            proj.encoderNbyKManual = eye(proj.nBasesSource);
        end
    end

    methods
        function [decoderKbyN, encoderNbyK, proj] = computeProjectionCoefficients(proj, pset, varargin)
            decoderKbyN = proj.decoderKbyNManual;
            encoderNbyK = proj.encoderNbyKManual;
            
            assert(size(decoderKbyN, 2) == pset.nBases, 'Decoder matrix must be K by N where N matches pset.nBases');
            assert(size(encoderNbyK, 1) == pset.nBases, 'Encoder matrix must be N by K where N matches pset.nBases');
            assert(size(decoderKbyN, 1) == size(encoderNbyK, 2), 'Encoder matrix (N x K) and decoder matrix (K x N) differ on size K');
        end
    end

end
