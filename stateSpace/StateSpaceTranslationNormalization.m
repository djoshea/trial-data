classdef StateSpaceTranslationNormalization
% This class stores the offsets and multipliers which can be applied to a
% PopulationTrajectorySet. This constitutes a translation (subtraction off
% of each basis) followed by a normalization (multiplication of each basis
% by a scalar)

    properties
        % nBases x 1 numeric vector of offsets to ADD to each basis
        translationByBasis
        
        % string description of the translation (e.g. 'mean-subtracted')
        translationDescription = '';
        
        % nBases x 1 numeric vector of offsets to MULTIPLY to each basis
        normalizationByBasis
        
        % string description of the normalization (e.g. 'variance-normalized')
        normalizationDescription = '';
    end
    
    properties(Dependent,SetAccess=protected)
        nBases
        
        isTranslation
        
        isNormalization

        % this is a shortcut to ~isnan(translationByBasis) & ~isnan(normalizationByBasis)
        % it's primary use is to ensure consistency when taken from one
        % pset and applied to another pset which doesn't have the same set
        % of baes marked valid 
        basisValid 
        
        translationByBasisNonNaN % replace NaN with 0
        normalizationByBasisNonNaN % replace NaN with 0
    end
    
    methods(Access=protected) % don't call constructor, use factory methods
        function obj = StateSpaceTranslationNormalization()
            % should be built using factory methods below
        end
    end
    
    methods 
        function v = get.nBases(obj)
            v = numel(obj.translationByBasis);
        end
        
        function tf = get.isTranslation(obj)
            tf = any(obj.translationByBasis ~= 0);
        end
        
        function tf = get.isNormalization(obj)
            tf = any(obj.normalizationByBasis ~= 1);
        end
        
        function mask = get.basisValid(obj)
            mask = ~isnan(obj.translationByBasis) & ~isnan(obj.normalizationByBasis);
        end
            
        function t = get.translationByBasisNonNaN(obj)
            t = obj.translationByBasis;
            t(isnan(t)) = 0;
        end
           
        function n = get.normalizationByBasisNonNaN(obj)
            n = obj.normalizationByBasis;
            n(isnan(n)) = 0;
        end
    end
    
    methods % methods a subclass might wish to override
        function obj = assertValid(obj)
            % if overriding be sure to call this superclass method
            obj.warnIfNoArgOut(nargout);
            
            assert(isempty(obj.translationByBasis) || (isvector(obj.translationByBasis) && isnumeric(obj.translationByBasis)), ...
                'translationByBasis must be numeric vector');
            assert(isempty(obj.normalizationByBasis) || (isvector(obj.normalizationByBasis) && isnumeric(obj.normalizationByBasis)), ...
                'normalizationByBasis must be numeric vector');
            assert(length(obj.normalizationByBasis) == length(obj.translationByBasis), ...
                'normalizationByBasis and translationByBasis must have the same length (nBases)');
            
            assert(ischar(obj.normalizationDescription), ...
                'normalizationDescription must be string');
            assert(ischar(obj.translationDescription), ...
                'translationDescription must be string');
            
            assert(~any(obj.normalizationByBasis == 0), 'normalizationByBasis cannot contain 0 values');
            
            obj.translationByBasis = makecol(obj.translationByBasis);
            obj.normalizationByBasis = makecol(obj.normalizationByBasis);
        end
        
        function desc = getDescription(obj)
            % generate a textual description of the applied translation and
            % normalization
            if isempty(obj.translationDescription)
                if isempty(obj.normalizationDescription)
                    desc = 'none';
                else
                    desc = obj.normalizationDescription;
                end
            elseif isempty(obj.normalizationDescription)
                desc = obj.translationDescription;
            else
                desc = sprintf('%s, %s', obj.translationDescription, obj.normalizationDescription);
            end
        end
        
        function printOneLineDescription(obj)
            tcprintf('inline', '{yellow}%s: {none}%s\n', ...
                class(obj), obj.getDescription());
        end
        
        function disp(obj)
            obj.printOneLineDescription();
            fprintf('\n');
            builtin('disp', obj);
        end
        
        function convertedBasisUnits = convertBasisUnits(obj, basisUnits)
            % used to specify a new name for each basis if you'd like
            %
            % convert nBases x 1 cellstr basisUnits to nBases x 1 cellstr
            % convertedBasisUnits. Be sure to return a column vector
            % cellstr
            if ~obj.isTranslation && ~obj.isNormalization
                convertedBasisUnits = basisUnits;
            else
                convertedBasisUnits = cellfun(@(s) sprintf('normalized %s', s), ...
                    basisUnits, 'UniformOutput', false);
            end
        end
    end
    
    methods(Sealed)
        function obj = filterBases(obj, mask)
            % be sure to call this superclass method if overriding
            obj.warnIfNoArgOut(nargout);
            obj.translationByBasis = obj.translationByBasis(mask);
            obj.normalizationByBasis = obj.normalizationByBasis(mask);
        end
        
        function obj = setBasesInvalid(obj, mask)
            obj.warnIfNoArgOut(nargout);
            obj.translationByBasis(mask) = NaN;
            obj.normalizationByBasis(mask) = NaN;
        end
        
        function obj = combineWith(obj, varargin)
            % be sure to call this superclass method if overriding and to
            % not break the basic behavior
            for i = 1:numel(varargin)
                new = varargin{i};
                obj.translationByBasis = obj.translationByBasis + new.translationByBasis;
                obj.normalizationByBasis = obj.normalizationByBasis .* new.normalizationByBasis;
                if isempty(obj.translationDescription)
                    obj.translationDescription = new.translationDescription;
                elseif ~isempty(new.translationDescription)
                        obj.translationDescription = [obj.translationDescription, ', ', new.translationDescription];
                end
                if isempty(obj.normalizationDescription)
                    obj.normalizationDescription = new.normalizationDescription;
                elseif ~isempty(new.normalizationDescription)
                    obj.normalizationDescription = [obj.normalizationDescription, ', ', new.normalizationDescription];
                end
            end
        end
        
        function tr = getTranslationOnly(obj)
            tr = StateSpaceTranslationNormalization(); 
            tr.translationByBasis = obj.translationByBasis;
            tr.normalizationByBasis = ones(obj.nBases, 1);
        end

        function inv = getInverse(obj)
            obj.warnIfNoArgOut(nargout);
            % in the inverse, the order of translation, then normalization
            % persists, so we must derive the inverse as follows.
            % 
            % the forward transformation looks like y = f(x) = (x-t) / n.
            % the inverse transform looks like:
            % x = f'(y) 
            %   = (y-ti) * ni 
            %   = ((x-t)*n - ti) * ni
            %   = (x*n - t*n - ti) * ni
            % Letting ti = -t*n and ni = 1/n, we have:
            %   = (x*n - t*n - (-t*n)) * (1/n)
            %   = (x*n) / n
            %   = x
                        
            inv = StateSpaceTranslationNormalization(); 
            inv.translationByBasis = -obj.translationByBasis .* obj.normalizationByBasis;
            inv.normalizationByBasis = obj.normalizationByBasis.^(-1);
        
            % string description of the normalization (e.g. 'variance-normalized')
            if ~isempty(obj.translationDescription)
                inv.translationDescription = ['invert ' obj.translationDescription];
            end
            if ~isempty(obj.normalizationDescription)
                inv.normalizationDescription = ['invert ' obj.normalizationDescription];
            end
        end
    end
        
    methods(Static)
        % subclasses should redefine this method to construct their own
        % normalizer
        function obj = buildFromPopulationTrajectorySet(pset)
            obj = StateSpaceTranslationNormalization.buildIdentityForPopulationTrajectorySet(pset);
        end
    end
       
    methods(Static, Sealed)
        function obj = buildIdentityForPopulationTrajectorySet(pset)
            obj = StateSpaceTranslationNormalization.buildManual(...
                zeros(pset.nBases, 1), ones(pset.nBases, 1));
        end
        
        function obj = buildIdentityManual(nBases)
            obj = StateSpaceTranslationNormalization.buildManual(...
                zeros(nBases, 1), ones(nBases, 1));
        end
        
        function obj = buildManual(varargin)
            p = inputParser();
            p.addRequired('translationByBasis', @isvector);
            p.addRequired('normalizationByBasis', @isvector);
            p.addParameter('translationDescription', 'manual offset', @ischar);
            p.addParameter('normalizationDescription', 'manual normalization', @ischar);
            p.parse(varargin{:});
            
            obj = StateSpaceTranslationNormalization();
            
            obj.translationByBasis = p.Results.translationByBasis;
            if all(obj.translationByBasis == 0)
                obj.translationDescription = '';
            else
                obj.translationDescription = p.Results.translationDescription;
            end
            
            obj.normalizationByBasis = p.Results.normalizationByBasis;
            if all(obj.normalizationByBasis == 1)
                obj.normalizationDescription = '';
            else
                obj.normalizationDescription = p.Results.normalizationDescription;
            end
            
            obj = obj.assertValid();
        end
        
        function obj = concatenate(trNormCell)
            % be sure to call this superclass method if overriding and to
            % not break the basic behavior
            
            emptyMask = cellfun(@isempty, trNormCell);
            if all(emptyMask)
                obj = [];
                return;
            end
            trNormCell = trNormCell(~emptyMask);
            
            translations = cellfun(@(t) t.translationByBasis, trNormCell, 'UniformOutput', false);
            translationByBasis = cat(1, translations{:});
            norms = cellfun(@(t) t.normalizationByBasis, trNormCell, 'UniformOutput', false);
            normByBasis = cat(1, norms{:});
            
            obj = StateSpaceTranslationNormalization.buildManual(translationByBasis, normByBasis);
        end
        
        function tf = checkEqual(varargin)
            tf = true;
            first = varargin{1};
            for i = 2:numel(varargin)
                if ~all(TrialDataUtilities.Data.isequaltol(first.translationByBasis, varargin{i}.translationByBasis)) || ...
                   ~all(TrialDataUtilities.Data.isequaltol(first.normalizationByBasis, varargin{i}.normalizationByBasis))
                    tf = false;
                    continue;
                end
            end
        end
        
        function obj = combine(varargin)
            obj = varargin{1}.combineWith(varargin(2:end));
        end
    end
    
    methods(Sealed) % apply normalization and/or normalization to data as a vector or cell array
        function newData = applyTranslationNormalizationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                offCell = num2cell(repmat(obj.translationByBasis, sz));
                normCell = num2cell(repmat(obj.normalizationByBasis, sz));
                newData = cellfun(@(x, offset, mult) (x + offset) * mult, data, ...
                    offCell, normCell, 'UniformOutput', false);
            else
                newData = bsxfun(@times, bsxfun(@plus, data, obj.translationByBasis), obj.normalizationByBasis);
            end
        end
        
        function newData = applyTranslationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                offCell = num2cell(repmat(obj.translationByBasis, sz));
                newData = cellfun(@(x, offset) x + offset, data, ...
                    offCell, 'UniformOutput', false);
            else
                newData = bsxfun(@plus, data, obj.translationByBasis);
            end
        end
        
        function newData = applyNormalizationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                normCell = num2cell(repmat(obj.normalizationByBasis, sz));
                newData = cellfun(@(x, norm) x .* norm, data, ...
                    normCell, 'UniformOutput', false);
            else
                newData = bsxfun(@times, data, obj.normalizationByBasis);
            end
        end
        
        function newData = undoTranslationNormalizationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                offCell = num2cell(repmat(obj.translationByBasis, sz));
                normCell = num2cell(repmat(obj.normalizationByBasis, sz));
                newData = cellfun(@(x, offset, mult) (x ./ mult) - offset, data, ...
                    offCell, normCell, 'UniformOutput', false);
            else
                newData = bsxfun(@minus, bsxfun(@rdivide, data, obj.normalizationByBasis), obj.translationByBasis);
            end
        end
        
        function newData = undoTranslationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                offCell = num2cell(repmat(obj.translationByBasis, sz));
                newData = cellfun(@(x, offset) x - offset, data, ...
                    offCell, 'UniformOutput', false);
            else
                newData = bsxfun(@minus, data, obj.translationByBasis);
            end
        end
        
        function newData = undoNormalizationToData(obj, data)
            assert(size(data, 1) == obj.nBases, ...
                'Data must be nBases along dimension 1');
            sz = size(data);
            sz(1) = 1;
            
            if iscell(data)
                normCell = num2cell(repmat(obj.normalizationByBasis, sz));
                newData = cellfun(@(x, norm) x ./ norm, data, ...
                    normCell, 'UniformOutput', false);
            else
                newData = bsxfun(@times, data, obj.normalizationByBasis);
            end
        end
    end
    
    methods(Sealed, Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~ishandle(obj)
                message = sprintf('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
                expr = sprintf('debug(''%s'')', message);
                evalin('caller', expr); 
            end
        end
    end
    
end
