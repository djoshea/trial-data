classdef StateSpaceTranslationNormalization
% This class stores the offsets and multipliers which can be applied to a
% PopulationTrajectorySet. This constitutes a normalization (subtraction off
% of each basis) followed by a normalization (multiplication of each basis
% by a scalar).

    properties(SetAccess=protected)
        % nBases x 1 numeric vector of offsets to ADD to each basis
        translationByBasis
        
        % string description of the translation (e.g. 'mean-subtracted')
        translationDescription
        
        % nBases x 1 numeric vector of offsets to ADD to each basis
        normalizationByBasis
        
        % string description of the normalization (e.g. 'variance-normalized')
        normalizationDescription
    end
    
    properties(Dependent)
        nBases
        
        isTranslation
        
        isNormalization
    end
    
    methods(Access=protected)
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
        
        function obj = assertValid(obj)
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
        
    end
    
    methods % subclasses might wish to override
        function desc = getDescription(obj)
            % generate a textual description of the applied translation and
            % normalization
            if isempty(obj.translationDescription)
                if isempty(obj.normalizationDescription)
                    desc = '';
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
            % convert nBases x 1 cellstr basisUnits to nBases x 1 cellstr
            % convertedBasisUnits
            if ~obj.isTranslation && ~obj.isNormalization
                convertedBasisUnits = basisUnits;
            else
                convertedBasisUnits = cellfun(@(s) sprintf('norm %s', s), ...
                    basisUnits, 'UniformOutput', false);
            end
        end
        
        function obj = filterBases(obj, mask)
            obj.warnIfNoArgOut(nargout);
            obj.translationByBasis = obj.translationByBasis(mask);
            obj.normalizationByBasis = obj.normalizationByBasis(mask);
        end
        
        function obj = combineWith(obj, varargin)
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
    end
        
    methods(Static)
        % subclasses should redefine this method to construct their own
        % normalizer
        function obj = buildFromPopulationTrajectorySet(pset)
            obj = StateSpaceTranslationNormalization.buildIdentityForPopulationTrajectorySet(pset);
        end
        
        function obj = buildIdentityForPopulationTrajectorySet(pset)
            obj = StateSpaceTranslationNormalization.buildManual(...
                zeros(pset.nBases, 1), ones(pset.nBases, 1));
        end
    end
       
    methods(Static, Sealed)
        function obj = buildManual(varargin)
            p = inputParser();
            p.addRequired('translationByBasis', @isvector);
            p.addRequired('normalizationByBasis', @isvector);
            p.addParamValue('translationDescription', 'manual offset', @ischar);
            p.addParamValue('normalizationDescription', 'manual normalization', @ischar);
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
                newData = bsxfun(@rdivide, bsxfun(@plus, data, obj.translationByBasis), obj.normalizationByBasis);
            end
        end
        
        function newData = applyTranslationToData(obj, data)
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, offset) (x + offset), data, ...
                    num2cell(obj.translationByBasis), ...
                    'UniformOutput', false);
            else
                newData = bsxfun(@plus, data, obj.translationByBasis);
            end
        end
        
        function newData = applyNormalizationToData(obj, data)
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, mult) x * mult, data, ...
                    num2cell(obj.normalizationByBasis), ...
                    'UniformOutput', false);
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
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, offset) (x - offset), data, ...
                    num2cell(obj.translationByBasis), ...
                    'UniformOutput', false);
            else
                newData = bsxfun(@minus, data, obj.translationByBasis);
            end
        end
        
        function newData = undoNormalizationToData(obj, data)
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, mult) x ./ mult, data, ...
                    num2cell(obj.normalizationByBasis), ...
                    'UniformOutput', false);
            else
                newData = bsxfun(@rdivide, data, obj.normalizationByBasis);
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