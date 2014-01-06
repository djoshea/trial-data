classdef StateSpaceTranslationNormalization
% This class stores the offsets and multipliers which can be applied to a
% PopulationTrajectorySet. This constitutes a translation (subtraction off
% of each basis) followed by a normalization (multiplication of each basis
% by a scalar).

    properties
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
    end
    
    methods
        function obj = StateSpaceTranslationNormalization()
            obj.translationDescription = 'manual offset';
            obj.normalizationDescription = 'manual normalization';
        end
        
        function v = get.nBases(obj)
            v = numel(obj.translationByBasis);
        end
        
        function desc = getDescription(obj)
            % generate a textual description of the applied translation and
            % normalization
            desc = sprintf('%s, %s', obj.translationDescription, obj.normalizationDescription);
        end
        
        function convertedBasisUnits = convertBasisUnits(basisUnits)
            % convert nBases x 1 cellstr basisUnits to nBases x 1 cellstr
            % convertedBasisUnits
            convertedBasisUnits = cellfun(@(s) sprintf('norm %s', s), ...
                basisUnits, 'UniformOutput', false);
        end
    end
    
     methods(Static)
        function obj = buildFromPopulationTrajectorySet(pset)
            % determine the coefficients for the translation and
            % normalization using data in a PopulationTrajectorySet
            % this will be more usefully defined in subclasses
            obj.translationByBasis = zerovec(pset.nBases);
            obj.normalizationByBasis = onesvec(pset.nBases);
        end
    end
    
    methods % apply translation and/or normalization to data as a vector or cell array
        function newData = applyTranslationNormalizationToData(obj, data)
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, offset, mult) (x + offset) * mult, data, ...
                    num2cell(obj.translationByBasis), num2cell(obj.normalizationByBasis), ...
                    'UniformOutput', false);
            else
                newData = (data + obj.translationByBasis) .* obj.normalizationByBasis;
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
                newData = data + obj.translationByBasis;
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
                newData = data .* obj.normalizationByBasis;
            end
        end
        
        function newData = undoTranslationNormalizationToData(obj, data)
            assert(isvector(data) && numel(data) == obj.nBases, ...
                'Data must be nBases vector');
            if iscell(data)
                newData = cellfun(@(x, offset, mult) (x ./ mult) - offset, data, ...
                    num2cell(obj.translationByBasis), num2cell(obj.normalizationByBasis), ...
                    'UniformOutput', false);
            else
                newData = (data ./ obj.normalizationByBasis) - obj.translationByBasis;
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
                newData = data - obj.translationByBasis;
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
                newData = data ./ obj.normalizationByBasis;
            end
        end
    end
    
   
    
end