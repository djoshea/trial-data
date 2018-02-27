classdef PropertyShapeMeta < matlab.mixin.CustomDisplay
    % a lightweight metadata class for describing the meaning of each dimension in a nested
    % array, i.e. one in which the top levels are cell arrays and the bottom level is a matrix,
    % For example, an A x B x C cell array of D x E matrices would be {{'A','B','C'}, {'D','E'}}
    %
    % PropertyShapeMeta can also specify custom handling of slicing and concatenation operations on
    % specific dimensions. Typically, this can be handled automatically based on the shape of the value.
    % However it is also possible to specify custom handling functions for specific dimensions, even if these dimensions
    % are not included in dimsByLevel. This is useful if the property holds an class instance that internally has data
    % along that dimension. For example, the StateSpaceTranslationNormalization class manages arrays of length 'N' internally.
    % So a property like pset.translationNormalization that holds a scalar StateSpaceTranslationNormalization can handle slicing
    % along 'N' using a custom function, e.g. 'selectBases'.
    properties
        dimsByLevel cell = {} % will be cell of cellstr, each top level corresponds to a dim
        leafClass char = 'float';
        attr struct = struct(); % contains any additional metadata applying to the whole array

        customDims cell = {}; % list of dimension names handled via some custom function
        customSelectAlongDimFn function_handle = @(meta, value, dimNames, masksByDim, varargin) error('Custom slicing function not specified');
        customConcatenateAlongDimFn function_handle = @(meta, valueCell, dimName) error('Custom concatenation function not specified');
    end

    properties(Dependent)
        depth
        flatDimNames
    end

    methods
        function meta = PropertyShapeMeta(dimsByLevel, leafClass, varargin)
            p = inputParser();
            p.addParameter('customDims', {}, @iscellstr);
            p.addParameter('customSelectAlongDimFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('customConcatenateAlongDimFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            assert(iscell(dimsByLevel) && (isvector(dimsByLevel) || isempty(dimsByLevel)));
            if iscellstr(dimsByLevel) && ~isempty(dimsByLevel)
                dimsByLevel = {dimsByLevel};
            end

            dimsByLevel = makecol(dimsByLevel);
            for i = 1:numel(dimsByLevel)
                assert(iscellstr(dimsByLevel{i}) && isvector(dimsByLevel{i}));
                dimsByLevel{i} = makecol(dimsByLevel{i});
            end
            meta.dimsByLevel = dimsByLevel;

            % check for dim name uniqueness
            dimNames = cat(1, meta.dimsByLevel{:});
            assert(numel(dimNames) == numel(unique(dimNames)));

            meta.leafClass = leafClass;

            % setup custom dim handling if specified
            meta.customDims = p.Results.customDims;
            if ~isempty(p.Results.customSelectAlongDimFn)
                meta.customSelectAlongDimFn = p.Results.customSelectAlongDimFn;
            end
            if ~isempty(p.Results.customConcatenateAlongDimFn)
                meta.customConcatenateAlongDimFn = p.Results.customConcatenateAlongDimFn;
            end

            meta.attr = p.Unmatched;
        end

        function str = getDescription(meta)
            strCell = cell(1, meta.depth);
            for iL = 1:meta.depth
                if iL < meta.depth
                    type = 'cell of ';
                elseif meta.depth > 1
                    type = 'arrays';
                else
                    type = 'array';
                end
                strCell{iL} = sprintf('%s %s', strjoin(meta.dimsByLevel{iL}, ' x '), type);
            end
            str = strjoin(strCell, '');
        end

        function v = get.depth(meta)
            v = numel(meta.dimsByLevel);
        end

        function v = get.flatDimNames(meta)
            v = cat(1, meta.dimsByLevel{:});
        end

        function [tf, level, dim] = hasDimByName(meta, dimName)
            [level, dim] = meta.findDimByName(dimName);
            tf = ~isnan(level);
        end

       function [level, dim] = findDimByName(meta, dimName)
            if ischar(dimName)
                dimName = {dimName};
            end
            N = numel(dimName);
            [level, dim] = deal(nan(N, 1));

            for iL = 1:meta.depth
                [tf, which] = ismember(dimName, meta.dimsByLevel{iL});
                level(tf) = iL;
                dim(tf) = which(tf);
            end
        end

        function value = getAttrWithDefault(meta, fld, default)
            if isfield(meta.attr, fld)
                value = meta.attr.(fld);
            else
                value = default;
            end
        end
    end

    methods % Data transformation methods
        function [data, wasUpdated] = selectAlongDimByName(meta, data, dimNames, masks)
            % selects along dimension dimName with mask
            % if dimName is char, mask is a logical or numerical vector
            % if dimName is vector, mask is a cell of mask vectors for each dimName
            % This method will call the
            if ischar(dimNames)
                dimNames = {dimNames};
                masks = {masks};
            end
            
            maskCustom = ismember(dimNames, meta.customDims);
            custom_dimNames = dimNames(maskCustom);
            custom_masks = masks(maskCustom);
            dimNames = dimNames(~maskCustom);
            masks = masks(~maskCustom);

            wasUpdated = false;
            
            if any(~maskCustom)
                [level, dim] = meta.findDimByName(dimNames);
                if any(~isnan(level))
                    data = meta.selectAlongLevelDim(data, level, dim, masks);
                    wasUpdated = true;
                end
            end
            if any(maskCustom)
                data = meta.customSelectAlongDimFn(meta, data, custom_dimNames, custom_masks);
                wasUpdated = true;
            end
        end

        function [data, wasUpdated] = selectAlongLevelDim(meta, data, level, dim, masks)
            % this method performs a slicing operation on specific levels/dims. Typically the client would call
            % selectAlongDimByName which looksup the level/dims by the name of the dimension and also calls any custom functions
            if ~iscell(masks)
                masks = {masks};
            end
            data = filterInner(data, 1);
            wasUpdated = numel(level) > 0;

            function d = filterInner(d, levelThis)
                % apply any filtering that needs to happen at this level, and then recurse on deeper levels
                idx = find(level == levelThis);
                if any(idx)
                    d = selectAlongDimension(d, dim(idx), masks(idx));
                end
                if levelThis < meta.depth && any(level > levelThis)
                    for iD = 1:numel(d)
                        d{iD} = filterInner(d{iD}, levelThis+1);
                    end
                end
            end

            function t = selectAlongDimension(t, dim, select)
                sz = size(t);
                if numel(sz) < numel(dim)
                    sz = cat(2, sz, ones(1, numel(dim) - numel(sz)));
                end
                maskByDim = cell(numel(sz), 1);
                for iD = 1:numel(sz)
                    which = find(dim == iD, 1);
                    if ~isempty(which)
                        maskByDim{iD} = select{which};
                    else
                        maskByDim{iD} = true(sz(iD), 1);
                    end
                end
                t = t(maskByDim{:});
            end
        end

        function data = assignValueMaskedSelectionAlongDimByName(meta, data, dimNames, masks, value, clearCells)
            [level, dim] = meta.findDimByName(dimNames);
            data = meta.assignValueMaskedSelectionAlongLevelDim(data, level, dim, masks, value, clearCells);
        end

        function data = assignValueMaskedSelectionAlongLevelDim(meta, data, level, dim, masks, value, clearCells)
            % for a specific dimension at the specified level and dim, assign value into positions selected by mask
            % if level is the last level, value will be used. If level is not the last level (i.e. this level is cells of a deeper level)
            % clearCells==true will replace those cell elements contents with {}
            % clearCells==false will replace those cell elements contents at the last level with value

            assert(isscalar(value));
            if ~iscell(masks)
                masks = {masks};
            end

            data = assignInner(data, 1);

            function d = assignInner(d, levelThis)
                % apply any filtering that needs to happen at this level, and then recurse on deeper levels
                idx = find(level == levelThis);

                if levelThis < meta.depth
                    % not the last level
                    if any(idx)
                        masksByDim = buildMasksByDim(d, dim(idx), masks(idx));
                        if clearCells
                            % we'll be setting these positions in the cell to empty
                            [d{masksByDim{:}}] = deal({});
                        else
                            d(masksByDim{:}) = setFinalLevelToValue(d(masksByDim{:}), value);
                        end
                    end
                    if any(level > levelThis)
                        for iD = 1:numel(d)
                            d{iD} = assignInner(d{iD}, levelThis+1);
                        end
                    end
                else
                    if any(idx)
                        masksByDim = buildMasksByDim(d, dim(idx), masks(idx));
                        d(masksByDim{:}) = value;
                    end
                end
            end

            function maskByDim = buildMasksByDim(t, dim, select)
                sz = size(t);
                if numel(sz) < numel(dim)
                    sz = cat(2, sz, ones(1, numel(dim) - numel(sz)));
                end
                maskByDim = cell(numel(sz), 1);
                for iD = 1:numel(sz)
                    which = find(dim == iD, 1);
                    if ~isempty(which)
                        maskByDim{iD} = select{which};
                    else
                        maskByDim{iD} = true(sz(iD), 1);
                    end
                end
            end

            function t = setFinalLevelToValue(t, value)
                if iscell(t)
                    for iT = 1:numel(t)
                        t{iT} = setFinalLevelToValue(t{iT}, value);
                    end
                else
                    t(:) = value;
                end
            end
        end

        function data = applyFnToLastLevel(meta, data, fn)
            data = inner(data, 1);

            function d = inner(d, levelThis)
                if levelThis == meta.depth
                    d = fn(d);
                else
                    for i = 1:numel(d)
                        d{i} = inner(d{i}, levelThis+1);
                    end
                end
            end
        end

        function out = emptyFromSizes(meta, sizeStruct, fill)
            if nargin < 3
                fill = NaN;
            end

            function sz = sizesFromDimNames(dimNames)
                sz = nan(numel(dimNames), 1);
                for iD = 1:numel(sz)
                    sz(iD) = sizeStruct.(dimNames{iD});
                end
                if numel(sz) == 1
                    sz(2) = 1;
                end
            end

            prev = repmat(fill, sizesFromDimNames(meta.dimsByLevel{meta.depth})');
            for iL = meta.depth-1:-1:1
                this = cell(sizesFromDimNames(meta.dimsByLevel{iL})');
                this(:) = {prev};
                prev = this;
            end

            out = prev;
        end

        function vals = sizeValuesFromData(meta, data)
            vals = struct();
            vals = inner(vals, data, 1);

            function vals = inner(vals, d, level)
                for iD = 1:numel(meta.dimsByLevel{level})
                    vals.(meta.dimsByLevel{level}{iD}) = size(d, iD);
                end
                if level < meta.depth
                    vals = inner(vals, d{1}, level+1);
                end
            end
        end
    end

    methods(Static)
        function sizes = deepSize(data)
            if iscell(data)
                sizes = cat(1, {size(data)}, PropertyShapeMeta.deepSize(data{1}));
            else
                sizes = {size(data)};
            end
        end
        
        function metaStruct = filterMetaStruct(metaStruct, keepFn, nested)
            % given a struct consisting of other structs and PropertyShapeMeta objects, keep those fields where
            % keepFn(propMeta, prop) evaluates to true
            if nargin < 3
                nested = false;
            end
            metaStruct = TrialDataUtilities.Struct.filterFields(metaStruct, keepFn, nested);
        end
    end

    methods(Access=protected)
        function header = getHeader(meta)
            if ~isscalar(meta)
                header = getHeader@matlab.mixin.CustomDisplay(meta);
            else
                header = sprintf('%s: %s\n  attr: %s', class(meta), meta.getDescription(), structToString(meta.attr, ' '));
            end

            function str = structToString(s, separator)
                % given a struct s with string or numeric vector values, convert to string
                fields = fieldnames(s);
                if isempty(fields)
                    str = '';
                    return;
                end
                vals = structfun(@convertToString, s);

                str = strjoin(cellfun(@(fld, val) [fld '=' val], fields, vals, 'UniformOutput', false), separator);

                return;

                function str = convertToString(v)
                    if ischar(v)
                        str = v;
                    elseif isempty(v)
                        str = '[]';
                    elseif isnumeric(v) || islogical(v)
                        str = mat2str(v);
                    elseif iscellstr(v)
                        str = ['{', strjoin(v, ','), '}'];
                    elseif isobject(v)
                        if ismethod(v, 'char')
                            str = char(v);
                        elseif ismethod(v, 'describe')
                            str = describe(v);
                        else
                            str = class(v);
                        end
                    else
                        error('Could not convert struct field value');
                    end

                    str = {str};
                end

            end
        end
    end

end
