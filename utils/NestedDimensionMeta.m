classdef NestedDimensionMeta < matlab.mixin.CustomDisplay
    % a lightweight metadata class for describing the meaning of each dimension in a nested
    % array, i.e. one in which the top levels are cell arrays and the bottom level is a matrix,
    % For example, an A x B x C cell array of D x E matrices would be
    properties
        dimsByLevel cell = {} % will be cell of cellstr, each top level corresponds to a dim
        attr struct = struct(); % contains any additional metadata applying to the whole array
    end

    properties(Dependent)
        depth
        flatDimNames
    end

    methods
        function ndm = NestedDimensionMeta(dimsByLevel, varargin)
            p = inputParser();
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
            ndm.dimsByLevel = dimsByLevel;

            % check for dim name uniqueness
            dimNames = cat(1, ndm.dimsByLevel{:});
            assert(numel(dimNames) == numel(unique(dimNames)));

            ndm.attr = p.Unmatched;
        end
        
        function str = getDescription(ndm)
            strCell = cell(1, ndm.depth);
            for iL = 1:ndm.depth
                if iL < ndm.depth
                    type = 'cell of ';
                elseif ndm.depth > 1
                    type = 'arrays';
                else
                    type = 'array';
                end
                strCell{iL} = sprintf('%s %s', strjoin(ndm.dimsByLevel{iL}, ' x '), type);
            end
            str = strjoin(strCell, '');
        end

        function v = get.depth(ndm)
            v = numel(ndm.dimsByLevel);
        end

        function v = get.flatDimNames(ndm)
            v = cat(1, ndm.dimsByLevel{:});
        end
        
        function [tf, level, dim] = hasDimByName(ndm, dimName)
            [level, dim] = ndm.findDimByName(dimName);
            tf = ~isnan(level);
        end

       function [level, dim] = findDimByName(ndm, dimName)
            if ischar(dimName)
                dimName = {dimName};
            end
            N = numel(dimName);
            [level, dim] = deal(nan(N, 1));
            
            for iL = 1:ndm.depth
                [tf, which] = ismember(dimName, ndm.dimsByLevel{iL});
                level(tf) = iL;
                dim(tf) = which(tf);
            end
        end

        function data = filterDimByName(ndm, data, dimName, masks)
            % selects along dimension dimName with mask
            % if dimName is char, mask is a logical or numerical vector
            % if dimName is vector, mask is a cell of mask vectors for each dimName
            [level, dim] = ndm.findDimByName(dimName);
            data = ndm.filterLevelDim(data, level, dim, masks);
        end
        
        function data = filterLevelDim(ndm, data, level, dim, masks)
            if ~iscell(masks)
                masks = {masks};
            end
            data = filterInner(data, 1);
            
            function d = filterInner(d, levelThis)
                % apply any filtering that needs to happen at this level, and then recurse on deeper levels
                idx = find(level == levelThis);
                if any(idx)
                    d = selectAlongDimension(d, dim(idx), masks(idx));
                end
                if levelThis < ndm.depth && any(level > levelThis)
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
        
        function data = assignValueMaskedSelectionAlongDimByName(ndm, data, dimName, masks, value, clearCells)
            [level, dim] = ndm.findDimByName(dimName);
            data = ndm.assignValueMaskedSelectionAlongLevelDim(data, level, dim, masks, value, clearCells);
        end
        
        function data = assignValueMaskedSelectionAlongLevelDim(ndm, data, level, dim, masks, value, clearCells)
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
                
                if levelThis < ndm.depth 
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
        
        function data = applyFnToLastLevel(ndm, data, fn)
            data = inner(data, 1);
            
            function d = inner(d, levelThis)
                if levelThis == ndm.depth
                    d = fn(d);
                else
                    for i = 1:numel(d)
                        d{i} = inner(d{i}, levelThis+1);
                    end
                end
            end 
        end

        function out = emptyFromSizes(ndm, sizeStruct, fill)
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
            
            prev = repmat(fill, sizesFromDimNames(ndm.dimsByLevel{ndm.depth})');
            for iL = ndm.depth-1:-1:1
                this = cell(sizesFromDimNames(ndm.dimsByLevel{iL})');
                this(:) = {prev};
                prev = this;
            end
            
            out = prev;
        end
        
        function vals = sizeValuesFromData(ndm, data)
            vals = struct();
            vals = inner(vals, data, 1);
            
            function vals = inner(vals, d, level)
                for iD = 1:numel(ndm.dimsByLevel{level})
                    vals.(ndm.dimsByLevel{level}{iD}) = size(d, iD);
                end
                if level < ndm.depth
                    vals = inner(vals, d{1}, level+1);
                end
            end
        end
    end
    
    methods(Static)
        function sizes = deepSize(data)
            if iscell(data)
                sizes = cat(1, {size(data)}, NestedDimensionMeta.deepSize(data{1}));
            else
                sizes = {size(data)};
            end
        end
    end
    
    methods(Access=protected)
        function header = getHeader(ndm)
            if ~isscalar(ndm)
                header = getHeader@matlab.mixin.CustomDisplay(ndm);
            else
                header = sprintf('%s: %s\n  attr: %s', class(ndm), ndm.getDescription(), structToString(ndm.attr, ' '));
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
