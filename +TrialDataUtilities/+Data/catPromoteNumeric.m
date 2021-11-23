function out = catPromoteNumeric(dim, varargin)

empty = cellfun(@(x) isempty(x), varargin);
types = cellfun(@(x) string(class(x)), varargin);

if any(empty) && ~all(empty)
    indNonEmpty = find(~empty, 1, 'first');
    types(empty) = types(indNonEmpty);
end

if numel(unique(types)) == 1
    out = cat(dim, varargin{~empty});
else
    if any(types == "double")
        outtype = "double";
    elseif any(types == "single")
        outtype = "single";
    else
        % all integer types
        if any(endsWith(types, "64"))
            bits = 64;
        elseif any(endsWith(types, "32"))
            bits = 32;
        elseif any(endsWith(types, "16"))
            bits = 16;
        elseif any(endsWith(types, "8"))
            bits = 8;
        else
            error('Not sure what type to use for concatenation');
        end
        
        if any(startsWith(types, "int"))
            outtype = "int" + string(bits);
        else
            assert(all(startsWith(types, "uint")));
            outtype = "uint" + string(bits);
        end
    end
    
    varargin = cellfun(@(x) cast(x, outtype), varargin(~empty), 'UniformOutput', false);
    out = cat(dim, varargin{:});
end