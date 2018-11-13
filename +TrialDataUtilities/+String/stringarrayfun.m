function out = stringarrayfun(varargin)

    outCell = arrayfun(varargin{:}, 'UniformOutput', false);
    empty = cellfun(@isempty, outCell);

    outCell(empty) = {""};

    out = string(outCell);

end