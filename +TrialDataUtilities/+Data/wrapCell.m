function d = wrapCell(d)
    if ~iscell(d)
        d = {d};
    end
end