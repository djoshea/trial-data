function ceval = evalColorMapAt(cmap, at)
    if isa(cmap, 'function_handle')
        cmap = cmap(1000);
    end
    N = size(cmap, 1);
    ceval = interp1((0:N-1) / (N-1), cmap, at);
end