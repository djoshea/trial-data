function ceval = evalColorMapAt(cmap, at, clim)
    if isa(cmap, 'function_handle')
        cmap = cmap(1000);
    end
    if nargin < 3
        clim = [0 1];
    end
    
    clim = cast(clim, 'like', cmap);
    if clim(2) == clim(1)
        clim(2) = clim(2) + 0.1;
    end
    
    N = size(cmap, 1);
    at = cast(at, 'like', cmap);
    at = clamp(at, clim(1), clim(2));
    
    cmapOrigAt = (0:N-1) / (N-1) * (clim(2)-clim(1)) + clim(1);
    cmap = TrialDataUtilities.Color.convert('RGB->LCH', cmap);
    ceval = interp1(cmapOrigAt, cmap, at);
    ceval = TrialDataUtilities.Color.convert('LCH>RGB', ceval);
    
    ceval = max(min(ceval, 1), 0);
end

function val = clamp(val, lo, hi)
    val(val < lo) = lo;
    val(val > hi) = hi;
end

