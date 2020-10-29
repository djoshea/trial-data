function ceval = interpColorMap(cmap, nvals, varargin)
    ceval = TrialDataUtilities.Color.evalColorMapAt(cmap, linspace(0, 1, nvals), [0 1], varargin{:});
end