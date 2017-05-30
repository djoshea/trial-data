function cmap = blendcolors(colors, n, varargin)
    p = inputParser();
    p.addParameter('interval', [0 1], @isvector);
    p.KeepUnmatched = false;
    p.parse(varargin{:});
    interval = p.Results.interval;
    
    cmap = interp1(linspace(0, 1, size(colors, 1)), colors, linspace(interval(1), interval(2), n));

end
