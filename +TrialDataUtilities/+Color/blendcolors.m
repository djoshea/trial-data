function cmap = blendcolors(colors, n, varargin)
    p = inputParser();
    p.addParameter('inputSpace', 'rgb', @ischar);
    p.addParameter('interval', [0 1], @isvector);
    p.KeepUnmatched = false;
    p.parse(varargin{:});
    interval = p.Results.interval;
    
    colorLuv = TrialDataUtilities.Color.convert(sprintf('%s->Luv', p.Results.inputSpace), colors);
    
    cmapLuv = interp1(linspace(0, 1, size(colors, 1)), colorLuv, linspace(interval(1), interval(2), n));
    
    cmap = TrialDataUtilities.Color.convert('Luv->rgb', cmapLuv);

end
