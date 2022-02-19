function cmap = blendcolors(colors, nOrEvalAt, varargin)
    p = inputParser();
    p.addParameter('inputSpace', 'rgb', @ischar);
    p.addParameter('locations', linspace(0, 1, size(colors, 1)), @isvector);
    p.addParameter('interval', [], @isvector);
    p.addParameter('forceEvalAt', false, @islogical);
    p.KeepUnmatched = false;
    p.parse(varargin{:});
    interval = p.Results.interval;
    locations = p.Results.locations;
    if isempty(interval)
        interval = [min(locations) max(locations)];
    end
    
    colorLuv = TrialDataUtilities.Color.convert(sprintf('%s->Luv', p.Results.inputSpace), colors);
    
    if isscalar(nOrEvalAt) && nOrEvalAt > 1 && ~p.Results.forceEvalAt
        n = nOrEvalAt;
        cmapLuv = interp1(locations, colorLuv, linspace(interval(1), interval(2), n));
    else
        evalAt = nOrEvalAt;
        cmapLuv = interp1(locations, colorLuv, evalAt);
    end
    
    cmap = TrialDataUtilities.Color.convert('Luv->rgb', cmapLuv);

end
