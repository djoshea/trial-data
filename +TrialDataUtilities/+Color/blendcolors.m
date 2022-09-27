function cmap = blendcolors(colors, nOrEvalAt, varargin)
    p = inputParser();
    p.addParameter('interpSpace', 'Luv', @ischar);
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
    
    useOkLab = strcmpi(p.Results.interpSpace, 'oklab');
    
    if useOkLab
        if ~strcmp(p.Results.inputSpace, 'rgb')
            colors = TrialDataUtilities.Color.convert(sprintf('%s->rgb', p.Results.inputSpace), colors);
        end
        colorInterpSpace = TrialDataUtilities.Color.rgb2oklab(colors);
    else
        colorInterpSpace = TrialDataUtilities.Color.convert(sprintf('%s->%s', p.Results.inputSpace, p.Results.interpSpace), colors);
    end

    if isscalar(nOrEvalAt) && nOrEvalAt > 1 && ~p.Results.forceEvalAt
        n = nOrEvalAt;
        cmapInterpSpace = interp1(locations, colorInterpSpace, linspace(interval(1), interval(2), n));
    else
        evalAt = nOrEvalAt;
        cmapInterpSpace = interp1(locations, colorInterpSpace, evalAt);
    end
    
    if useOkLab
        cmap = TrialDataUtilities.Color.oklab2rgb(cmapInterpSpace);
    else
        cmap = TrialDataUtilities.Color.convert(sprintf('%s->rgb', p.Results.interpSpace), cmapInterpSpace);
    end

end
