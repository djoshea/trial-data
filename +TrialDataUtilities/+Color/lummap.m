function map = lummap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('hue', 0, @isscalar); % between 0 and 1
    p.addParameter('saturation', 0.6, @isscalar);
    p.addParameter('luminanceLimits', [0.2 0.8], @isscalar);
    p.parse(varargin{:});
    
    if ~isscalar(n)
        lumVals = makecol(n);
        n = numel(lumVals);
    else
        lims = p.Results.luminanceLimits;
        lumVals = linspace(lims(1), lims(2), n)';
    end
    hsl = [p.Results.hue * ones(n, 1), p.Results.saturation * ones(n, 1), lumVals];
    map = TrialDataUtilities.Color.hsl2rgb(hsl);
end

