function map = hslmap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('fracHueShift', 0.001, @issscalar); % between 0 and 1
    p.addParameter('saturation', 0.7, @isscalar);
    p.addParameter('luminance', 0.65, @isscalar);
    p.addParameter('fracHueSpan', 1, @isscalar);
    p.parse(varargin{:});
    
    hues = 360 * mod(p.Results.fracHueShift + circspace(0, p.Results.fracHueSpan, n)', 1);    
    hsl = [hues, p.Results.saturation * ones(n, 1), p.Results.luminance * ones(n,1)];
    map = colorspace('HSL->RGB', hsl);
end

