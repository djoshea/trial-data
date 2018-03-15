function map = hslmap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('saturation', 0.7, @isscalar);
    p.addParameter('luminance', 0.65, @isscalar);
    p.addParameter('hueSpan', 360*[-0.01 0.8], @isvector);
    p.parse(varargin{:});
    
    if nargin < 1
        n = 100;
    end
    
    hueSpan = p.Results.hueSpan;
    if mod(hueSpan(1), 1) == mod(hueSpan(2), 1)
        hues = mod(circspace(hueSpan(1), hueSpan(2), n)', 360);   
    else
        hues = mod(linspace(hueSpan(1), hueSpan(2), n)', 360);
    end
     
    hsl = [hues, p.Results.saturation * ones(n, 1), p.Results.luminance * ones(n,1)];
    map = colorspace('HSL->RGB', hsl);
end

function v = circspace(d1, d2, n)
% v = circspace(d1, d2, n)
% like linspace, except considers d1 == d2 in a circular axis

    if nargin == 2
        n = 100;
    else
        n = floor(double(n));
    end

    delta = (d2-d1)/n;
    v = linspace(d1, d2-delta, n);

end