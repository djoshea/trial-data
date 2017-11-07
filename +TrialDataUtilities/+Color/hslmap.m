function map = hslmap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('saturation', 1, @isscalar);
    p.addParameter('luminance', 0.6, @isscalar);
    p.addParameter('hueSpan', [-0.01 0.8], @isvector);
    p.parse(varargin{:});
    
    if nargin < 1
        n = 100;
    end
    
    hueSpan = p.Results.hueSpan;
    if mod(hueSpan(1), 1) == mod(hueSpan(2), 1)
        hues = 360 * mod(circspace(hueSpan(1), hueSpan(2), n)', 1);   
    else
        hues = 360 * mod(linspace(hueSpan(1), hueSpan(2), n)', 1);  
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