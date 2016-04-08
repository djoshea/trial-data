function map = hslmap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('fracHueShift', 0.001, @issscalar); % between 0 and 1
    p.addParameter('saturation', 0.6, @isscalar);
    p.addParameter('luminance', 0.65, @isscalar);
    p.addParameter('fracHueSpan', 1, @isscalar);
    p.parse(varargin{:});
    
    hues = 360 * mod(p.Results.fracHueShift + circspace(0, p.Results.fracHueSpan, n)', 1);    
    hsl = [hues, p.Results.saturation * ones(n, 1), p.Results.luminance * ones(n,1)];
    map = colorspace('HSL->RGB', hsl);
%    map = hsl2rgb(hsl);
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

function rgb=hsl2rgb(hsl_in)
    %Converts Hue-Saturation-Luminance Color value to Red-Green-Blue Color value
    %
    %Usage
    %       RGB = hsl2rgb(HSL)
    %
    %   converts HSL, a M [x N] x 3 color matrix with values between 0 and 1
    %   into RGB, a M [x N] X 3 color matrix with values between 0 and 1
    %
    %See also rgb2hsl, rgb2hsv, hsv2rgb

    % (C) Vladimir Bychkovsky, June 2008
    % written using: 
    % - an implementation by Suresh E Joel, April 26,2003
    % - Wikipedia: http://en.wikipedia.org/wiki/HSL_and_HSV

    hsl=reshape(hsl_in, [], 3);

    H=hsl(:,1);
    S=hsl(:,2);
    L=hsl(:,3);

    lowLidx=L < (1/2);
    q=(L .* (1+S) ).*lowLidx + (L+S-(L.*S)).*(~lowLidx);
    p=2*L - q;
    hk=H; % this is already divided by 360

    t=zeros([length(H), 3]); % 1=R, 2=B, 3=G
    t(:,1)=hk+1/3;
    t(:,2)=hk;
    t(:,3)=hk-1/3;

    underidx=t < 0;
    overidx=t > 1;
    t=t+underidx - overidx;

    range1=t < (1/6);
    range2=(t >= (1/6) & t < (1/2));
    range3=(t >= (1/2) & t < (2/3));
    range4= t >= (2/3);

    % replicate matricies (one per color) to make the final expression simpler
    P=repmat(p, [1,3]);
    Q=repmat(q, [1,3]);
    rgb_c= (P + ((Q-P).*6.*t)).*range1 + ...
            Q.*range2 + ...
            (P + ((Q-P).*6.*(2/3 - t))).*range3 + ...
            P.*range4;

    rgb_c=round(rgb_c.*10000)./10000; 
    rgb=reshape(rgb_c, size(hsl_in));
end
