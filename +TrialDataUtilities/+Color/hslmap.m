function map = hslmap(n, hfirst, sat, lum)
% map = hslmap(n, sat, lum)
    
    if nargin < 2
        hfirst = 0.01;
    end
    if nargin < 3
        sat = 0.6;
    end
    if nargin < 4
        lum = 0.65;
    end

    hsl = [circspace(hfirst, 360, n)', sat*ones(n, 1), lum*ones(n,1)];
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