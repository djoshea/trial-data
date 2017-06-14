function equalizeAxisLimits(hax, which, varargin)

p = inputParser();
p.addOptional('sameExtent', false, @islogical);
p.parse(varargin{:});

if nargin < 2
    which = 'xy';
end
doX = ismember('x', which);
doY = ismember('y', which);
doZ = ismember('z', which);

sameExtent = p.Results.sameExtent;
if isscalar(sameExtent)
    sameExtent = repmat(sameExtent, 3, 1);
end

argsX = {};
argsY = {};
argsZ = {};

if doX
    xl = cell2mat(get(hax, 'XLim'));
    if sameExtent(1)
        xl = [min(xl(:, 1)) max(xl(:, 2))];
    else
        xm = max(abs(xl(:)));
        xl = [-xm xm];
    end
    argsX = {'XLim', xl};
end

if doY
    yl = cell2mat(get(hax, 'YLim'));
    if sameExtent(2)
        yl = [min(yl(:, 1)) max(yl(:, 2))];
    else
        ym = max(abs(yl(:)));
        yl = [-ym ym];
    end
    argsY = {'YLim', yl};
end

if doZ
    zl = cell2mat(get(hax, 'ZLim'));
    if sameExtent(3)
        zl = [min(zl(:, 1)) max(zl(:, 2))];
    else
        zm = max(abs(zl(:)));
        zl = [-zm zm];
    end
    argsZ = {'ZLim', zl};
end

set(hax, argsX{:}, argsY{:}, argsZ{:});

