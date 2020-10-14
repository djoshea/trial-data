function h = identityLine(varargin)

p = inputParser();
p.addParameter('slope', 1, @isscalar);
p.addParameter('offset', 0, @isscalar);
p.KeepUnmatched = true;
p.parse(varargin{:});

xl = xlim();

xd = diff(xl);
xv = xl;
xv(1) = xv(1) - 5*xd;
xv(2) = xv(2) + 5*xd;
yv = p.Results.slope*xv + p.Results.offset;

h = plot(xv, yv, ...
    'Color', [0.5 0.5 0.5], p.Unmatched, 'XLimInclude', 'off', 'YLimInclude', 'off');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

end
