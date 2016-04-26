function h = identityLine(varargin)

xl = xlim();

xd = diff(xl);
xv = xl;
xv(1) = xv(1) - 5*xd;
xv(2) = xv(2) + 5*xd;
yv = xv;

h = plot(xv, yv, ...
    'Color', [0.5 0.5 0.5], varargin{:}, 'XLimInclude', 'off', 'YLimInclude', 'off');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

end