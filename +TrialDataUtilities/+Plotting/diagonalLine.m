function h = diagonalLine(varargin)

xl = xlim();
yl = ylim();

xd = diff(xl);
yd = diff(yl);
if xd > yd
    lim(1) = xl(1) - 5*xd;
    lim(2) = xl(2) + 5*xd;
else
    lim(1) = yl(1) - 5*yd;
    lim(2) = yl(2) + 5*yd;
end

h = plot(lim, lim, 'Color', [0.5 0.5 0.5], varargin{:}, 'XLimInclude', 'off', 'YLimInclude', 'off');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');