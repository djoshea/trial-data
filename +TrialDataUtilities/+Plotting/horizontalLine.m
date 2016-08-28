function h = horizontalLine(yVal, varargin)

xl = xlim();

xd = diff(xl);
xl(1) = xl(1) - 5*xd;
xl(2) = xl(2) + 5*xd;

h = gobjects(numel(yVal), 1);
for i = 1:numel(yVal)
    h(i) = plot(xl, [yVal(i), yVal(i)], ...
        'Color', [0.5 0.5 0.5], varargin{:}, 'XLimInclude', 'off');
    set(get(get(h(i), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end