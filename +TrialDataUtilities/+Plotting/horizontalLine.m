function h = horizontalLine(yVal, varargin)

p = inputParser();
p.addParameter('xlim', [NaN NaN], @isvector);
p.addParameter('Parent', gca);
p.KeepUnmatched = true;
p.parse(varargin{:});

xl = p.Results.xlim;

% auto xlim where xl is nan
xlim_auto = xlim(p.Results.Parent);
xd = diff(xlim_auto);
xv = xlim_auto;
xv(1) = xv(1) - 5*xd;
xv(2) = xv(2) + 5*xd;

if ~isnan(xl(1))
    xv(1) = xl(1);
end
if ~isnan(xl(2))
    xv(2) = xl(2);
end

h = gobjects(numel(yVal), 1);
for i = 1:numel(yVal)
    h(i) = plot(xv, [yVal(i), yVal(i)], ...
        'Color', [0.5 0.5 0.5], p.Unmatched, 'Parent', p.Results.Parent, 'XLimInclude', 'off');
    set(get(get(h(i), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end