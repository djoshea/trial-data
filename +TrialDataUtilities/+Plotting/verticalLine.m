function h = verticalLine(xVal, varargin)

p = inputParser();
p.addParameter('axh', gca, @ishandle);
p.KeepUnmatched = true;
p.parse(varargin{:});

h = gobjects(numel(xVal), 1);
yl = ylim();
yd = diff(yl);
yl(1) = yl(1) - 5*yd;
yl(2) = yl(2) + 5*yd;

for i = 1:numel(xVal)
    h(i) = plot([xVal(i), xVal(i)], yl, ...
        'Color', [0.5 0.5 0.5], 'YLimInclude', 'off', ...
        'Parent', p.Results.axh, p.Unmatched);
    set(get(get(h(i), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end