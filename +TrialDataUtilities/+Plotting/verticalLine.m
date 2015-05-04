function h = verticalLine(xVal, varargin)

p = inputParser();
p.addParameter('axh', gca, @ishandle);
p.KeepUnmatched = true;
p.parse(varargin{:});

h = nan(numel(xVal), 1);
if strcmp(get(p.Results.axh, 'YScale'), 'log')
    yv = [1e-100 1e100];
else
    yv = [-1e10 1e10];
end
for i = 1:numel(xVal)
    h(i) = plot([xVal(i), xVal(i)], yv, ...
        'Color', [0.5 0.5 0.5], 'YLimInclude', 'off', ...
        'Parent', p.Results.axh, p.Unmatched);
    set(get(get(h(i), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end