function h = identityLine(varargin)

p = inputParser();
p.addParameter('slope', 1, @isscalar);
p.addParameter('offset', 0, @isscalar);
p.addParameter('Parent', gca);
p.addParameter('xlim', [NaN NaN], @isvector);
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
yv = p.Results.slope*xv + p.Results.offset;

h = plot(xv, yv, ...
    'Color', [0.5 0.5 0.5], p.Unmatched, 'Parent', p.Results.Parent, 'XLimInclude', 'off', 'YLimInclude', 'off');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

end
