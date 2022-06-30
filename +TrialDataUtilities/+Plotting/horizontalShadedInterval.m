function h = horizontalShadedInterval(mat, varargin)
% mat should be nIntervals by 2 

if isvector(mat) && numel(mat) == 2
    mat = makerow(mat);
end
assert(size(mat, 2) == 2);

p = inputParser();
p.addParameter('Parent', gca, @ishandle);
p.addParameter('Color', [0.3 0.3 0.3], @(x) true);
p.addParameter('alpha', 0.5, @isscalar);
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.parse(varargin{:});

xl = xlim();
xd = diff(xl);
xl(1) = xl(1) - 5*xd;
xl(2) = xl(2) + 5*xd;

X = repmat([xl(1); xl(2); xl(2); xl(1)], 1, size(mat, 1)); % fill width
Y = [mat(:, 1)'; mat(:, 1)'; mat(:, 2)'; mat(:, 2)'];

h =  fill(X, Y, p.Results.Color, ...
    'XLimInclude', 'off', 'EdgeColor', 'none', ...
    'FaceAlpha', p.Results.alpha, ...
    'Parent', p.Results.Parent, p.Unmatched);

TrialDataUtilities.Plotting.hideInLegend(h);
end