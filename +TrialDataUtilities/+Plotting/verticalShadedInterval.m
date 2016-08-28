function h = verticalShadedInterval(mat, varargin)
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

yl = ylim();
yd = diff(yl);
yl(1) = yl(1) - 5*yd;
yl(2) = yl(2) + 5*yd;

X = [mat(:, 1)'; mat(:, 1)'; mat(:, 2)'; mat(:, 2)'];
Y = repmat([yl(1); yl(2); yl(2); yl(1)], 1, size(mat, 1));
    
h =  fill(X, Y, p.Results.Color, ...
    'YLimInclude', 'off', 'EdgeColor', 'none', ...
    'FaceAlpha', p.Results.alpha, ...
    'Parent', p.Results.Parent, p.Unmatched);

TrialDataUtilities.Plotting.hideInLegend(h);
end