function [hLines, traceCenters] = plotStackedTraces(tvec, mat, varargin)
% [hLines, traceCenters] = plotStacked(tvec, mat, varargin)
% Plots the columns of mat stacked vertically on top of each other
% options include all options for plot(...) plus:
%   normalize: [false] scale each signal to the same range before plotting
%   spacingFraction: [1.02] space each trace by this fraction of the previous
%     trace's range

p = inputParser();
p.addParamValue('normalize', false, @islogical);
p.addParamValue('intercalate', true, @islogical);
p.addParamValue('spacingFraction', 1.05, @isscalar);
p.addParamValue('colorByTrace', [], @ismatrix);
p.KeepUnmatched = true;
p.parse(varargin{:});

% subtract the min so each trace starts at zero
mat = fliplr(mat);
matShift = bsxfun(@minus, mat, nanmin(mat, [], 1));

if p.Results.normalize
    matShift = bsxfun(@rdivide, matShift, range(matShift, 1));
end

ranges = range(matShift, 1);

% figure out where each trace should start
if p.Results.intercalate
    % the offsets should be set so that pairs of points across
    % all time points should satisfy:
    %   trace1(t) * spacingFraction <= trace(2) + offset
    % or:
    %   offset = max_t (trace1(t) * spacingFraction - trace2(t)
    
    deltas = matShift(:, 1:end-1) * p.Results.spacingFraction - matShift(:, 2:end);
    maxDeltas = nanmax(deltas, [], 1);
    traceOffsets = [0, cumsum(maxDeltas)];
    
else 
    rangesPadded = ranges * (p.Results.spacingFraction);
    cs = cumsum(rangesPadded);
    traceOffsets = [0, cs(1:end-1)];
end

traceCenters = (traceOffsets + ranges / 2)';

matShift = bsxfun(@plus, matShift, traceOffsets);

if isempty(p.Results.colorByTrace)
    hLines = plot(tvec, matShift, '-', 'Color', 'k', p.Unmatched);
else
    nTraces = size(matShift, 2);
    hLines = nanvec(nTraces);
    hold on;
    for i = 1:nTraces
        hLines(i) = plot(tvec, matShift(:, i), '-', 'Color', p.Results.colorByTrace(i, :), p.Unmatched);
    end
end
ylim([0 traceOffsets(end) + ranges(end)]);

end