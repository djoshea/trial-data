function [ymean, ysem, xvals] = averageMultipleTimeseriesByTimeseries(xinCell, yinCell, varargin)
% given 2 timeseries vectors, x and y, resample y as a function of x, i.e.
% average the values of y based on bins of x. Provide 'n' (default 1000) to
% linear space bins over the range of x, or provide 'xvals' vector to
% specify the values directly

p = inputParser();
p.addParamValue('n', 1000, @isscalar);
p.addParamValue('xvals', [], @isvector);
p.parse(varargin{:});

% build a concatenated list of x values, y values, and which-trial indices
xcat = cat(1, xinCell{:});
[ycat, trialIdx] = TensorUtils.catWhich(1, yinCell{:});

% figure out xvals bins
if ~isempty(p.Results.xvals)
    xvals = p.Results.xvals;
else
    % compute linear space bins over the full range
    n = p.Results.n;
    xmin = nanmin(xcat);
    xmax = nanmax(xcat);
    xvals = linspace(xmin, xmax, n);
end

nTrials = numel(xinCell);
nBins = numel(xvals);

% compute histogram over x values from all trials
[freq, xbin] = histc(xcat, xvals);

% build up a sparse data matrix that is 
% nBins x nTrials
dataByTrial = accumarray([xbin, trialIdx], ycat, [nBins, nTrials], ...
    @nanmean, NaN, false);

% lineraly interpolate to fill in missing bin values on each trial 
for iTrial = 1:nTrials
    d = dataByTrial(:, iTrial);
    mask = ~isnan(d);
    if nnz(mask) <= 2
        continue;
    end
    dataByTrial(:, iTrial) = interp1(xvals(mask), d(mask), xvals, 'linear');
end

ymean = nanmean(dataByTrial, 2)';
ysem = nansem(dataByTrial, 2)';


