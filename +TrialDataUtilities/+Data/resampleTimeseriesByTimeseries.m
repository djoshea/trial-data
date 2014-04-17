function [ymean, xvals] = resampleTimeseriesByTimeseries(xin, yin, varargin)
% given 2 timeseries vectors, x and y, resample y as a function of x, i.e.
% average the values of y based on bins of x. Provide 'n' (default 1000) to
% linear space bins over the range of x, or provide 'xvals' vector to
% specify the values directly

p = inputParser();
p.addParamValue('n', 1000, @isscalar);
p.addParamValue('xvals', [], @isvector);
p.parse(varargin{:});

if ~isempty(p.Results.xvals)
    xvals = p.Results.xvals;
else
    n = p.Results.n;
    xvals = linspace(nanmin(xin), nanmax(xin), n);
end

[~, xbin] = histc(xin, xvals);
% numel(x) x numel(unique(x))
ysp = sparse(1:numel(xin), xbin, yin);
ymean = full(sum(ysp)./sum(ysp~=0)); % each column in the matrix corresponds to a "index"
%ysem = full(std(ysp)./sqrt(sum(ys

ysem = [];