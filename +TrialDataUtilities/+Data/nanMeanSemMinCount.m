function [m, se, n, stdev] = nanMeanSemMinCount(x, dim, minCount)
% [mean, sem, count, stdev] = nanmeanMinCount(x,dim,minCount)
% Computes mean value along dimension dim, ignoring NaNs and 
% marking as Nan when the number of non-nan values is below minCount. Also
% computes standard deviation.

if nargin == 2
    minCount = 1;
end

nans = isnan(x);

% auto-choose dim as first non singleton dimension
if ~exist('dim', 'var') || isempty(dim)
    if numel(x) == 1
        dim = 1;
    else
        dim = find(size(x) > 1, 1, 'first');
    end
end
    
% count up non-NaNs along dim
n = sum(~nans,dim);
tooFew = n < minCount;

% automatically inserts NaNs where too few trials
nThresh = n;
nThresh(tooFew) = NaN;

% Sum up non-NaNs, and divide by the number of non-NaNs.
x0 = x; x0(nans) = 0;
m = sum(x0,dim) ./ nThresh;

% compute standard deviation 
stdev = nanstd(x, [], dim);
se = stdev ./ sqrt(nThresh);

% m(tooFew) = NaN; % not necessary since nThresh will be NaN here
stdev(tooFew) = NaN;
se(tooFew) = NaN;
end

