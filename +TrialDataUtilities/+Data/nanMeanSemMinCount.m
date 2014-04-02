function [m, se, n] = nanMeanSemMinCount(x, dim, minCount)
% [mean, sem, count] = nanmeanMinCount(x,dim,minCount)
% Computes mean value along dimension dim, ignoring NaNs and 
% marking as Nan when the number of non-nan values is below minCount

if nargin == 2
    minCount = 1;
end

nans = isnan(x);
x(nans) = 0;

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
m = sum(x,dim) ./ nThresh;
se = std(x, [], dim) ./ sqrt(nThresh);

%m(tooFew) = NaN;
end

