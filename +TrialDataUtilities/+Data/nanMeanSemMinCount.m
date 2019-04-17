function [m, se, n, stdev] = nanMeanSemMinCount(x, dim, minCount, minRatio, w)
% [mean, sem, count, stdev] = nanmeanMinCount(x,dim,minCount)
% Computes mean value along dimension dim, ignoring NaNs and 
% marking as Nan when the number of non-nan values is below minCount. Also
% computes standard deviation.

if nargin < 3
    minCount = 1;
end
if nargin < 5
    w = ones(size(x));
end
if isvector(w)
    w = TensorUtils.orientVectorAlongDim(w, dim);
    w = w .* ones(size(x));
end
assert(isequal(size(x), size(w)));

invalid = isnan(x) | isnan(w) | w == 0;
x(invalid) = 0;
w(invalid) = 0;

% auto-choose dim as first non singleton dimension
if ~exist('dim', 'var') || isempty(dim)
    if numel(x) == 1
        dim = 1;
    else
        dim = find(size(x) > 1, 1, 'first');
    end
end

% count up non-NaNs along dim
n = sum(~invalid,dim);
tooFew = n < minCount;

if exist('minRatio', 'var') && ~isnan(minRatio) && minRatio > 0
    minCountFromRatio = ceil(minRatio * size(x, dim));
    tooFew = tooFew | n < minCountFromRatio;
end

% automatically inserts NaNs where too few trials
nThresh = n;
nThresh(tooFew) = NaN;

% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(x .* w, dim) ./ sum(w, dim);

% compute variance
v = TrialDataUtilities.Stats.varWeightedCorrected(x, w, dim);

% compute standard deviation 
stdev = sqrt(v);
se = stdev ./ sqrt(nThresh);

m(tooFew) = NaN;
stdev(tooFew) = NaN;
se(tooFew) = NaN;
end

