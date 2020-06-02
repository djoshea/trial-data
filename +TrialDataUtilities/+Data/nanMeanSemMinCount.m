function [m, se, n, stdev] = nanMeanSemMinCount(x, dim, minCount, minRatio, varargin)
% [mean, sem, count, stdev] = nanmeanMinCount(x,dim,minCount)
% Computes mean value along dimension dim, ignoring NaNs and 
% marking as Nan when the number of non-nan values is below minCount. Also
% computes standard deviation.

if ~exist('minCount', 'var') || isempty(minCount)
    minCount = 1;
end

p = inputParser();
p.addParameter('weights', [], @isnumeric);
p.addParameter('assumePoissonStatistics', false, @islogical);
p.addParameter('poissonCountMultipliers', [], @isnumeric);
p.parse(varargin{:});

assumePoissonStatistics = p.Results.assumePoissonStatistics;

w = p.Results.weights;
if isempty(w)
    w = ones(size(x));
end
if isvector(w) 
    if size(w, dim) ~= size(x, dim)
        w = TensorUtils.orientVectorAlongDim(w, dim);
    end
    w = w .* ones(size(x));
end
assert(isequal(size(x), size(w)));

if assumePoissonStatistics
    % make poissonCountMultipliers same size as x
    poissonCountMultipliers = p.Results.poissonCountMultipliers;
    if isempty(poissonCountMultipliers)
        error('poissonCountMultipliers required with assumePoissonStatistics');
    end
    if isvector(poissonCountMultipliers) 
        if size(poissonCountMultipliers, dim) ~= size(x, dim)
            poissonCountMultipliers = TensorUtils.orientVectorAlongDim(poissonCountMultipliers, dim);
        end
        poissonCountMultipliers = poissonCountMultipliers .* ones(size(x));
    end
end

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

if exist('minRatio', 'var') && ~isempty(minRatio) && ~isnan(minRatio) && minRatio > 0
    minCountFromRatio = ceil(minRatio * size(x, dim));
    tooFew = tooFew | n < minCountFromRatio;
end

% automatically inserts NaNs where too few trials
nThresh = n;
nThresh(tooFew) = NaN;

% Sum up non-NaNs, and divide by the number of non-NaNs.
w_sum = sum(w, dim);
m = sum(x .* w, dim) ./ w_sum;

% compute variance
if assumePoissonStatistics
    % assume that variance == mean, before any multiplication is performed. so if mean is (mult * count), var is mult^2 * count  = x * mult
    % with weights, we have, mean = wnorm1 * mult1 * count1 + wnorm1 * mult2 * count2 + ..
    % then var = wnorm1^2 * mult1^2 * count1 + wnorm2^2 * mult2^2 * count2 =  wnorm1^2 * mult1 * x1 + wnorm2^2 * mult2 * x2
    v = sum(w.* poissonCountMultipliers .* x, dim) ./ (sum(w, dim) - 1);
else
    % compute frequency-weighted sample variance
    v = TrialDataUtilities.Stats.varWeightedCorrected(x, w, dim, true);
end

% compute standard deviation 
stdev = sqrt(v);
se = stdev ./ sqrt(nThresh);

m(tooFew) = NaN;
stdev(tooFew) = NaN;
se(tooFew) = NaN;

end

