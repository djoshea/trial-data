function gmed = geometric_median(X, median_over_dim, vector_space_dim)
% gm = geometric_median(X, median_over_dim, vector_space_dim)
% Calculate geometric median over axis median_over_dim. Each point is a point in size(X, vector_space_dim) space.
% Each dimension not included in median_over_dim or vector_space_dim will be looped over.

if nargin < 2
    median_over_dim = 1;
end
if nargin < 3
    vector_space_dim = 2;
end

ndimX = ndims(X);
do_permute = median_over_dim ~= 1 || vector_space_dim ~= 2;
if do_permute
    X = TensorUtils.shiftdimToFirstDim(X, [vector_space_dim, ]);
end
sz = size(X);

gmed = nan([sz(1) 1 sz(3:end)], like=X);
nPages = prod(sz(3:end));
for iP = 1:nPages
    gmed(:, :, iP) = gmed_inner(X(:, :, iP));
end

if do_permute
    gmed = TensorUtils.unshiftdimToFirstDim(gmed, [vector_space_dim, median_over_dim], ndimX);
end

end

function gmed = gmed_inner(X)
    % X is dims x num_observations
    options = optimset('Display','off');
    vecnorm_omitnan = @(A, p, dim) sum(abs(A).^p, dim, 'omitnan').^(1./p);
    gmed = fminsearch(@(y) sum(vecnorm_omitnan(X-y, 2, 2), 1, 'omitnan'), mean(X, 2, 'omitnan'), options);
end
