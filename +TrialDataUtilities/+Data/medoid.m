function [m, idx] = medoid(X, median_over_dim, vector_space_dim)
% 
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
X = TensorUtils.shiftdimToFirstDim(X, [median_over_dim, vector_space_dim]);
sz = size(X);

m = nan([1 sz(2:end)], like=X);
idx = nan([1 1 sz(3:end)], like=X);

for iP = 1:nPages
    [m(:, :, iP), idx(:, :, iP)] = medoid_inner(X(:, :, iP));
end

m = TensorUtils.unshiftdimToFirstDim(m, [median_over_dim, vector_space_dim], ndimX);
idx = TensorUtils.unshiftdimToFirstDim(idx, [median_over_dim, vector_space_dim], ndimX);

end

function [m, idx] = medoid_inner(X)
% 
% Calculate medoid
%
% Input:
%   X: d x n data matrix
%
% Output:
%   m: medoid
%   idx: index of medoids
%
% Written by Detang Zhong (detang.zhong@canada.ca). Modified by Dan O'Shea
%
% Demo:
%  X = [12,  9, 61, 76,  2, 17, 12, 11, 26,  0;
%        65, 72,  7, 64, 21, 92, 51, 48,  9, 65;
%        39,  7, 50, 56, 29, 79, 47, 45, 10, 52;
%        70, 12, 23, 97, 86, 14, 42, 90, 15, 16;
%        13,  7,  2, 47, 80, 53, 23, 59,  7, 15;
%        83,  2, 40, 12, 22, 75, 69, 61, 28, 53]
% 
% X =
% 
%     12     9    61    76     2    17    12    11    26     0
%     65    72     7    64    21    92    51    48     9    65
%     39     7    50    56    29    79    47    45    10    52
%     70    12    23    97    86    14    42    90    15    16
%     13     7     2    47    80    53    23    59     7    15
%     83     2    40    12    22    75    69    61    28    53
% 
% m = medoid(X)
% 
% m =
% 
%     12
%     51
%     47
%     42
%     23
%     69
%
%% Determine the size of X
[d,n] = size(X);
v = dot(X,X,1);
%% Calculate Euclidean distance matrix
% D = v+v'-2*(X'*X);                  
D = -2*(X'*X); 
D = bsxfun(@plus, D, v); 
D = bsxfun(@plus, D, v');
%% Reduce chance of numerical problems
D(sub2ind([n,n],1:n,1:n)) = 0; 
%% Get the medoid and its index
md = mean(D,2);
[~,idx] = min(md);
m = X(:,idx);
end
