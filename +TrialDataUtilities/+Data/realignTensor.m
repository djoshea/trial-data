function [data_aligned, tvec_aligned] = realignTensor(data, alignInds, varargin)

% size(delays) should match size(data) but be scalar along dims where all elements along that dimension share the same delay
%
% parameters:
%   fillMode: can be a scalar value like NaN or 0 to fill the edges of the
%     signals as they shift left or right. or can be 'hold' to hold the
%     edge values of the signals out to the end

p = inputParser();
p.addParameter('dt', 1, @isscalar); % for generating aligned time vector 
p.addParameter('dim', 2, @isscalar);
p.addParameter('keepValidOnly', false, @islogical); % slice out the invalid tails
p.addParameter('fillMode', NaN, @(x) isscalar(x) || ischar(x));
p.addParameter('progress', true, @islogical);
p.addParameter('message', 'Aligning tensor slices', @isstringlike);
p.parse(varargin{:});

% expand to size of data except along time dim
dim = p.Results.dim;
szData = size(data);
szAlign = TensorUtils.sizeNDims(alignInds, numel(szData));

% check size(align)
szAlignReq = szData;
szAlignReq(dim) = 1;
assert(all(szAlignReq == szAlign | szAlign == 1), 'Size of alignInds is not correct for data');

constantAlignDims = setdiff(find(szAlign == 1), dim); % dims where the delay is constant across all elements (e.g. channels which shift together)
otherDims = TensorUtils.otherDims(szData, union(dim, constantAlignDims));

hold = false;
if ischar(p.Results.fillMode)
    if strcmp(p.Results.fillMode, 'hold')
        hold = true;
        fill = NaN;
    else
        error('Unknown fillMode');
    end
else
    fill = p.Results.fillMode;
end

if p.Results.progress
    prog = ProgressBar('R', p.Results.message);
else
    prog.update = @(varargin) true;
    prog.finish = @(varargin) true;
end

% reshape to put time dim first (T), loop-over-dims second (R), other "channel" dims 3rd (C). 
% data is now T x R x C, alignInds is 1 x R
data = TensorUtils.reshapeByConcatenatingDims(data, {dim, otherDims, constantAlignDims});
alignInds = TensorUtils.reshapeByConcatenatingDims(alignInds, {dim, otherDims, constantAlignDims});
[T, R, C] = size(data);

% what are the existing edges at in aligned time
left_align_idx = 1 - alignInds;
right_align_idx = T - alignInds;

if p.Results.keepValidOnly
    new_left_idx = max(left_align_idx, [], 'all');
    new_right_idx = min(right_align_idx, [], 'all');
else
    new_left_idx = min(left_align_idx, [], 'all');
    new_right_idx = max(right_align_idx, [], 'all');
end

inds_aligned = new_left_idx : new_right_idx;
ind_new_zero = 1 - new_left_idx;
Tnew = numel(inds_aligned);
data_aligned = nan(Tnew, R, C, 'like', data);
data_aligned(:) = fill;

for r = 1:R
    this_aligned_inds = round((1:T) - alignInds(r));
    idxGrab = this_aligned_inds >= new_left_idx & this_aligned_inds <= new_right_idx;
    idxInsert = this_aligned_inds(idxGrab) + ind_new_zero;
    %mask = idxInsert > 1 & idxInsert < T;
    %idxInsert = idxInsert(mask);
    %idxGrab = idxGrab(mask);
    slice = data(:, r, :);
    data_aligned(idxInsert, r, :) = slice(idxGrab, :, :);
    
    if hold
        data_aligned(1:idxInsert(1)-1, r, :) = repmat(data_aligned(idxInsert(1), r, :), [numel(1:idxInsert(1)-1) 1 1]);
        data_aligned(idxInsert(end)+1:T, r, :) = repmat(data_aligned(idxInsert(end), r, :), [numel(idxInsert(end)+1:T) 1 1]);
    end
    prog.update(r);
end
prog.finish();

szDataAligned = szData;
szDataAligned(dim) = Tnew;
data_aligned = TensorUtils.undoReshapeByConcatenatingDims(data_aligned, {dim, otherDims, constantAlignDims}, szDataAligned);
tvec_aligned = inds_aligned * p.Results.dt;

end