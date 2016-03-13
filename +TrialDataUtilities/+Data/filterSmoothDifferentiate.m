function out = filterSmoothDifferentiate(mat, samplingRateHz, varargin)
% out = filterSmoothDifferentiate(mat, samplingRateHz, varargin)
% uses savitzy golay to compute smooth derivative ignoring nan values at
% the start and end of each trace along dim.
% parameters:
%   dim - time dimension, defaults to first non singleton dim
%   smoothing - number of samples over which to smooth, must be odd,
%       default 7
%   order - which derivative? default = 1
%   polynomial order - for savitzy golay filter

    p = inputParser;
    p.addParameter('dim', [], @(x) isempty(x) || isscalar(x));
    % these refer to parameters of the savitzy golay filter
    p.addParameter('smoothing', 7, @(x) isscalar(x) && mod(x, 2) == 1);
    p.addParameter('order', 1, @isscalar);
    p.addParameter('polynomialOrder', 2, @isscalar);
    p.parse(varargin{:});
    
    if iscell(mat)
        out = cellfun(@filterMat, mat, 'UniformOutput', false);
    else
        out = filterMat(mat);
    end
    
    function mat = filterMat(mat)
        dim = p.Results.dim;
        if isempty(dim)
            % auto choose first non singular dimension
            dim = TensorUtils.firstNonSingletonDim(mat);
        end
        
        if dim ~= 1
            % place dim at second dimension slot
            ndimsOrig = ndims(mat);
            mat = TensorUtils.shiftdimToFirstDim(mat, dim);
        end
        
        nCol = size(mat(:, :), 2);
        
        w = -samplingRateHz ^ p.Results.order;

        % for each column, filter from the first to last non-NaN ind
        emptyCol = nan(size(mat, 1), 1);
        for iC = 1:nCol
            thisCol = mat(:, iC);
            filtCol = emptyCol;
            colMask = ~isnan(thisCol);
            start = find(colMask, 1, 'first');
            stop = find(colMask, 1, 'last');
            
            filtCol(start:stop) = w .* TrialDataUtilities.Data.savitzkyGolayFilt(thisCol(start:stop), ...
                p.Results.polynomialOrder, p.Results.order, p.Results.smoothing);
            
            mat(:, iC) = filtCol;
        end

        if dim ~= 1
            % place dim at second dimension slot
            mat = TensorUtils.unshiftdimToFirstDim(mat, dim, ndimsOrig);
        end
    end
end