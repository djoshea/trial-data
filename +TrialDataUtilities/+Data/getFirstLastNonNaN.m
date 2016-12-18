function [idxFirst, idxLast] = getFirstLastNonNaN(data, dim)
% for tensor data, return the min and max time for which
% data is present on each trial. If data is a a matri

    % move the relevant dimension to dim 1
    [data, ndimOrig] = TensorUtils.shiftdimToFirstDim(data, dim);

    sz = size(data);
    sz(1) = 1;
    [idxFirst, idxLast] = deal(nan(sz));

    N = prod(sz);

    for i = 1:N
        mask = ~isnan(data(:, i));
        if any(mask)
            idxFirst(i) = find(mask, 1, 'first');
            idxLast(i) = find(mask, 1, 'last');
        end
    end

    % move the dimension back to its original place
    idxFirst = TensorUtils.unshiftdimToFirstDim(idxFirst, dim, ndimOrig);
    idxLast = TensorUtils.unshiftdimToFirstDim(idxLast, dim, ndimOrig);
    
end