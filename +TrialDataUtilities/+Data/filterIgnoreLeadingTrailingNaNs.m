function mat = filterIgnoreLeadingTrailingNaNs(B, A, mat)
    % matFilt = filterIgnoreLeadingTrailingNaNs(B, A, mat)
    
    if ~any(isnan(mat(:)))
        mat = filter(mat, B, A);
    else
        % for each column, filter from the first to last non-NaN ind
        nCol = size(mat(:, :), 2);
        emptyCol = zeros(size(mat, 1), 1);
        for iC = 1:nCol
            thisCol = mat(:, iC);
            filtCol = emptyCol;
            colMask = ~isnan(thisCol);
            start = find(colMask, 1, 'first');
            stop = find(colMask, 1, 'last');
            
            filtCol(start:stop) = filter(B, A, thisCol(start:stop));
            mat(:, iC) = filtCol;
        end
    end
end