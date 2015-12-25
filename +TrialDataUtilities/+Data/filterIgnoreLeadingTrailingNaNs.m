function mat = filterIgnoreLeadingTrailingNaNs(B, A, mat, varargin)
    % matFilt = filterIgnoreLeadingTrailingNaNs(B, A, mat)
    
    p = inputParser;
    % if true subtracts first sample from each signal which can help reduce onset transients for IIR filters
    p.addParameter('subtractFirstSample', false, @islogical);
    % if true, the first sample will be added back on afterwards
    p.addParameter('addBackFirstSample', false, @islogical);
    p.parse(varargin{:});
    
    if ~any(isnan(mat(:)))
        firstSample = mat(:, 1);
        if p.Results.subtractFirstSample
            mat = mat - firstSample;
        end
        mat = filter(mat, B, A);
        if p.Results.subtractFirstSample && p.Results.addBackFirstSample
            mat = mat + firstSample;
        end
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
            if p.Results.subtractFirstSample && ~isempty(start)
                firstSample = thisCol(start);
                thisCol = thisCol - firstSample;
            end
            filtCol(start:stop) = filter(B, A, thisCol(start:stop));
            if p.Results.subtractFirstSample && p.Results.addBackFirstSample && ~isempty(start)
                filtCol = filtCol + firstSample;
            end
            mat(:, iC) = filtCol;
        end
    end
end