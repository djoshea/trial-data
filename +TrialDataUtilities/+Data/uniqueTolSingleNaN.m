function [c, ia, ic] = uniqueTolSingleNaN(a, varargin)

[c, ia, ic] = uniquetol(a, varargin{:});

if istable(a) || (ismatrix(c) && size(c, 2) > 1)
    % table mode or 'rows' mode
    % here we replace the missing value in each column with an unused non-missing value and rerun unique
    % to see which rows we need to keep
    
    cc = c;
    has_nan = false;
    for col = 1:size(c, 2)
        if istable(a)
            col_vec = cc{:, col};
        else
            col_vec = cc(:, col);
        end
        mask_nan = ismissing(col_vec);
        
        if ~any(mask_nan)
            continue;
        end
        
        has_nan = true;
        
        if isstring(col_vec)
            cc{mask_nan, :} = matlab.lang.makeUniqueStrings("missing", col_vec(~mask_nan));
        elseif isnumeric(col_vec)
            val = Inf;
            while ismember(val, col_vec(~mask_nan))
                val = randn();
            end
            cc{mask_nan, :} = val;
        end
    end
    
    if has_nan
        [~, ~, icc] = unique(cc);
        drop_inds = setdiff(1:size(cc, 1), icc);
        ic = icc;
    else
        drop_inds = [];
    end
    
    c(drop_inds, :) = [];
    ia(drop_inds) = []; 
    
else
    % simple numeric vector mode
    mask_nan = ismissing(c);
    nan_inds = find(mask_nan);
    if numel(nan_inds) > 1
        drop_inds = nan_inds(2:end);
    else
        drop_inds = [];
    end

    if ~isempty(drop_inds)
        c(drop_inds, :) = [];
        ia(drop_inds) = [];
    
        mask_reassign = ismember(ic, drop_inds);
        ic(mask_reassign) = nan_inds(1);
    end

end



end