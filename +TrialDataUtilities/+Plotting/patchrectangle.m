function h = patchrectangle(x1, y1, x2, y2, varargin)

    p = inputParser;
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('z', 0, @isscalar);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    axh = p.Results.axh;

    x1 = makerow(x1);
    x2 = makerow(x2);
    y1 = makerow(y1);
    y2 = makerow(y2);
    
    X = [x1; x1; x2; x2];
    Y = [y1; y2; y2; y1];

    valid = ~all(isnan(X) | isnan(Y), 1);
    if ~any(valid), 
        h = NaN;
        return;
    end
    X = X(:, valid);
    Y = Y(:, valid);
    
    if p.Results.z ~= 0
        Z = p.Results.z * ones(size(X));
        h = patch(X, Y, Z, 'k', 'EdgeColor', 'none');
    else
        h = patch(X, Y, 'k', 'EdgeColor', 'none');
    end
    
    flds = fieldnames(p.Unmatched);
    for iF = 1:numel(flds)
        set(h, flds{iF}, p.Unmatched.(flds{iF}));
    end
end
