function h = patchcircle(xc, yc, diameter, varargin)

    p = inputParser;
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('z', 0, @isscalar);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    axh = p.Results.axh;

    assert(numel(xc) == numel(yc), 'Inputs must have same length');
    xc = makerow(xc);
    yc = makerow(yc);
    
    % generate patch coordinates for circles
    t = linspace(0, 2*pi, 20)';
    % these must be column vectors
    xv = diameter/2 * cos(t);
    yv = diameter/2 * sin(t);
    
    % scale diameter in each dimension from points to data units
    [xd, yd] = TrialDataUtilities.Plotting.getPointsToAxisDataScaling(axh);
    
    % these will be nVals x nPts
    X = bsxfun(@plus, xc, xv * xd);
    Y = bsxfun(@plus, yc, yv * yd);
    
    if p.Results.z ~= 0
        Z = p.Results.z * ones(size(X));
        h = patch(X, Y, Z, 'k', 'EdgeColor', 'none');
    else
        h = patch(X, Y, 'k', 'EdgeColor', 'none');
    end
    uistack(h, 'top');
    
    flds = fieldnames(p.Unmatched);
    for iF = 1:numel(flds)
        set(h, flds{iF}, p.Unmatched.(flds{iF}));
    end
end