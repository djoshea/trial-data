function h = patchsphere(xc, yc, zc, diameter, varargin)
% draw 3d spheres at locations in xc, yc, zc

    p = inputParser;
    p.addParamValue('axh', gca, @ishandle);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    axh = p.Results.axh;

    assert(numel(xc) == numel(yc), 'Inputs must have same length');
    assert(numel(xc) == numel(zc), 'Inputs must have same length');
   
    xc = makerow(xc);
    yc = makerow(yc);
    zc = makerow(zc);
    
    [xx, yy, zz] = sphere(10);
    
    % scale diameter in each dimension from points to data units
    [xd, yd, zd] = TrialDataUtilities.Plotting.getPointsToAxisDataScaling(axh, true);
    
    % these will be nVals x nPts
    h = nan(numel(xc), 1);
    for i = 1:numel(xc) 
        X = xx*xd*diameter + xc(i);
        Y = yy*yd*diameter + yc(i);
        Z = zz*zd*diameter + zc(i);
        h(i) = surf(axh, X, Y, Z, 'EdgeColor', 'none');
    end
    
    flds = fieldnames(p.Unmatched);
    for iF = 1:numel(flds)
        set(h, flds{iF}, p.Unmatched.(flds{iF}));
    end
end