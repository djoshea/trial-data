function ceval = evalColorMap2DAt(img, at1, at2, lim1, lim2, colorspace)
    % a 2D colormap is essentially a color image in RGB (or colorspace if
    % specified). intrpolation will be done in two dimensions to find the
    % right color. img should be n1 x n2 x 3

    assert(size(img, 3) == 3);
    
    if isscalar(at2)
        at2 = repmat(at2, numel(at1), 1);
    end
    if isscalar(at1)
        at1 = repmat(at1, numel(at2), 1);
    end
    at1 = makecol(at1);
    at2 = makecol(at2);
    
    assert(isvector(at1) && isvector(at2) && numel(at1) == numel(at2));
    
    if nargin < 4 || isempty(lim1)
        lim1 = [0 1];
    end
    if nargin < 5 || isempty(lim2)
        lim2 = [0 1];
    end
    
    N1 = size(img, 1);
    N2 = size(img, 2);
    
    x1 = linspace(lim1(1), lim1(2), N1)';
    x2 = linspace(lim2(1), lim2(2), N2)';
    
    gi =  {griddedInterpolant({x1, x2}, img(:, :, 1)); ...
           griddedInterpolant({x1, x2}, img(:, :, 2)); ...
           griddedInterpolant({x1, x2}, img(:, :, 3))};
    
    ceval = nan(numel(at1), 1);
    
    ceval(:, 1) = gi{1}(at1, at2);
    ceval(:, 2) = gi{2}(at1, at2);
    ceval(:, 3) = gi{3}(at1, at2);
    
    if nargin > 5 && ~strcmp(colorspace, 'RGB')
        ceval = TrialDataUtilities.Color.convert(sprintf('%s->RGB', colorspace), ceval);
    end
end