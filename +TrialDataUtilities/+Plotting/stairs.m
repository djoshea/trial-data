function h = stairs(x,y, varargin)
% like stairs, but extends the final X point one to the right
% stairs(x, y, ...)
%
% params:
%   color - width of stairs line, defaults to next color in ColorOrder
%   line width 


    p = inputParser();
    p.addParameter('color', [], @(x) true);
    p.addParameter('lineWidth', 1, @isscalar);
    p.addParameter('alpha', 1, @isscalar); % alpha for mean line, not supported yet
    p.addParameter('lastX', [], @isvector);
    p.addParameter('axh', gca, @ishandle);
    p.addParameter('clipping', 'on', @ischar);
    p.KeepUnmatched = true;
    p.CaseSensitive = false;
    p.parse(varargin{:});

    axh = p.Results.axh;
  
    if isempty(p.Results.lastX)
        x = makecol(x);
        lastX = x(end, :) + x(end, :) - x(end-1, :);
    else
        lastX = p.Results.lastX;
    end
    
    if isvector(x)
        x = makecol(x);
        y = makecol(y);
        assert(numel(x) == size(y, 1), 'X and Y must have same number of rows');
        assert(isscalar(lastX));
    else
        if isscalar(lastX);
            lastX = repmat(lastX, 1, size(x, 2));
        end
    end
    
    xWithLast = [x; lastX];
    yWithLast = [y; y(end)];
    
    if isempty(p.Results.color)
        colorArg = {};
    else
        colorArg = {'Color', p.Results.color};
    end
    
    h = stairs(xWithLast, yWithLast, colorArg{:}, 'LineWidth', p.Results.lineWidth, 'Parent', axh, 'Clipping', p.Results.clipping);
end
