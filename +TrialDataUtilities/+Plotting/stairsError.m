function [hMean, hError] = stairsError(x,y, varargin)
% like stairs, but with error regions.
% stairsError(x, y, yerror, ...) - plot +/- yerror 
% stairsError(x, y, ylo, yhi, ...) - plot confidence intervals
%
% for stairsError, the last value of x can be specified manually using
% parameter: 'lastX', otherwise it will be picked to match x(end)+(x(end)-x(end-1))
%
%
% params:
%   color - width of stairs line, defaults to next color in ColorOrder
%   line width 
%   errorColor - defaults to lighter version of color
%   errorAlpha - defaults to 0.8

    p = inputParser();
    p.addRequired('err1', @(x) isnumeric(x) && ismatrix(x));
    p.addOptional('err2', [], @(x) isnumeric(x) && ismatrix(x));
    p.addParameter('color', [], @(x) true);
    p.addParameter('errorStyle', 'fill', @ischar); % fill or stairs
    p.addParameter('errorColor', [], @(x) true);
    p.addParameter('errorAlpha', 0.5, @isscalar);
    p.addParameter('lineWidth', 1, @isscalar);
    p.addParameter('alpha', 1, @isscalar); % alpha for mean line, not supported yet
    p.addParameter('lastX', [], @isvector);
    p.addParameter('axh', gca, @ishandle);
    p.KeepUnmatched = true;
    p.CaseSensitive = false;
    p.parse(varargin{:});

    axh = p.Results.axh;
    
    if isempty(p.Results.err2)
        yhi = y + p.Results.err1;
        ylo = y - p.Results.err1;
    else
        ylo = p.Results.err1;
        yhi = p.Results.err2; 
    end
    
%     function tf = szeq(varargin)
%         szc = cellfun(@size, varargin{:}, 'UniformOutput', false);
%         tf = isequal(szc{:});
%     end

    if isvector(y)
        y = makecol(y);
        ylo = makecol(ylo);
        yhi = makecol(yhi);
    end
    
    if isempty(p.Results.lastX)
        x = makecol(x);
        lastX = x(end, :) + x(end, :) - x(end-1, :);
    else
        lastX = p.Results.lastX;
    end
    
%     assert(szeq(y,ylo,yhi), 'All inputs must have the same size and shape');
    if isvector(x)
        x = makecol(x);
        assert(numel(x) == size(y, 1), 'X and Y must have same number of rows');
        assert(isscalar(lastX));
    else
        if isscalar(lastX);
            lastX = repmat(lastX, 1, size(x, 2));
        end
    end
    
    if isempty(p.Results.color)
        colorArg = {};
    else
        colorArg = {'Color', p.Results.color};
    end
        
    % build polygons for fill 
    npts = size(x, 1);
    nlin = size(x, 2);
    % X Y
    % 1 H
    % 2 H
    % 2 H
    % 3 H
    % 3 L
    % 2 L
    % 2 L
    % 1 L
    
    % add on a tail value so that the last point gets plotted
    x = cat(1, x, lastX);
    y = cat(1, y, y(end, :));
    ylo = cat(1, ylo, ylo(end, :));
    yhi = cat(1, yhi, yhi(end, :));

    holdstate = ishold;
    hMean = stairs(x, y, colorArg{:}, 'LineWidth', p.Results.lineWidth);
    hold on;

    switch p.Results.errorStyle
        case 'fill'
            X = repelem(x, 2, 1);
            Xud = repelem(flipud(x(2:end-1, :)), 2, 1);
            X = circshift(cat(1, X, Xud), -1, 1);
            Y = cat(1, repelem(yhi(1:end-1, :), 2, 1), repelem(flipud(ylo(1:end-1)), 2, 1));

            hError = patch(X, Y, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', p.Results.errorAlpha, 'Parent', axh);
            for i = 1:nlin
                if isempty(p.Results.errorColor)
                    hError(i).FaceColor = hMean(i).Color;
                else
                    hError(i).FaceColor = p.Results.errorColor;
                end
            end
            
        case 'stairs'
            hError = stairs(x, ylo, 'LineWidth', p.Results.lineWidth/2, 'Parent', axh);
            hold on;
            hError(:, 2) = stairs(x, yhi, 'LineWidth', p.Results.lineWidth/2, 'Parent', axh);
            
            for i = 1:nlin
                if isempty(p.Results.errorColor)
                    set(hError(i, :), 'Color', hMean(i).Color);
                else
                    set(hError(i, :), 'Color', p.Results.errorColor);
                end
            end
    end
 
    if ~holdstate
        hold off;
    end
    uistack(hMean, 'top');
    TrialDataUtilities.Plotting.hideInLegend(hError);
end
