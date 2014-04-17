function [xd, yd, zd] = getPointsToAxisDataScaling(axh, mode3d)
% multiply xd,yd,zd by a measurement in points to get data units

    if nargin < 1
        axh = gca;
    end
    if nargin < 2 || mode3d
        % 2-d mode
        %figh = getParentFigure(axh);
        axUnits = get(axh, 'Units');
        switch axUnits
            case 'inches'
                utop = 72;
                apos = get(axh, 'Position');
            case 'centimeters'
                utop = 72 / 2.54;
                apos = get(axh, 'Position');
            case 'points'
                utop = 1;
                apos = get(axh, 'Position');
            case {'normalized', 'characters', 'pixels'}
                set(axh, 'Units', 'points');
                utop = 1;
                apos = get(axh, 'Position');
                set(axh, 'Units', axUnits);
        end
               
        awPoints = apos(3) * utop;
        ahPoints = apos(4) * utop;
        
        xl = get(axh, 'XLim');
        yl = get(axh, 'YLim');
        awData = diff(xl);
        ahData = diff(yl);
        
        xd = awData / awPoints;
        yd = ahData / ahPoints;  
    end
end

function fig = getParentFigure(axh)
    % if the object is a figure or figure descendent, return the
    % figure. Otherwise return [].
    fig = axh;
    while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
      fig = get(fig,'parent');
    end
end
    