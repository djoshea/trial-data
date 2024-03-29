function [hl, hs] = errorshade(x, ym, ye, color, varargin)
% [hl, hs] = errorshade(x, ym, ye, color, varargin)

    p = inputParser();
    p.addParameter('axh', [], @(x) true);
    p.addParameter('showLine', true, @islogical);
    p.addParameter('lineArgs', {}, @iscell);
    p.addParameter('lineAlpha', 1, @isscalar);
    p.addParameter('Clipping', true, @islogical);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    x = makecol(x);
    ym = makecol(ym);
    ye = makecol(ye);
    
    y1 = ym - ye;
    y2 = ym + ye;
    
    if isempty(p.Results.axh)
        axh = newplot;
    else
        axh = p.Results.axh;
    end
   
    hold(axh, 'on');
    hs = TrialDataUtilities.Plotting.errorshadeInterval(x, y1, y2, color, 'axh', axh, 'Clipping', p.Results.Clipping, p.Unmatched);
    
    if p.Results.showLine
        if p.Results.lineAlpha < 1
            % use patchline for drawing translucent lines
            hl = TrialDataUtilities.Plotting.patchline(x, ym, ...
               'EdgeColor', color, 'EdgeAlpha', p.Results.lineAlpha, ...
               'z', z, p.Results.lineArgs{:});
            TrialDataUtilities.Plotting.hideInLegend(hl);
        else
            % use plot for opaque lines
%             if z ~=0
%                 zv = z*ones(size(v));
%                 hl = plot(x, ym, zv, 'Color', color, 'Parent', axh, p.Results.lineArgs{:});
%             else
                hl = plot(x, ym, 'Color', color, 'Parent', axh, 'Clipping', p.Results.Clipping, p.Results.lineArgs{:});
%             end
        end
    else
        hl = [];
    end
end
