function [hl, hs] = errorshadeCI(x, ym, yCI, color, varargin)
% [hl, hs] = errorshade(x, ym, ye, color, varargin)

    p = inputParser();
    p.addParameter('axh', [], @(x) true);
    p.addParameter('showLine', true, @islogical);
    p.addParameter('lineArgs', {}, @iscell);
    p.addParameter('lineAlpha', 1, @isscalar);
    p.addParameter('DisplayName', "", @isstringlike);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    y1 = yCI(:, 1);
    y2 = yCI(:, 2);
    
    if isempty(p.Results.axh)
        axh = newplot;
    else
        axh = p.Results.axh;
    end
    
    washolding = ishold(axh);
   
    hold(axh, 'on');
    hs = TrialDataUtilities.Plotting.errorshadeInterval(x, y1, y2, color, 'axh', axh, p.Unmatched);
    
    if p.Results.showLine
        hl = plot(x, ym, 'Color', color, 'Parent', axh, p.Results.lineArgs{:}, 'DisplayName', p.Results.DisplayName);
        TrialDataUtilities.Plotting.setLineOpacity(hl, p.Results.lineAlpha);
        hl = [];
    end
    
    if ~washolding
        hold(axh, 'off');
    end
end
