function [hMean, hError] = errorlineCI(x,y,yci, varargin)
% plots line plots with vertical symmetric errorbars 

    p = inputParser();
    p.addParameter('axh', gca, @ishandle);
    p.addParameter('Color', 'k', @(x) true);
    p.addParameter('LineWidth', 1, @isscalar);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    x = makecol(x);
    y = makecol(y);
    lo = makecol(yci(:, 1));
    hi = makecol(yci(:, 2));
   
    washolding = ishold;
    hMean = plot(x,y, '-', 'Color', p.Results.Color, ...
        'LineWidth', p.Results.LineWidth, p.Unmatched);
    hold on;
    hError = plot([x'; x'], [lo'; hi'], ...
        'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
    TrialDataUtilities.Plotting.hideInLegend(hError);
    
    if ~washolding
        hold off;
    end
    
end
