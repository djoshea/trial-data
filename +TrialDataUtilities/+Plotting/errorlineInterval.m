function [hMean, hError] = errorlineInterval(x,y,lo,hi, varargin)
% plots line plots with vertical symmetric errorbars 

    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('Color', 'k', @(x) true);
    p.addParamValue('LineWidth', 1, @isscalar);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    x = makecol(x);
    y = makecol(y);
    lo = makecol(lo);
    hi = makecol(hi);
   
    axh = p.Results.axh;
    if isempty(axh)
        axh = gca;
    end

    hMean = plot(x,y, '-', 'Parent', axh, 'Color', p.Results.Color, ...
        'LineWidth', p.Results.LineWidth, p.Unmatched);
    hold(axh, 'on');
    hError = plot([x'; x'], [lo'; hi'], 'Parent', axh, ...
        'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
    TrialDataUtilities.Plotting.hideInLegend(hError);
end
