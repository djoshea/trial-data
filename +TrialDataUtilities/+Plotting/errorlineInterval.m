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
        axh_args = {};
        hold_args = {};
        newplot;
    else
        hold_args = {axh};
        axh_args = {'Parent', axh};
    end

    hMean = plot(x,y, '-', axh_args{:}, 'Color', p.Results.Color, ...
        'LineWidth', p.Results.LineWidth, p.Unmatched);
    hold(hold_args{:}, 'on');
    hError = plot([x'; x'], [lo'; hi'], axh_args{:}, ...
        'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
    TrialDataUtilities.Plotting.hideInLegend(hError);
    hold(hold_args{:}, 'off');
    
end
