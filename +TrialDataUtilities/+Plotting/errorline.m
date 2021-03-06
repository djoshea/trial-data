function [hMean, hError] = errorline(x,y,e, varargin)
% plots line plots with vertical symmetric errorbars 

    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('Color', 'k', @(x) true);
    p.addParamValue('LineWidth', 1, @isscalar);
    p.addParameter('MarkerSize', 0, @isscalar);
    p.parse(varargin{:});

    x = makecol(x);
    y = makecol(y);
    e = makecol(e);
   
    axh = p.Results.axh;
    if isempty(axh)
        axh = gca;
    end

    if p.Results.MarkerSize == 0
        hMean = plot(x,y, '-', 'Parent', axh, 'Color', p.Results.Color, ...
            'LineWidth', p.Results.LineWidth);
    else
        hMean = plot(x,y, '-o', 'Parent', axh, 'MarkerSize', p.Results.MarkerSize, 'Color', p.Results.Color, ...
            'MarkerFaceColor', p.Results.Color, ...
            'LineWidth', p.Results.LineWidth);
    end
    hold(axh, 'on');
    hError = plot([x'; x'], [y'-e'; y'+e'], 'Parent', axh, ...
        'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
    TrialDataUtilities.Plotting.hideInLegend(hError);
end
