function [hMean, hEllipse, ellipseX, ellipseY] = plotCovarianceEllipse(x, y, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('Color', 'k', @(x) ischar(x) || isvector(x));
    p.addParamValue('LineWidth', 2, @isscalar);
    p.addParamValue('MarkerSize', 8, @isscalar);
    p.parse(varargin{:});
    
    axh = p.Results.axh;
    
    mask = ~isnan(x) & ~isnan(y);
    x = x(mask);
    y = y(mask);
    C = cov(x,y);
    muX = mean(x);
    muY = mean(y);

    n=100; % Number of points around ellipse
    t=0:pi/n:2*pi; % angles around a circle

    [eigvec,eigval] = eig(C); % Compute eigen-stuff
    xy = [cos(t'),sin(t')] * sqrt(eigval) * eigvec'; % Transformation
    ellipseX = xy(:,1) + muX;
    ellipseY = xy(:,2) + muY;
    
    hMean = plot(axh, muX, muY, 'Marker', 'o', 'MarkerSize', p.Results.MarkerSize, ...
        'MarkerFaceColor', p.Results.Color, ...
        'MarkerEdgeColor', 'none');

    hEllipse = plot(axh, ellipseX, ellipseY, '-', ...
        'LineWidth', p.Results.LineWidth, ...
        'Color', p.Results.Color);
    
    TrialDataUtilities.Plotting.hideInLegend([hMean, hEllipse]);
end