function [hEllipse, ellipseX, ellipseY] = plotEllipseFromCovariance(C, varargin)
    p = inputParser();
    p.addParameter('center', [0 0]);
    p.addParameter('dims', [1 2]);
    p.addParameter('axh', gca, @ishandle);
    p.addParameter('Color', 'k', @(x) ischar(x) || isvector(x));
    p.addParameter('LineWidth', 2, @isscalar);
    p.addParameter('LineStyle', '-');
    p.addParameter('alpha', 1, @isscalar);
    p.addParameter('MarkerSize', 8, @isscalar);
    p.parse(varargin{:});
    
    axh = p.Results.axh;
    
    n=100; % Number of points around ellipse
    t=0:pi/n:2*pi; % angles around a circle

    dims = p.Results.dims;
    C = C(dims, dims);

    [eigvec,eigval] = eig(C); % Compute eigen-stuff
    xy = [cos(t'),sin(t')] * sqrt(eigval) * eigvec'; % Transformation
    ellipseX = xy(:,1) + p.Results.center(1);
    ellipseY = xy(:,2) + p.Results.center(2);
    
    hEllipse = plot(axh, ellipseX, ellipseY, '-', ...
        'LineWidth', p.Results.LineWidth, ...
        'Color', p.Results.Color, 'LineStyle', p.Results.LineStyle);

    if p.Results.alpha < 1
        TrialDataUtilities.Plotting.setLineOpacity(hEllipse, p.Results.alpha);
    end
    
    TrialDataUtilities.Plotting.hideInLegend([hEllipse]);
end