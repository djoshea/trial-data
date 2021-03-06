function [hPoints, hError] = plotScatterErrorBars(xData, yData, varargin)

iscolor = @(x) ischar(x) || (isnumeric(x) && isvector(x) && numel(x) == 3);
p = inputParser();
% specify these:
p.addParameter('xError', [], @(x) true);
p.addParameter('yError', [], @(x) true);

% or these:
p.addParameter('xConfHigh', [], @(x) true);
p.addParameter('yConfHigh', [], @(x) true);
p.addParameter('xConfLow', [], @(x) true);
p.addParameter('yConfLow', [], @(x) true);

p.addParameter('axh', [], @(x) isempty(x) || ishandle(x));
p.addParameter('edgeColor', 'none', iscolor);
p.addParameter('edgeAlpha', 0.6, @isscalar);
p.addParameter('color', [0.1 0.1 0.1], @(x) isempty(x) || iscolor(x));
p.addParameter('alpha', 0.8, @isscalar);
p.addParameter('markerSize', 3, @isscalar);
p.addParameter('errorLineWidth', 1, @isscalar);
p.addParameter('errorLineAlpha', 0.8, @isscalar);
p.addParameter('errorLineColor', [], @(x) isempty(x) || iscolor(x));
p.CaseSensitive = false;
p.parse(varargin{:});

axh = p.Results.axh;
if isempty(axh) 
    axh = newplot();
end

xData = xData(:);
yData = yData(:);

% default to no error
xErrorLow = xData;
xErrorHigh = xData;
% use confidence intervals if specified
if ~isempty(p.Results.xConfHigh)
    xErrorHigh = p.Results.xConfHigh(:);
    xErrorLow = p.Results.xConfLow(:);
% use error bars if specified
elseif ~isempty(p.Results.xError)
    xErrorHigh = xData + p.Results.xError(:);
    xErrorLow = xData - p.Results.xError(:);
end

yErrorLow = yData;
yErrorHigh = yData;
if ~isempty(p.Results.yConfHigh)
    yErrorHigh = p.Results.yConfHigh(:);
    yErrorLow = p.Results.yConfLow(:);
elseif ~isempty(p.Results.yError)
    yErrorHigh = yData + p.Results.yError(:);
    yErrorLow = yData - p.Results.yError(:);
end

xLine = [makerow(xErrorLow), makerow(xData); makerow(xErrorHigh), makerow(xData)];
yLine = [makerow(yData), makerow(yErrorLow); makerow(yData), makerow(yErrorHigh)];

if isempty(p.Results.errorLineColor)
    errorLineColor = 0.5 * p.Results.color + 0.5;
else
    errorLineColor = p.Results.errorLineColor;
end

hold(axh, 'on');
hError = line(xLine, yLine, 'LineWidth', p.Results.errorLineWidth, ...
    'Color', errorLineColor, 'Parent', axh);
TrialDataUtilities.Plotting.setLineOpacity(hError, p.Results.errorLineAlpha);
TrialDataUtilities.Plotting.hideInLegend(hError);

% draw points on top
hPoints = plot(xData, yData, 'o', 'Parent', axh, ...
    'MarkerSize', p.Results.markerSize, ...
    'MarkerFaceColor', p.Results.color, ...
    'MarkerEdgeColor', p.Results.edgeColor);
TrialDataUtilities.Plotting.setMarkerOpacity(hPoints, p.Results.alpha, p.Results.edgeAlpha);


end


