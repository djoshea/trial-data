function [hPoints, hError] = plotScatterErrorBars(xData, yData, varargin)

iscolor = @(x) isstringlike(x) || (isnumeric(x) && isvector(x) && numel(x) == 3);
p = inputParser();
% specify these:
p.addParameter('xError', [], @(x) true);
p.addParameter('yError', [], @(x) true);

% or these (N x 2)
p.addParameter('yCI', [], @(x) true);
p.addParameter('xCI', [], @(x) true);

% or these:
p.addParameter('xConfHigh', [], @(x) true);
p.addParameter('yConfHigh', [], @(x) true);
p.addParameter('xConfLow', [], @(x) true);
p.addParameter('yConfLow', [], @(x) true);

p.addParameter('axh', [], @(x) isempty(x) || ishandle(x));
p.addParameter('edgeColor', 'none', iscolor);
p.addParameter('edgeWidth', 1, @isscalar);
p.addParameter('edgeAlpha', 0.6, @isscalar);
p.addParameter('color', [0.1 0.1 0.1], @(x) isempty(x) || iscolor(x) || ismatrix(x));
p.addParameter('alpha', 0.8, @isscalar);
p.addParameter('markerSize', 20, @isscalar);
p.addParameter('errorLineWidth', 1, @isscalar);
p.addParameter('errorLineAlpha', 0.8, @isscalar);
p.addParameter('errorLineColor', [], @(x) isempty(x) || iscolor(x) || ismatrix(x));

p.addParameter('Clipping', 'on', @(x) true);
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
    
elseif ~isempty(p.Results.xCI)
    xErrorLow = p.Results.xCI(:, 1);
    xErrorHigh = p.Results.xCI(:, 2);
    
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
    
elseif ~isempty(p.Results.yCI)
    yErrorLow = p.Results.yCI(:, 1);
    yErrorHigh = p.Results.yCI(:, 2);
 
elseif ~isempty(p.Results.yError)
    yErrorHigh = yData + p.Results.yError(:);
    yErrorLow = yData - p.Results.yError(:);
end

assert(numel(yErrorHigh) == numel(yData));
assert(numel(yErrorLow) == numel(yData));
assert(numel(xErrorHigh) == numel(xData));
assert(numel(xErrorLow) == numel(xData));

xLine = [makerow(xErrorLow), makerow(xData); makerow(xErrorHigh), makerow(xData)];
yLine = [makerow(yData), makerow(yErrorLow); makerow(yData), makerow(yErrorHigh)];

if isempty(p.Results.errorLineColor)
    errorLineColor = 0.5 * p.Results.color + 0.5;
else
    errorLineColor = p.Results.errorLineColor;
end

N = numel(xData);
if size(errorLineColor, 1) == 1
    errorLineColor = repmat(errorLineColor, N, 1);
end

if size(errorLineColor, 1) == N
    errorLineColor = repmat(errorLineColor, 2, 1);
end

hold(axh, 'on');
hError = line(xLine, yLine, 'LineWidth', p.Results.errorLineWidth, 'Parent', axh);


for iH = 1:2*N
    hError(iH).Color = errorLineColor(iH, :);
end
    

TrialDataUtilities.Plotting.setLineOpacity(hError, p.Results.errorLineAlpha);
TrialDataUtilities.Plotting.hideInLegend(hError);

% draw points on top
% hPoints = plot(xData, yData, 'o', 'Parent', axh, ...
%     'MarkerSize', p.Results.markerSize, ...
%     'MarkerFaceColor', p.Results.color, ...
%     'MarkerEdgeColor', p.Results.edgeColor);
% TrialDataUtilities.Plotting.setMarkerOpacity(hPoints, p.Results.alpha, p.Results.edgeAlpha);

hPoints = scatter(xData, yData, p.Results.markerSize, p.Results.color, 'filled', ...
    'MarkerEdgeColor', p.Results.edgeColor, 'MarkerFaceAlpha', p.Results.alpha, 'MarkerEdgeAlpha', p.Results.edgeAlpha, 'Clipping', p.Results.Clipping, 'LineWidth', p.Results.edgeWidth);

end


