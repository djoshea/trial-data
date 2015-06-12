function h = plotHistogram(varargin)
% Plots a histogram using stairs.
%
% plotHistogram('freq', countsToPlot, 'binEdges', binEdges, ...)
% plotHistogram('values', rawValueVector, 'binEdges', binEdges, ...)
% plotHistogram('values', rawValueVector, 'nBins', nBins, ...)
%
% Parameters:
% axh : axis to plot into
% color : line color
% lineWidth: line width 
% lineAlpha : line alpha
% fillColor
% fillAlpha
iscolor = @(x) ischar(x) || (isnumeric(x) && isvector(x) && numel(x) == 3);
p = inputParser();
p.addParameter('nBins', 100, @isscalar);
p.addParameter('binEdges', [], @(x) isempty(x) || isvector(x));
p.addParameter('freq', [], @(x) isempty(x) || isvector(x));
p.addParameter('values', [], @(x) isempty(x) || isvector(x));
p.addParameter('axh', [], @(x) isempty(x) || ishandle(x));
p.addParameter('color', [0 0 0], iscolor);
p.addParameter('lineAlpha', 1, @isscalar);
p.addParameter('fillColor', [0.5 0.5 0.5], @(x) isempty(x) || iscolor(x));
p.addParameter('fillAlpha', 0.5, @isscalar);
p.addParameter('lineWidth', 1, @isscalar);
p.CaseSensitive = false;
p.parse(varargin{:});

if isempty(p.Results.binEdges)
    assert(~isempty(p.Results.values), 'Must provide ''values'' parameter');
    
    values = p.Results.values;
    binEdges = linspace(nanmin(values), nanmax(values), p.Results.nBins + 1);
    freq = histc(values, binEdges);
else
    binEdges = p.Results.binEdges;
    if ~isempty(p.Results.freq)
        freq = p.Results.freq;
    elseif ~isempty(p.Results.values)
        freq = histc(p.Results.values, binEdges);
    end
end

axh = p.Results.axh;
if isempty(axh) 
    axh = newplot();
end


binLast = binEdges(end) + binEdges(end) - binEdges(end-1);
binEdges = [makecol(binEdges); binLast];
freqPadded = [makecol(freq); 0];
[xPts, yPts] = stairs(binEdges, freqPadded);
xPts = [xPts(1); xPts; xPts(end)];
yPts = [0; yPts; 0];

if ~isempty(p.Results.fillColor) && ~strcmp(p.Results.fillColor, 'none')
    h.fill = patch(xPts, yPts, p.Results.fillColor, 'FaceAlpha', p.Results.fillAlpha, ...
        'EdgeColor', 'none', 'Parent', axh);
    hold(axh, 'on');
end

if ~isempty(p.Results.color) && ~strcmp(p.Results.color, 'none')
%     h.line = stairs(xPts, yPts, '-', 'LineWidth', p.Results.lineWidth, ...
%         'Color', p.Results.color, ...
%         'Parent', axh);

    [n,nc] = size(yPts);
    ndx = [1:n;1:n];
    y2 = yPts(ndx(1:2*n-1),:);
    if size(xPts,2)==1,
        x2 = xPts(ndx(2:2*n),ones(1,nc));
    else
        x2 = xPts(ndx(2:2*n),:);
    end
    h.line = plot(x2, y2, '-', 'LineWidth', p.Results.lineWidth, ...
        'Color', p.Results.color, ...
        'Parent', axh);
    
    SaveFigure.setLineOpacity(h.line, p.Results.lineAlpha);
end


