function plotHistogram(binEdges, freq, varargin)

iscolor = @(x) ischar(x) || (isnumeric(x) && isvector(x) && numel(x) == 3);
p = inputParser();
p.addParamValue('axh', [], @(x) isempty(x) || ishandle(x));
p.addParamValue('color', [0 0 0], iscolor);
p.addParamValue('fillColor', [0.5 0.5 0.5], @(x) isempty(x) || iscolor(x));
p.addParamValue('fillAlpha', 0.5, @isscalar);
p.addParamValue('lineWidth', 1, @isscalar);
p.CaseSensitive = false;
p.parse(varargin{:});

axh = p.Results.axh;
if isempty(axh) 
    axh = newplot();
end

[xPts, yPts] = stairs(binEdges, freq);
xPts = [xPts(1); xPts; xPts(end)];
yPts = [0; yPts; 0];

if ~isempty(p.Results.fillColor) && ~strcmp(p.Results.fillColor, 'none')
    patch(xPts, yPts, p.Results.fillColor, 'FaceAlpha', p.Results.fillAlpha, ...
        'EdgeColor', 'none', 'Parent', axh);
    hold(axh, 'on');
end

if ~isempty(p.Results.color) && ~strcmp(p.Results.color, 'none')
    stairs(xPts, yPts, '-', 'LineWidth', p.Results.lineWidth, ...
        'Color', p.Results.color, ...
        'Parent', axh);
end


