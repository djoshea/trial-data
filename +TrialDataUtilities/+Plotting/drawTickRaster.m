function drawTickRaster(timesCell, varargin)

    p = inputParser();
    p.addParamValue('xOffset', 0, @isscalar);
    p.addParamValue('yOffset', 0, @isscalar);
    p.addParamValue('rowHeight', 1, @isscalar);
    p.addParamValue('tickHeight', 0.99, @isscalar);
    p.parse(varargin{:});

    nRows = numel(timesCell);
end

