function hLine = drawTickRaster(timesCell, varargin)

    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('color', 'k', @(x) ischar(x) || isvector(x));
    p.addParamValue('lineWidth', 1, @isscalar);
    p.addParamValue('xOffset', 0, @isscalar);
    p.addParamValue('yOffset', 0, @isscalar);
    p.addParamValue('rowHeight', 1, @isscalar);
    p.addParamValue('tickHeight', 0.99, @isscalar);
    
    p.parse(varargin{:});
    
    rowHeight = p.Results.rowHeight;
    tickHeight = p.Results.tickHeight;

    nTrials = numel(timesCell);
    
    % build line commands
    XByTrial = cell(1, nTrials);
    YByTrial = cell(1, nTrials);
    for iE = 1:nTrials 
        if ~isempty(timesCell{iE})
            XByTrial{iE} = repmat(makerow(timesCell{iE}), 3, 1);
            XByTrial{iE}(3, :) = NaN;
            YByTrial{iE} = repmat([-rowHeight*(iE-1); -rowHeight*(iE-1)-tickHeight; NaN], 1, numel(timesCell{iE}));
        end
    end
            
    X = cell2mat(XByTrial) + p.Results.xOffset;
    Y = cell2mat(YByTrial) + p.Results.yOffset;
            
    % filter within time limits?
    if ~isempty(X)
        hLine = plot(X(:), Y(:), 'Parent', p.Results.axh, 'Color', p.Results.color, ...
            'LineWidth', p.Results.lineWidth);
    else
        hLine = NaN;
    end
end

