function hLine = drawTickRaster(timesCell, varargin)

    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('color', 'k', @(x) ischar(x) || isvector(x));
    p.addParamValue('lineWidth', 1, @isscalar);
    p.addParamValue('xOffset', 0, @isscalar);
    p.addParamValue('yOffset', 0, @isscalar);
    p.addParamValue('rowHeight', 1, @isscalar);
    p.addParamValue('tickHeight', 0.99, @isscalar);
    p.addParameter('alpha', 1, @isscalar);
    p.addParameter('waveCell', [], @(x) isempty(x) || iscell(x));
    p.addParameter('waveformTimeRelative', [], @(x) isempty(x) || isvector(x));
    p.addParameter('waveScaleHeight', 1, @isscalar);
    p.addParameter('waveScaleTime', 1, @isscalar);
    %p.addParameter('waveOffsetY', 0, @isscalar); % offset will be applied before waveScaleHeight
    
    p.parse(varargin{:});
    
    rowHeight = p.Results.rowHeight;
    tickHeight = p.Results.tickHeight;

    nTrials = numel(timesCell);
    
    if ~isempty(p.Results.waveCell)
        % plotting waveforms where the ticks would normally be
        tvec = makerow(p.Results.waveformTimeRelative * p.Results.waveScaleTime);
        
        [waveYByTrial, waveTByTrial] = cellvec(nTrials);
        for iE = 1:nTrials
            if isempty(p.Results.waveCell{iE})
                waveYByTrial{iE} = zeros(0, 1);
                waveTByTrial{iE} = zeros(0, 1);
            else
                waves = (p.Results.waveCell{iE} - nanmin(p.Results.waveCell{iE}(:))) * p.Results.waveScaleHeight - rowHeight*iE;
                tvecMat = bsxfun(@plus, repmat(tvec, size(waves, 1), 1), timesCell{iE});
                if size(waves, 2) < size(tvecMat, 2)
                    tvecMat = tvecMat(:, 1:size(waves, 2));
                end
                % cat the waveforms into one long column with NaNs inserted
                % between
                waveYByTrial{iE} = TensorUtils.flatten(cat(2, waves,   nan(size(waves, 1),1) )');
                waveTByTrial{iE} = TensorUtils.flatten(cat(2, tvecMat, nan(size(waves, 1),1) )');
            end
        end
        
        X = cat(1, waveTByTrial{:}) + p.Results.xOffset;
        Y = cat(1, waveYByTrial{:}) + p.Results.yOffset; 
        
        % filter within time limits?
        if ~isempty(X)
            hLine = plot(X(:), Y(:), 'Parent', p.Results.axh, 'Color', p.Results.color, ...
                'LineWidth', p.Results.lineWidth);
            if p.Results.alpha < 1
                TrialDataUtilities.Plotting.setLineOpacity(hLine, p.Results.alpha);
            end
        else
            hLine = NaN;
        end
    else
        % draw vertical ticks
        
        % build line commands
        XByTrial = cell(1, nTrials);
        YByTrial = cell(1, nTrials);
        for iE = 1:nTrials 
            if ~isempty(timesCell{iE})
                XByTrial{iE} = repmat(makerow(timesCell{iE}), 2, 1);
                YByTrial{iE} = repmat([-rowHeight*(iE-1); -rowHeight*(iE-1)-tickHeight], 1, numel(timesCell{iE}));
            end
        end

        X = cell2mat(XByTrial) + p.Results.xOffset;
        Y = cell2mat(YByTrial) + p.Results.yOffset;
        
        % filter within time limits?
        if ~isempty(X)
            hLine = line(X, Y, 'Parent', p.Results.axh, 'Color', p.Results.color, ...
                'LineWidth', p.Results.lineWidth);
            if p.Results.alpha < 1
                TrialDataUtilities.Plotting.setLineOpacity(hLine, p.Results.alpha);
            end
        else
            hLine = NaN;
        end

    end
    
    
end

