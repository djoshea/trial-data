function hLine = drawTickRaster(timesCell, varargin)

    p = inputParser();
    p.addParameter('axh', [], @(x) true);
    p.addParameter('color', 'k', @(x) ischar(x) || isvector(x));
    p.addParameter('lineWidth', 1, @isscalar);
    p.addParameter('xOffset', 0, @isscalar);
    p.addParameter('yOffset', 0, @isscalar);
    p.addParameter('rowHeight', 1, @isscalar);
    p.addParameter('tickHeight', 0.99, @isscalar);
    p.addParameter('alpha', 1, @isscalar);
    p.addParameter('waveCell', [], @(x) isempty(x) || iscell(x));
    p.addParameter('waveformTimeRelative', [], @(x) isempty(x) || isvector(x));
    p.addParameter('normalizeWaveforms', true, @islogical); % if true, treat as unnormalized and scale them, if false, treat the values as normalized to [0 1]
    p.addParameter('waveScaleHeight', 1, @isscalar);
    p.addParameter('waveScaleTime', 1, @isscalar);
    p.addParameter('quick', false, @islogical); % use ugly dots instead
    p.addParameter('MarkerSize', 3, @isscalar); % only with quick
    
    p.addParameter('datatipMetaBySpike', {}, @(x) iscell(x) || isstruct(x)); % cell same size as timesCell of structs with same length as timesCell{.}
    p.addParameter('datatipMetaByTrial', {}, @isstruct); % nTrials x 1 struct of per-trial meta data
    
    %p.addParameter('waveOffsetY', 0, @isscalar); % offset will be applied before waveScaleHeight
    
    p.parse(varargin{:});
    
    axh = p.Results.axh;
    if isempty(axh)
        axh = gca;
    end
    
    rowHeight = p.Results.rowHeight;
    tickHeight = p.Results.tickHeight;
    nTrials = numel(timesCell);
    datatipMetaBySpike = p.Results.datatipMetaBySpike;
    if iscell(datatipMetaBySpike)
        datatipMetaBySpike = cat(1, datatipMetaBySpike{:});
    end
    datatipMetaByTrial = p.Results.datatipMetaByTrial;
    
    if p.Results.normalizeWaveforms && ~isempty(p.Results.waveCell)
        minWave = nanmin(cellfun(@(w) nanmin(w(:)), p.Results.waveCell));
        maxWave = nanmax(cellfun(@(w) nanmax(w(:)), p.Results.waveCell));
    end
    
    timesCell = makecol(timesCell);
    
    nSpikesByTrial = cellfun(@numel, timesCell);
    emptyMask = nSpikesByTrial == 0;
    
    if ~isempty(p.Results.waveCell) && ~p.Results.quick
        % plotting waveforms where the ticks would normally be
        tvec = makerow(p.Results.waveformTimeRelative * p.Results.waveScaleTime);
        
        [waveYByTrial, waveTByTrial] = cellvec(nTrials);
        for iE = 1:nTrials
            if isempty(p.Results.waveCell{iE})
                waveYByTrial{iE} = zeros(0, 1);
                waveTByTrial{iE} = zeros(0, 1);
            else
                waves = p.Results.waveCell{iE};
                if p.Results.normalizeWaveforms
                    waves = (waves - minWave)/ (maxWave-minWave) * p.Results.waveScaleHeight;
                else
                    % regard as [0 1] scaled already (for common scaling
                    % purposes) and just multiply by height
                    waves = waves * p.Results.waveScaleHeight;
                end
                waves = waves - rowHeight*iE;
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
        X = X(:);
        Y = Y(:);
        
        % for metadata plotting
        nPointsPerSpike = numel(tvec);
        
        % filter within time limits?
        if ~isempty(X)
            hLine = plot(X, Y, 'Parent', axh, 'Color', p.Results.color, ...
                'LineWidth', p.Results.lineWidth);
            if p.Results.alpha < 1
                TrialDataUtilities.Plotting.setLineOpacity(hLine, p.Results.alpha);
            end
        else
            hLine = gobjects(1,1);
        end
        
    elseif ~p.Results.quick
        % draw vertical ticks

        % build line commands by intercalating a NaN in between each spike (upper, lower, NaN, upper, lower, NaN, ...)
        XByTrial = cell(1, nTrials);
        YByTrial = cell(1, nTrials);
        for iE = 1:nTrials 
            if ~emptyMask(iE)
                XByTrial{iE} = [repmat(makerow(timesCell{iE}), 2, 1); nan(1, numel(timesCell{iE}))];
                YByTrial{iE} = repmat([-rowHeight*(iE-1); -rowHeight*(iE-1)-tickHeight; NaN], 1, numel(timesCell{iE}));
            end
        end

        X = cell2mat(XByTrial(~emptyMask)) + p.Results.xOffset;
        Y = cell2mat(YByTrial(~emptyMask)) + p.Results.yOffset;
        
        % turn into columns so only one 1 is plotted
        X = X(:);
        Y = Y(:);
        nPointsPerSpike = 3;
        
        % filter within time limits?
        if ~isempty(X)
            hLine = plot(X, Y, 'Parent', axh, 'Color', p.Results.color, ...
                'LineWidth', p.Results.lineWidth);
            if p.Results.alpha < 1
                TrialDataUtilities.Plotting.setLineOpacity(hLine, p.Results.alpha);
            end
        else
            hLine = gobjects(1,1);
        end
 
    else
        % draw dots quickly
        % build line commands
        XByTrial = cell(1, nTrials);
        YByTrial = cell(1, nTrials);
        for iE = 1:nTrials 
            if ~isempty(timesCell{iE})
                XByTrial{iE} = makecol(timesCell{iE});
                YByTrial{iE} = repmat(-rowHeight*(iE-1) - tickHeight/2, numel(timesCell{iE}), 1);
            end
        end

        X = cat(1, XByTrial{:}) + p.Results.xOffset;
        Y = cat(1, YByTrial{:}) + p.Results.yOffset;
        nPointsPerSpike = 1;

        hLine = plot(X, Y, '.', 'Color', p.Results.color, 'MarkerSize', p.Results.markerSize);
    end
    
    if ~isempty(datatipMetaBySpike) || ~isempty(datatipMetaByTrial)
        catc = @(x) cat(1, x{:});
        trialInd = catc(arrayfun(@(ind, n) repmat(ind, n, 1), (1:nTrials)', nSpikesByTrial, 'UniformOutput', false));
%         indWithinTrial = catc(arrayfun(@(n) (1:n)', nSpikesByTrial, 'UniformOutput', false));
        trialInd = repelem(trialInd, nPointsPerSpike, 1);
%         indWithinTrial = repelem(indWithinTrial, nPointsPerSpike, 1);
    end
    
    if ~isempty(datatipMetaBySpike) && ~isa(hLine, 'matlab.graphics.GraphicsPlaceholder')
        flds = fieldnames(datatipMetaBySpike);
        for iF = 1:numel(flds)
            valBySpike = cat(1, datatipMetaBySpike.(flds{iF}));
            r = dataTipTextRow(label, valBySpike);
            hLine.DataTipTemplate.DataTipRows(end+1) = r;
        end
    end
    
    if ~isempty(datatipMetaByTrial) && ~isa(hLine, 'matlab.graphics.GraphicsPlaceholder')
        flds = fieldnames(datatipMetaByTrial);
        for iF = 1:numel(flds)
            valByTrial = cat(1, datatipMetaByTrial.(flds{iF}));
            valBySpike = valByTrial(trialInd);
            r = dataTipTextRow(flds{iF}, valBySpike);
            hLine.DataTipTemplate.DataTipRows(end+1) = r;
        end
    end
   
end

