function [traceCenters, hLines] = plotStackedTraces(tvec, mat, varargin)
% [traceCenters] = plotStacked(tvec, mat, varargin)
% Plots the columns of mat stacked vertically on top of each other
%
% For mat as numeric tensor, tvec is vector matching size(mat, 1)
%   Dim 1 of mat is time, dim 2 is traces to stack vertically.
%   If dim 3 has length > 1, these traces will be superimposed at the same
%   location with the colors set by colorMap
%
% For mat as cell matrix, tvec is cell with same size
%   Dim 1 of mat is traces to be stacked vertically
%   Dim 2 is traces to be superimposed
%   Inside each cell time is dim 1, multiple traces with same color is cols
%
% options include all options for plot(...) plus:
%   normalize: [false] scale each signal to the same range before plotting
%   spacingFraction: [1.02] space each trace by this fraction of the previous
%     trace's range

p = inputParser();
p.addParamValue('evenSpacing', false, @islogical);
p.addParamValue('normalize', false, @islogical);
p.addParamValue('intercalate', false, @islogical);
p.addParamValue('spacingFraction', 1.2, @isscalar);
p.addParamValue('colormap', @TrialDataUtilities.Color.hslmap, @(x) isa(x, 'function_handle') || ismatrix(x));
p.addParamValue('maintainScaleSuperimposed', true, @islogical); % when superimposing multiple traces, keep the relative size and offset between the superimposed traces
p.addParameter('labels', {}, @isvector);
p.addParameter('labelsLinesWithinEachTrace', {}, @iscell);
p.addParameter('showLabels', 'auto', @(x) islogical(x) || ischar(x));
p.addParameter('clickable', true, @islogical);
p.KeepUnmatched = true;
p.parse(varargin{:});

if ~iscell(mat)
    % all traces share common time vector
    nTraces = size(mat, 2);
    nSuperimposed =size(mat, 3);
    nTime = size(mat, 1);
    if isempty(tvec)
        tvec = makecol(1:nTime);
    end
    assert(numel(tvec) == nTime, 'Time vector must match size(mat, 1)');
    
    if ~isfloat(mat)
        mat = single(mat);
    end
else
    % everytraces has different time vector
    nTraces = size(mat, 1);
    nSuperimposed = size(mat, 2);
end
   
% construct labels by traces
if isempty(p.Results.labels)
    labels = arrayfun(@num2str, 1:nTraces, 'UniformOutput', false);
elseif isnumeric(p.Results.labels)
    labels = arrayfun(@num2str, p.Results.labels, 'UniformOutput', false);
else
    labels = p.Results.labels;
end
labels = makecol(labels);

% construct labels for each superimposed line within each trace
if isempty(p.Results.labelsLinesWithinEachTrace)
    labelsLinesWithinEachTrace = arrayfun(@num2str, 1:nSuperimposed, 'UniformOutput', false);
else
    labelsLinesWithinEachTrace = p.Results.labelsLinesWithinEachTrace;
end

if ~iscell(mat)
    % invert the order of the traces so the first is plotted at the top
    if p.Results.maintainScaleSuperimposed
        % subtract the min so each group of traces has min == 0
        minEachGroup = TensorUtils.nanminMultiDim(mat, [1 3]);
        matShift = bsxfun(@minus, mat, minEachGroup);

        if p.Results.normalize
            maxEachGroup = TensorUtils.nanmaxMultiDim(matShift, [1 3]);
            matShift = bsxfun(@rdivide, matShift, maxEachGroup);
        end
    else
        % subtract the min so each trace has min at zero
        matShift = bsxfun(@minus, mat, nanmin(mat, [], 1));

        if p.Results.normalize
            matShift = bsxfun(@rdivide, matShift, range(matShift, 1));
        end
    end

    % compute the max range each row
    ranges = TensorUtils.nanmaxMultiDim(matShift, [1 3]);

    % figure out where each trace should start
    if p.Results.evenSpacing
        % trace(k) will be offset by spacing * k
        if p.Results.intercalate
            deltas = matShift(:, 2:end, :) * p.Results.spacingFraction - matShift(:, 1:end-1, :);
            maxDeltas = TensorUtils.nanmaxMultiDim(deltas, [1 3]); % max over time and superimposed traces
            traceOffsets = (nTraces-1:-1:0) * nanmax(maxDeltas);
        else
            rangesPadded = makerow(ranges * (p.Results.spacingFraction));
            traceOffsets = (nTraces-1:-1:0) * nanmax(rangesPadded);
        end
    else
        if p.Results.intercalate
            % the offsets should be set so that pairs of points across
            % all time points should satisfy:
            %   traceN+1(t) * spacingFraction <= traceN(t) + offset
            % or:
            %   offset = max_t (traceN+1(t) * spacingFraction - traceN(t))

            deltas = matShift(:, 2:end, :) * p.Results.spacingFraction - matShift(:, 1:end-1, :);
            maxDeltas = TensorUtils.nanmaxMultiDim(deltas, [1 3]); % max over time and superimposed traces
            cs = fliplr(cumsum(fliplr(maxDeltas)));
            traceOffsets = [cs, 0];

        else 
            rangesPadded = makerow(ranges * (p.Results.spacingFraction));
            cs = fliplr(cumsum(fliplr(rangesPadded)));
            traceOffsets = [cs(2:end), 0];
        end
    end
    traceCenters = (traceOffsets + ranges / 2)';

    matShift = bsxfun(@plus, matShift, traceOffsets);

    % expand colormap to be exactly nSuperimposed long
    map = TrialDataUtilities.Plotting.expandWrapColormap(p.Results.colormap, nSuperimposed);

    if nSuperimposed == 1
        % plot simultaneously
        hLines = plot(tvec, matShift, '-', 'Color', map(1, :), p.Unmatched);
    else
        set(gca, 'ColorOrder', map); % the map will superimpose automatically
        hold on;

        % here we arrange so that all traces are stacked vertically but that
        % all positions along dim 2 are grouped together 
        matShiftCat = TensorUtils.reshapeByConcatenatingDims(matShift, {1 [3 2]});

        hLines = plot(tvec, matShiftCat, '-', p.Unmatched);
        hLines = reshape(hLines, [nSuperimposed nTraces])';
    end

else
    % cell mode - everyone has a different time vector
    if ~iscell(tvec)
        tvec = repmat({tvec}, nTraces, nSuperimposed);
    end
    
    map = TrialDataUtilities.Plotting.expandWrapColormap(p.Results.colormap, nSuperimposed);
    
    % invert the order of the traces so the first is plotted at the top
%     mat = flipud(mat);
%     labels = flipud(labels);
%     tvec = flipud(tvec);

    if p.Results.maintainScaleSuperimposed
        % subtract the min so each group trace has min at zero
        minEachRow = double(nanmin(cellfun(@(mat) nanminNanEmpty(mat(:)), mat), [], 2));
        minCell = num2cell(repmat(minEachRow, 1, size(mat, 2)));
        cellShift = cellfun(@(mat, min) mat - min, mat, minCell, 'UniformOutput', false);

        if p.Results.normalize
            maxEachRow = nanmax(cellfun(@(mat) nanmax(mat(:)), mat), [], 2);
            maxCell = num2cell(repmat(maxEachRow, 1, size(mat, 2)));
            cellShift = cellfun(@(matShift, max) matShift ./ max, cellShift, maxCell, 'UniformOutput', false);
        end
    else
        
        % subtract the min so each trace has min at zero
        cellShift = cellfun(@(mat) bsxfun(@minus, mat, nanmin(mat, [], 1)), mat, 'UniformOutput', false);

        if p.Results.normalize
            cellShift = cellfun(@(matShift) bsxfun(@rdivide, matShift, range(matShift, 1)), 'UniformOutput', false);
        end
    end

    % compute the max range each row
    ranges = max(cellfun(@(matShift) nanmaxNanEmpty(matShift(:), [], 1), cellShift), [], 2);

    ranges(isnan(ranges)) = 0;
    
    % can't intercalate without doing time consuming interpolation
    rangesPadded = ranges * (p.Results.spacingFraction);
    cs = flipud(cumsum(flipud(rangesPadded)));
    traceOffsets = [cs(2:end); 0];

    traceOffsets = traceOffsets - min(traceOffsets);
    traceCenters = (traceOffsets + ranges / 2)';

    hLines = cell(nTraces, nSuperimposed);
    %lineDescriptionsCell = cell(nTraces, nSuperimposed);
    for iT = 1:nTraces
        for iS = 1:nSuperimposed
            matShift = cellShift{iT, iS} + traceOffsets(iT);
            if isempty(tvec{iT, iS})
                tvecThis = 1:numel(matShift);
            else
                tvecThis = tvec{iT, iS};
            end

            hLines{iT, iS} = plot(tvecThis, matShift, '-', 'Color', map(iS, :), p.Unmatched);
            hold on;
        end
    end
end

ylim([0 traceOffsets(1) + ranges(1)]);

if strcmp(p.Results.showLabels, 'auto')
    showLabels = nTraces < 25;
else
    showLabels = p.Results.showLabels;
end

if showLabels
    au = AutoAxis(gca);
    spans = [makerow(traceOffsets); makerow(traceOffsets + ranges)];
    au.addLabeledSpan('y', 'span', spans, 'label', labels);
    au.addAutoAxisX();
    au.xlabel('Time');
    au.update();
end

if p.Results.clickable
    % build line handle descriptions
    lineDescriptions = cell(nTraces, nSuperimposed);
    if iscell(hLines)
        for iT = 1:nTraces
            for iS = 1:nSuperimposed
                set(hLines{iT, iS}, 'Description', sprintf('Trace %s, Line %s', labels{iT}, labelsLinesWithinEachTrace{iS}));
            end
        end
        hvec = cat(1, hLines{:});
    else
        for iT = 1:nTraces
            for iS = 1:nSuperimposed
                set(hLines(iT, iS), 'Description', sprintf('Trace %s, Line %s', labels{iT}, labelsLinesWithinEachTrace{iS}));
            end
        end
        hvec = hLines(:);
    end

    TrialDataUtilities.Plotting.makeClickableShowDescription(hvec);
end

% hold(ifelse(oldhold, 'on', 'off'));
hold off;

end

function r = nanmaxNanEmpty(v1, varargin)
    if isempty(v1) && numel(varargin) >= 1 && isempty(varargin{1})
        r = NaN;
    else
        r = nanmax(v1, varargin{:});
    end
end

