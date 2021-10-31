function [traceCenters, hLines] = plotStackedTraces(tvec, data, varargin)
% [traceCenters] = plotStacked(tvec, data, varargin)
% Plots the columns of data stacked vertically on top of each other
%
% For data as numeric tensor, tvec is vector matching size(data, 1)
%   Dim 1 of data is time, dim 2 is traces to stack vertically [nTraces]
%   If dim 3 has length > 1, these traces will be superimposed at the same
%   location with the colors set by colorMap [nTracesSuperimposed]
%
% For data as cell matrix, tvec is cell with same size
%   Dim 1 of data is traces to be stacked vertically [nTraces]
%   Dim 2 is traces to be superimposed [nTracesSuperimposed]
%   Inside each cell time is dim 1, multiple traces with same color can be plotted over cols
%
% nTraces is 
%
% options include all options for plot(...) plus:
%   normalize: [false] scale each signal to the same range before plotting
%   spacingFraction: [1.02] space each trace by this fraction of the previous
%     trace's range

% TODO : fix intercalate - probably not working anymore

p = inputParser();
p.addParameter('data_ci', [], @(x) isnumeric(x) || iscell(x) || isempty(x)); % time x stack x superimposed dims x 2 (one more dim than data), or cell with the same
p.addParameter('ci_alpha', 0.5, @isscalar);
p.addParameter('evenSpacing', false, @islogical); % the vertical space allocated to each stacked trace is the same?
p.addParameter('normalize', false, @islogical); % the vertical height of each trace is normalized? or in original data units
p.addParameter('gain', 1, @isscalar);
p.addParameter('intercalate', false, @islogical); % the traces should be squished together as close as possible without touching
p.addParameter('spacingFraction', 1.2, @isscalar); % the gap between each trace / the height of those traces
p.addParameter('colormap', [], @(x) isempty(x) || isa(x, 'function_handle') || ismatrix(x)); % for superimposed traces 
p.addParameter('colormapStacked', [], @(x) isempty(x) || isa(x, 'function_handle') || ismatrix(x)); % for stacked traces, only used if colormap is empty
p.addParameter('maintainScaleSuperimposed', true, @islogical); % when superimposing multiple traces, keep the relative size and offset between the superimposed traces
p.addParameter('labelPrefix', '', @isstringlike);
p.addParameter('labels', {}, @(x) isempty(x) || isvector(x)); % labels over nTraces for the y axis
p.addParameter('labelRotation', 0, @isvector);
p.addParameter('labelsSuperimposed', {}, @isstringlike); % labels over the nSuperimposed traces, for clickable descriptions
p.addParameter('legendUniqueLabelsOnly', false, @islogical);
p.addParameter('showLabels', 'auto', @(x) islogical(x) || ischar(x)); % show the labels on the left axis, slow if too many traces, 'auto' is true if nTraces < 25
p.addParameter('clickable', false, @islogical); % make each trace clickable and show a description
p.addParameter('timeUnits', '', @ischar); 
p.addParameter('timeScaleBar', false, @islogical); % use scale bar instead of tick bridge for time axis?
p.addParameter('dataUnits', [], @(x) ischar(x) || (isvector(x) && iscellstr(x))); % either a string describing units for all traces, or a nTraces x 1 cell of units for each set of traces running vertically
p.addParameter('verticalScaleBarHideLabel', false, @islogical);
p.addParameter('showVerticalScaleBars', false, @(x) islogical(x) || ischar(x)); % show intelligent y axis scale bars on the right hand side
p.addParameter('showDataRanges', false, @(x) islogical(x) || ischar(x)); % show intelligent y axis scale bars on the right hand side
p.addParameter('showSpanLines', true, @islogical);
p.addParameter('dataRangeFormat', '%.4g', @ischar);
p.addParameter('baseline', [], @(x) true);
p.addParameter('quick', false, @islogical);
% p.addParameter('lineStyle', '-', @ischar);
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.parse(varargin{:});

axh = newplot;

data_ci = p.Results.data_ci;
ci_alpha = p.Results.ci_alpha;

if ~iscell(data)
    % all traces share common time vector
    nTraces = size(data, 2);
    nSuperimposed =size(data, 3);
    nTime = size(data, 1);
    if isempty(tvec)
        tvec = makecol(0:nTime-1);
    end
    assert(numel(tvec) == nTime, 'Time vector must match size(data, 1)');
    
    if ~isfloat(data)
        data = single(data);
    end
    
else
    % everytraces has different time vector
    nTraces = size(data, 1);
    nSuperimposed = size(data, 2);
    emptyMask = cellfun(@isempty, data);
    data(emptyMask) = {NaN}; % prevents errors with cellfun later
    
    if isempty(data_ci)
        data_ci = cell(size(data));
        [data_ci{:}] = deal(NaN);
    end
end
   
% construct labels by traces
if isempty(p.Results.labels)
    labels = arrayfun(@(x) sprintf("%s%d", p.Results.labelPrefix, x), 1:nTraces);
elseif isnumeric(p.Results.labels)
    labels = arrayfun(@string, p.Results.labels);
else
    labels = string(p.Results.labels);
end
labels = makecol(labels);

% construct labels for each superimposed line within each trace
if isempty(p.Results.labelsSuperimposed)
    labelsLinesWithinEachTrace = arrayfun(@string, 1:nSuperimposed);
else
    labelsLinesWithinEachTrace = string(p.Results.labelsSuperimposed);
end

show_trace_within_legend = true(nTraces, 1);
show_superimposed_within_legend = true(nSuperimposed, 1);

if p.Results.legendUniqueLabelsOnly
    show_trace_within_legend(:) = false;
    [~, ia] = unique(labels);
    show_trace_within_legend(ia) = true;

    show_superimposed_within_legend(:) = false;
    [~, ia] = unique(labelsLinesWithinEachTrace);
    show_superimposed_within_legend(ia) = true;
end

if ~iscell(data)
    has_ci = ~isempty(data_ci);
    
    % invert the order of the traces so the first is plotted at the top
    if p.Results.maintainScaleSuperimposed
        rangeDims = [1 3];
    else
        rangeDims = 1;
    end
    
    % subtract the min so each group of traces has min == 0
    if has_ci
        data_lo = data_ci(:, :, :, 1);
        data_hi = data_ci(:, :, :, 2);
        
        data_lo_scaling = min(data, data_ci(:, :, :, 1), 'omitnan');
        data_hi_scaling = max(data, data_ci(:, :, :, 2), 'omitnan');
    else
        data_lo = data;
        data_hi = data;
        
        data_lo_scaling = data;
        data_hi_scaling = data;
    end

    minEachGroup = min(data_lo_scaling, [], rangeDims, 'omitnan');
    dataLowOrig = minEachGroup;
    data = data - minEachGroup;
    data_lo = data_lo - minEachGroup;
    data_hi = data_hi - minEachGroup;
    data_lo_scaling = data_lo_scaling - minEachGroup;
    data_hi_scaling = data_hi_scaling - minEachGroup;

    rangesOrig = max(data_hi_scaling, [], rangeDims, 'omitnan');

    if p.Results.normalize
        norms = rangesOrig;
        data = data ./ norms;
        data_lo = data_lo ./ norms;
        data_hi = data_hi ./ norms;
        data_lo_scaling = data_lo_scaling ./ norms;
        data_hi_scaling = data_hi_scaling ./ norms;
    else
        norms = onesvec(nTraces);
    end

    % compute the max range each row
    rangesNorm = max(data_hi_scaling, [], [1 3], 'omitnan');
    rangesNorm(isnan(rangesNorm)) = 0;

    % figure out where each trace should start
    if p.Results.evenSpacing
        % trace(k) will be offset by spacing * k
        if p.Results.intercalate
            deltas = data_hi_scaling(:, 2:end, :) * p.Results.spacingFraction - data_lo_scaling(:, 1:end-1, :);
            maxDeltas = max(deltas, [1 3], 'omitnan') / p.Results.gain; % max over time and superimposed traces
            traceOffsets = (nTraces-1:-1:0) * max(maxDeltas, [], 'omitnan');
        else
            rangesPadded = makerow(rangesNorm / p.Results.gain * (p.Results.spacingFraction));
            traceOffsets = (nTraces-1:-1:0) * max(rangesPadded, [], 'omitnan');
        end
    else
        if p.Results.intercalate
            % the offsets should be set so that pairs of points across
            % all time points should satisfy:
            %   traceN+1(t) * spacingFraction <= traceN(t) + offset
            % or:
            %   offset = max_t (traceN+1(t) * spacingFraction - traceN(t))

            deltas = data_hi_scaling(:, 2:end, :) * p.Results.spacingFraction - data_lo_scaling(:, 1:end-1, :);
            maxDeltas = max(deltas, [1 3], 'omitnan') / p.Results.gain; % max over time and superimposed traces
            cs = fliplr(cumsum(fliplr(maxDeltas)));
            traceOffsets = [cs, 0];
        else 
            rangesPadded = makerow(rangesNorm / p.Results.gain * (p.Results.spacingFraction));
            cs = fliplr(cumsum(fliplr(rangesPadded)));
            traceOffsets = [cs(2:end), 0];
        end
    end  

    % time x stack x superimposed
    data = data + traceOffsets;
    data_lo = data_lo + traceOffsets;
    data_hi = data_hi + traceOffsets;

    % expand colormap to be exactly nSuperimposed long
    [map, colorByStack] = getColormap(p.Results.colormapStacked, nTraces, p.Results.colormap, nSuperimposed);

    if nSuperimposed == 1 && ~colorByStack
        % plot simultaneously
        TrialDataUtilities.Plotting.errorshadeInterval(tvec, data_lo, data_hi, map(1, :), 'axh', axh);
        hold on;
        hLines = plot(tvec, data, '-', 'Color', map(1, :), 'MarkerFaceColor', map(1, :), 'MarkerEdgeColor', map(1, :), 'axh', axh, p.Unmatched);
    else
        % here we arrange so that all traces are stacked vertically but that
        % all positions along dim 2 are grouped together 
        matShiftCat = TensorUtils.reshapeByConcatenatingDims(data, {1 [3, 2]});
        
        if has_ci
            for iSuper = 1:nSuperimposed
                for iStack = 1:nTraces
                    if colorByStack
                        color_this = map(iStack, :);
                    else
                        color_this = map(iSuper, :);
                    end
                    TrialDataUtilities.Plotting.errorshadeInterval(tvec, data_lo(:, iStack, iSuper), data_hi(:, iStack, iSuper), color_this, 'axh', axh, 'alpha', ci_alpha);
                    hold on;
                end
            end
        end
        hLines = plot(tvec, matShiftCat, '-', p.Unmatched);
        hLines = reshape(hLines, [nSuperimposed nTraces])';
        
        if colorByStack
            for iT = 1:nTraces
                set(hLines(iT, :), 'Color', map(iT, :), 'MarkerFaceColor', map(iT, :), 'MarkerEdgeColor', map(iT, :));
                if show_trace_within_legend(iT)
                    TrialDataUtilities.Plotting.showFirstInLegend(hLines(iT, :), labels{iT});
                else
                    TrialDataUtilities.Plotting.hideInLegend(hLines(iT, :));
                end
            end
        else
            for iS = 1:nSuperimposed
                set(hLines(:, iS), 'Color', map(iS, :), 'MarkerFaceColor', map(iS, :), 'MarkerEdgeColor', map(iS, :));
                if show_superimposed_within_legend(iS)
                    TrialDataUtilities.Plotting.showFirstInLegend(hLines(:, iS), labelsLinesWithinEachTrace{iS});
                else
                    TrialDataUtilities.Plotting.hideInLegend(hLines(:, iS));
                end
            end
        end
    end

else
    % cell mode - everyone has a different time vector
    if ~iscell(tvec)
        tvec = repmat({tvec}, nTraces, nSuperimposed);
    end
    
    [map, colorByStack] = getColormap(p.Results.colormapStacked, nTraces, p.Results.colormap, nSuperimposed);
        
    % invert the order of the traces so the first is plotted at the top
%     data = flipud(data);
%     labels = flipud(labels);
%     tvec = flipud(tvec);

    if p.Results.maintainScaleSuperimposed
        % subtract the min so each group trace has min at zero
        minEachRow = nanmin(cellfun(@(data) double(nanminNanEmpty(data(:))), data), [], 2);
        dataLowOrig = minEachRow;
        minCell = num2cell(repmat(minEachRow, 1, size(data, 2)));
        cellShift = cellfun(@(data, min) data - min, data, minCell, 'UniformOutput', false);

        rangesOrig = nanmax(cellfun(@(data) double(nanmax(data(:))), cellShift), [], 2);
        if p.Results.normalize
            maxCell = num2cell(repmat(rangesOrig, 1, size(data, 2)));
            cellShift = cellfun(@(matShift, max) double(matShift ./ max), cellShift, maxCell, 'UniformOutput', false);
            norms = rangesOrig;
        else
            norms = onesvec(nTraces);
        end
    else
        % subtract the min so each trace has min at zero
        dataLowOrig = nanmin(data, [], 1);
        cellShift = cellfun(@(data) bsxfun(@minus, data, dataLowOrig), data, 'UniformOutput', false);

        if p.Results.normalize
            cellShift = cellfun(@(matShift) bsxfun(@rdivide, matShift, range(matShift, 1)), 'UniformOutput', false);
        end
    end

    % compute the max range each row
    rangesNorm = max(cellfun(@(matShift) double(nanmaxNanEmpty(matShift(:), [], 1)), cellShift), [], 2);

    rangesNorm(isnan(rangesNorm)) = 0;
    
    % can't intercalate without doing time consuming interpolation
    rangesPadded = rangesNorm * (p.Results.spacingFraction);
    cs = flipud(cumsum(flipud(rangesPadded)));
    traceOffsets = [cs(2:end); 0];

    traceOffsets = traceOffsets - min(traceOffsets);
%     traceCenters = (traceOffsets + rangesNorm / 2)';

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

            if colorByStack
                color = map(iT, :);
            else
                color = map(iS, :);
            end
            hLines{iT, iS} = plot(tvecThis, matShift, '-', 'Color', color, ...
                'MarkerFaceColor', color, 'MarkerEdgeColor', color, p.Unmatched);
            hold on;
        end
    end
end


traceCenters = (traceOffsets + rangesNorm / 2)';
traceHighs = traceOffsets + rangesNorm;
traceLows = traceOffsets;
dataHighOrig = dataLowOrig + rangesOrig;

ylim([nanmin(traceLows) nanmax(traceHighs) + eps]);

if strcmp(p.Results.showLabels, 'auto')
    showLabels = nTraces < 25;
else
    showLabels = p.Results.showLabels;
end

if isempty(p.Results.timeUnits)
    xlab = 'Time';
else
    xlab = sprintf('Time (%s)', p.Results.timeUnits);
end

if ~p.Results.quick
    au = AutoAxis(gca);

    au.xUnits = p.Results.timeUnits;
    if p.Results.timeScaleBar
        au.addAutoScaleBarX();
    else
        au.addAutoAxisX();
        au.xlabel(xlab);
    end

    if showLabels
        if ~isempty(p.Results.colormapStacked)
            colormapStacked = TrialDataUtilities.Plotting.expandWrapColormap(p.Results.colormapStacked, numel(traceLows));
        else
            colormapStacked = [0 0 0];
        end
        spans = [makerow(traceLows); makerow(traceHighs)];
        au.addLabeledSpan('y', 'span', spans, 'label', labels, 'color', colormapStacked, ...
            'rotation', p.Results.labelRotation, 'showSpanLines', p.Results.showSpanLines);
    end

    hold off;

    if p.Results.clickable
        % build line handle descriptions
    %     lineDescriptions = cell(nTraces, nSuperimposed);
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

    if strcmp(p.Results.showVerticalScaleBars, 'auto')
        showScaleBars = nTraces < 25 && p.Results.maintainScaleSuperimposed;
    else
        showScaleBars = p.Results.showVerticalScaleBars;
    end

    % add scale bars to right side of axis?
    if showScaleBars
        if ~p.Results.maintainScaleSuperimposed
            warning('Cannot show vertical scale bars when maintainScaleSuperimposed is false');
        else
            au = AutoAxis(gca);

            dataUnits = p.Results.dataUnits;
            if iscell(dataUnits)
                u = unique(dataUnits);
                commonDataUnits = numel(u) == 1;
                if commonDataUnits
                    dataUnits = u{1};
                end
            elseif isempty(dataUnits) || ischar(dataUnits)
                commonDataUnits = true;
            else
                commonDataUnits = false;
                dataUnits = '';
            end

            if ~p.Results.normalize && commonDataUnits
                % common scaling, all same units, show single auto scale bar

                au.yUnits = dataUnits;
                au.addAutoScaleBarY();

            elseif p.Results.normalize && commonDataUnits
                % need one scale bar for each of nTraces but all in the same
                % units
                % but we can label only one scale bar at the bottom and have the others be the same length
                % the actual data size of each will
                % be set by the largest amplitude signal, which has the
                % largest norm

                actualValue = TrialDataUtilities.Plotting.closestNiceNumber(min(rangesOrig), 'down');
                for iT = 1:nTraces
                    plottedValue = actualValue / norms(iT);
                    if iT == nTraces
                        % label the bottom trace
                        label = sprintf('%g %s', actualValue, dataUnits);
                    else
                        label = '';
                    end
                    au.addScaleBar('y', 'length', plottedValue, 'manualLabel', label, ...
                        'manualPositionAlongAxis', traceLows(iT), 'hideLabel', p.Results.verticalScaleBarHideLabel);
                end

            else
                % channels all have different units, show each with its own
                % labeled scale bar
                for iT = 1:nTraces
                    actualValue = TrialDataUtilities.Plotting.closestNiceNumber(rangesOrig(iT), 'down');
                    plottedValue = actualValue / norms(iT);
                    if iscell(dataUnits)
                        units = dataUnits{iT};
                    else
                        units = dataUnits;
                    end
                    label = sprintf('%g %s', actualValue, units);
                    au.addScaleBar('y', 'length', plottedValue, 'manualLabel', label, ...
                        'manualPositionAlongAxis', traceLows(iT), 'hideLabel', p.Results.verticalScaleBarHideLabel);
                end
            end

        end
    end

    if strcmp(p.Results.showDataRanges, 'auto')
        showYRanges = nTraces < 25 && p.Results.maintainScaleSuperimposed;
    else
        showYRanges = p.Results.showDataRanges;
    end

    axis tight
    
    ax = gca;
    ax.TickDir = 'out';
    if size(map, 1) > 1
        ax.ColorSpace.Colormap = map;
        ax.CLim = [1 size(map, 1)];
    end

    % show y extents as tick bridges on the right hand side
    if showYRanges && ~showScaleBars
        if ~p.Results.maintainScaleSuperimposed
            warning('Cannot show data ranges when maintainScaleSuperimposed is false');
        else
            au = AutoAxis(gca);

            dataUnits = p.Results.dataUnits;

            for iT = 1:nTraces
                if iscell(dataUnits)
                    units = dataUnits{iT};
                else
                    units = dataUnits;
                end
                labelHigh = sprintf([p.Results.dataRangeFormat, ' %s'], dataHighOrig(iT), units);
                labelLow  = sprintf([p.Results.dataRangeFormat, ' %s'], dataLowOrig(iT), units);

                au.addTickBridge('y', 'tick', [traceLows(iT) traceHighs(iT)], ...
                    'tickLabel', {labelLow; labelHigh}, 'tickAlignment', {'bottom', 'top'}, ...
                    'otherSide', true);
            end

            au.axisMarginRight = 2;
        end
    end

    hold off;
    au.update();
end

axis tight

if nargout == 0
    clear traceCenters hLines;
end

end

function r = nanmaxNanEmpty(v1, varargin)
    if isempty(v1) && numel(varargin) >= 1 && isempty(varargin{1})
        r = NaN;
    else
        r = nanmax(v1, varargin{:});
    end
end

function [map, colorByStack] = getColormap(cmapStacked, nStacked, cmapSuperimposed, nSuperimposed)
    if ~isempty(cmapSuperimposed)
        cmapFn = cmapSuperimposed;
        colorByStack = false;
        N = nSuperimposed;
        map = TrialDataUtilities.Plotting.expandWrapColormap(cmapFn, N);
    elseif ~isempty(cmapStacked)
        cmapFn = cmapStacked;
        colorByStack = true;
        N = nStacked;
        map = TrialDataUtilities.Plotting.expandWrapColormap(cmapFn, N);
    else
        % no colormap specified
        if nSuperimposed == 1
            map = [0 0 0];
            colorByStack = false;
        else
            map =  TrialDataUtilities.Plotting.expandWrapColormap(@TrialDataUtilities.Colormaps.linspecer, nSuperimposed);
            colorByStack = false;
        end
    end
end
