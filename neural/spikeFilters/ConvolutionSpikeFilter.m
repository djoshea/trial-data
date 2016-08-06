classdef ConvolutionSpikeFilter < SpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        % bin width that spikes will be binned into before filtering.
        % this is different from the sampling rate of the filtered firing
        % rate, which is determined by the time delta parameter passed
        % into .filterSpikeTrains by the caller
        binWidthMs = 1;
        binAlignmentMode = SpikeBinAlignmentMode.Acausal;
    end

    properties(Dependent)
        filter % cache the filter here so as to avoid recomputation
        indZero
    end

    methods 
        function f = get.filter(sf)
            f = sf.getFilter();
        end

        function ind = get.indZero(sf)
            [~, ind] = sf.getFilter();
        end
    end

    methods(Abstract)
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 0
        % doesn't need to be normalized, we'll take care of that
        [filt, indZero] = getFilter(sf);
    end
    
    methods
        function sf = ConvolutionSpikeFilter(varargin)
            p = inputParser;
            p.addParamValue('binWidthMs', 1, @isscalar);
            p.addParamValue('binAlignmentMode', SpikeBinAlignmentMode.Acausal, @(x) isa(x, 'SpikeBinAlignmentMode'));
            p.parse(varargin{:});

            sf.binWidthMs = p.Results.binWidthMs;
            sf.binAlignmentMode = p.Results.binAlignmentMode;
        end
        
        function plotFilter(sf)
            cla;
            [filt, indZero] = sf.getFilter(); %#ok<PROP>
            filt = filt ./ sum(filt);
            tvec = (1:numel(filt))' - indZero; %#ok<PROP>
            plot(tvec, filt, '.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
            xlabel('Time (ms)');
            ylabel('Impulse Response');
            
            hold on;
            TrialDataUtilities.Plotting.verticalLine(0, 'Color', [0.7 0.7 0.7]);
            h = gobjects(0, 1);
            if sf.getPreWindow() > 0
                h(end+1) = text(0, 0, '  Causal', 'HorizontalAlignment', 'left', 'BackgroundColor', 'none', 'VerticalAlignment', 'bottom', 'XLimInclude', 'on');
            end
            if sf.getPostWindow() > 0
                h(end+1) = text(0, 0, 'Acausal  ', 'HorizontalAlignment', 'right', 'BackgroundColor', 'none', 'VerticalAlignment', 'bottom', 'XLimInclude', 'on');
           end
            hold off;
            ax = AutoAxis.replace();
            ax.addAnchor(AutoAxis.AnchorInfo(h, AutoAxis.PositionType.Bottom, gca, AutoAxis.PositionType.Bottom, 0.1));
            ax.update();
        end
    end

    methods(Access=protected)     
        function checkTimeDeltaOkay(sf, timeDelta)
            % doesn't make sense to sample more finely than the
            % spikeBinWidth used
            assert(timeDelta >= sf.binWidthMs, 'TimeDelta is shorter than spike filter binWidthMs. Filter output should not be sampled more finely than the spike bin resolution');
        end
        
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPreWindow(sf)
            filtSize = length(sf.filter);
            % this gives us the right number of ms for the spike bins to
            % the left of zero, as well as the extra ms needed for the t=0
            % bin
            t = (filtSize - sf.indZero)*sf.binWidthMs - sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.binWidthMs);
        end

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPostWindow(sf)
            % this gives us the right number of ms for the spike bins to
            % the right of zero, as well as the extra ms needed for the t=0
            % bin
            t = (sf.indZero - 1)*sf.binWidthMs + sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.binWidthMs); 
        end
        
        function tf = getIsCausal(sf) % allows subclasses to override
            tf = sf.getPostWindow() <= 0 && sf.binAlignmentMode == SpikeBinAlignmentMode.Causal;
        end

        % spikeCell is nTrains x 1 cell array of time points which will include 
        %   times in the preceding and postceding window
        function [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowByTrial, multiplierToSpikesPerSec)
            % build filter
            
            if isempty(spikeCell)
                rateCell = {};
                timeCell = {};
                return;
            end
            
            filt = sf.filter;
            % normalization is critical
            filt = filt ./ sum(filt);
            
            tPadPre = sf.preWindow;
            tPadPost = sf.postWindow;

            tMinByTrial = floor(tWindowByTrial(:, 1));
            tMaxByTrial = floor(tWindowByTrial(:, 2));
            
            % get time vector for bins including pre and post padding to accomodate filter
            % this depends on both the bin width and the bin alignment mode
            [~, tbinsForHistcByTrial] = SpikeBinAlignmentMode.generateMultipleBinnedTimeVectors(...
                    tMinByTrial-tPadPre, tMaxByTrial+tPadPost, sf.binWidthMs, sf.binAlignmentMode);
                
            % timeCell contains time vector without padding bins
            binOffsetStart = sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.binWidthMs);
            binOffsetStop = sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.binWidthMs);
            timeCell = SpikeBinAlignmentMode.generateMultipleBinnedTimeVectors(...
                tMinByTrial+binOffsetStart, tMaxByTrial+binOffsetStop, sf.binWidthMs, sf.binAlignmentMode);
            
            % filter via valid convolution, which automatically removes the padding
            nTrials = length(spikeCell);
            rateCell = cellvec(nTrials);
            for i = 1:nTrials
                if ~isempty(spikeCell{i})
                    countsPad = histc(spikeCell{i}, tbinsForHistcByTrial{i});
                    rateCell{i} = makecol(conv(countsPad(1:end-1), filt, 'valid') * multiplierToSpikesPerSec / sf.binWidthMs);
                else
                    rateCell{i} = zeros(size(timeCell{i}));
                end
            end
        end
    end
    
end
