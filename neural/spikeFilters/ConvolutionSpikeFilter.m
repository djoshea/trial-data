classdef ConvolutionSpikeFilter < SpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        % bin width that spikes may be binned into before filtering.
        % this is different from the sampling rate of the filtered firing
        % rate, which is determined by timeDelta.
        % This setting also determines the width of the bins specified by
        % the convolution filter
        binWidthMs = 1;
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
            p = inputParser();
            p.addOptional('binWidthMs', 1, @(x) isnumeric(x) && isscalar(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            sf = sf@SpikeFilter(p.Unmatched);
            sf.binWidthMs = p.Results.binWidthMs;
        end
        
        function plotFilter(sf)
            cla;
            [filt, indZero] = sf.getFilter(); %#ok<PROP>
            filt = filt ./ sum(filt);
            tvec = ((1:numel(filt))' - indZero) * sf.binWidthMs; %#ok<PROP>
            TrialDataUtilities.Plotting.verticalLine(0, 'Color', [0.7 0.7 0.7]);
            hold on;
            if true || isscalar(filt)
                stem(tvec, filt, 'filled', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
%             else
%                 plot(tvec, filt, '.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
            end
            xlabel('Time (ms)');
            ylabel('Impulse Response');
            
            h = gobjects(0, 1);
            if sf.getPreWindow() > 0
                h(end+1) = text(0, 0, '  Causal', 'HorizontalAlignment', 'left', 'BackgroundColor', 'none', 'VerticalAlignment', 'bottom', 'XLimInclude', 'on');
            end
            if sf.getPostWindow() > 0
                h(end+1) = text(0, 0, 'Acausal  ', 'HorizontalAlignment', 'right', 'BackgroundColor', 'none', 'VerticalAlignment', 'bottom', 'XLimInclude', 'on');
            end
            TrialDataUtilities.Plotting.verticalLine(0, 'Color', 'r');
            hold off;
            xlim([tvec(1) - sf.binWidthMs tvec(end) + sf.binWidthMs]);
            
            del = max(filt) - min(filt);
            yl = [min(filt) - del*0.05, max(filt) + del*0.05];
            ylim(yl);
            ax = AutoAxis.replace();
            
            ax.addAnchor(AutoAxis.AnchorInfo(h, AutoAxis.PositionType.Bottom, gca, AutoAxis.PositionType.Top, ax.axisPaddingTop));
            ax.update();
        end
    end

    % Time binning here is subtle. There are a few things to consider -
    % first, we want to provide a time vector that begins at start and
    % ends at stop, factoring in timeDelta (the ultimate rate
    % we are sampling at). So if time delta is 20, and start is -300,
    % we're going to be including a sample that effectively spans
    % -310:-290 (or -320:300 if Causal binning). Then we have to factor
    % in the spike filter itself, which will need extra bins to the
    % left and right in order to provide a sample at its center (e.g. a
    % 60 ms wide filter needs 30 ms pre and post). Lastly, we need to
    % factor in the time window required by spikeBinMs for each bin
    % that goes into the convolution step. 
    %
    % The spike data padding will be taken care of by the caller, and
    % will rely on getPadWindow().
    
    methods(Access=protected)      
        function w = getPadWindow(sf)
            % pre and post window are already in ms to accommodate the filter
            % then we accommodate both 
            bin = max(sf.binWidthMs, sf.timeDelta);
            w = [sf.preWindow - sf.binAlignmentMode.getBinStartOffsetForBinWidth(bin), ...
                sf.postWindow + sf.binAlignmentMode.getBinStopOffsetForBinWidth(bin)];
        end
        
        function checkSettingsOkay(sf)
            % doesn't make sense to sample more finely than the
            % spikeBinWidth used
            assert(sf.timeDelta >= sf.binWidthMs, 'timeDelta is shorter than spike filter binWidthMs. Filter output should not be sampled more finely than the spike bin resolution');
        end
        
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPreWindow(sf)
            % this gives us the right number of ms to include for the spike bins 
            % extending in the causal direction, which are the inds right
            % of zero in the filter (which is the impulse response).
            filtSize = length(sf.filter);
            t = (filtSize - sf.indZero) * sf.binWidthMs;
        end

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPostWindow(sf)
            % this gives us the right number of ms to include for the spike bins 
            % extending in the acausal direction, which are the inds left
            % of zero in the filter (which is the impulse response).
 
            t = (sf.indZero - 1)*sf.binWidthMs;
        end
        
        function tf = getIsCausal(sf) % allows subclasses to override
            tf = sf.getPostWindow() <= 0 && sf.binAlignmentMode == BinAlignmentMode.Causal;
        end

        % spikeCell is nTrains x 1 cell array of time points which will include 
        % times in the preceding and postceding padding window.
        
        function [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowByTrial, multiplierToSpikesPerSec, varargin)
            p = inputParser();
            p.addParameter('useTimeDelta', true, @islogical);
            p.parse(varargin{:});
            
            if p.Results.useTimeDelta
                timeDelta = sf.timeDelta;
            else
                % resampling will be done by caller later (typically on all
                % trials at once)
                timeDelta = sf.binWidthMs;
            end
            
            if isempty(spikeCell)
                rateCell = cell(size(spikeCell));
                timeCell = cell(size(spikeCell, 1), 0);
                return;
            end
            
            filt = sf.filter;
            % normalization is critical
            filt = filt ./ sum(filt);
            
             % get time vector for bins including pre and post padding to accomodate filter
            % this depends on both the bin width, time delta, and the bin
            % alignment mode and is computed in getPadWindow
            
            % Example: causal binning, binWidthMs = 1, timeDelta=20. To get a sample
            % at 0 we need binning to start at -19 (tMin = -20)
            % If binWidthMs = 20 and timeDelta = 1, to get a sample at 0 we
            % need binning to start at 0, which already spans -20:0.
            
            if timeDelta > sf.binWidthMs
            	window = [-sf.preWindow + sf.binAlignmentMode.getBinStartOffsetForBinWidth(timeDelta), ...
                         sf.postWindow + sf.binAlignmentMode.getBinStopOffsetForBinWidth(timeDelta)];
            else
                window = [-sf.preWindow sf.postWindow];
            end
            
            tMinByTrial = floor(tWindowByTrial(:, 1));
            tMaxByTrial = floor(tWindowByTrial(:, 2));
            
           [timeLabels, tbinsForHistcByTrial] = sf.binAlignmentMode.generateMultipleBinnedTimeVectors(...
                    tMinByTrial+window(1), tMaxByTrial+window(2), sf.binWidthMs); %#ok<ASGLU> % timeLabels is useful for debugging
                
            % then construct timeCell, which contains time vector without the extra samples
            % needed for the convolution, but WTIH the potential extra
            % samples needed for the input to resampling from binWidthMs to
            % timeDelta.
            
            
            if timeDelta > sf.binWidthMs
                % we need to set the bin limits wider to facilitate
                % resampling, as per the first example above.
                tPre = sf.binAlignmentMode.getBinStartOffsetForBinWidth(timeDelta);
                tPost = sf.binAlignmentMode.getBinStopOffsetForBinWidth(timeDelta);
            else
                tPre = 0;
                tPost = 0;
            end            
            
            % note that we don't use sf.timeDelta as the fourth argument
            % here because we're not binning to timeDelta below. We don't
            % convert to timeDelta until the resampling step at the end,
            % which will take care of the adjusted time limits then. The
            % tPre and tPost ensure that when the time limits are adjusted
            % they end up ultimately at tMinByTrial and tMaxByTrial.
            timeCell = sf.binAlignmentMode.generateMultipleBinnedTimeVectors(...
                tMinByTrial + tPre, tMaxByTrial + tPost, sf.binWidthMs, sf.binWidthMs);
            
            % filter via valid-region convolution, which automatically removes the padding
            rateCell = cellvec(size(spikeCell, 1));
            nTrials = size(rateCell, 1);
            nUnits = size(spikeCell(:, :), 2);
            for i = 1:nTrials
                nTimeThis = numel(timeCell{i});
                rateCell{i} = zeros(nTimeThis, nUnits); 
                for j = 1:nUnits
                    if ~isempty(spikeCell{i, j})
                        countsPad = histc(spikeCell{i, j}, tbinsForHistcByTrial{i}); % we drop the last bin from this since histc returns an last edge bin == last edge which we don't need
                        rateCell{i}(:, j) = makecol(conv(countsPad(1:end-1), filt, 'valid') * multiplierToSpikesPerSec / sf.binWidthMs);
                    elseif ~isnan(tMinByTrial(i))
                        % put in zeros
                        rateCell{i}(:, j) = zeros(nTimeThis, 1);
                    end
                end
            end
            
            % do the resampling to timeDelta bins
            if p.Results.useTimeDelta && ~TrialDataUtilities.Data.isequaltol(timeDelta, sf.binWidthMs)
                [rateCell, timeCell] = TrialDataUtilities.Data.resampleDataCellInTime(rateCell, timeCell, 'timeDelta', sf.timeDelta, ...
                    'timeReference', 0, 'binAlignmentMode', sf.binAlignmentMode, ...
                    'resampleMethod', sf.resampleMethod, 'uniformlySampled', true, ...
                    'tMinExcludingPadding', tMinByTrial, 'tMaxExcludingPadding', tMaxByTrial);
            end
        end
    end
    
end
