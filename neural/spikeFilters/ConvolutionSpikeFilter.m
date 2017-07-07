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
            p.addParameter('binWidthMs', 1, @isscalar);
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
        function w = getPadWindow(sf)
            w = [sf.preWindow - sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.binWidthMs), ...
                sf.postWindow + sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.binWidthMs)];
        end
        
        function checkSettingsOkay(sf)
            % doesn't make sense to sample more finely than the
            % spikeBinWidth used
            assert(sf.timeDelta >= sf.binWidthMs, 'timeDelta is shorter than spike filter binWidthMs. Filter output should not be sampled more finely than the spike bin resolution');
        end
        
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPreWindow(sf)
            % this gives us the right number of ms for the spike bins to
            % the left of zero, as well as the extra ms needed for the t=0
            % bin
            %t = (sf.indZero - 1)*sf.binWidthMs + sf.binAlignmentMode.getBinStopOffsetForBinWidth(max(sf.binWidthMs, sf.timeDelta)); 
            t = (sf.indZero - 1)*sf.binWidthMs;
        end

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPostWindow(sf)
            % this gives us the right number of ms for the spike bins to
            % the right of zero. We no longer need to include the extra ms needed for the t=0
            % bin
%             t = (sf.indZero - 1)*sf.binWidthMs + sf.binAlignmentMode.getBinStopOffsetForBinWidth(max(sf.binWidthMs, sf.timeDelta)); 
            t = (sf.indZero - 1)*sf.binWidthMs; 
        end
        
        function tf = getIsCausal(sf) % allows subclasses to override
            tf = sf.getPostWindow() <= 0 && sf.binAlignmentMode == BinAlignmentMode.Causal;
        end

        % spikeCell is nTrains x 1 cell array of time points which will include 
        % times in the preceding and postceding padding window.
        % padding here includes padding due to the spike filter (additional
        % time bins pre and post), as well as additional times included to
        % facilitate the binning itself with timeDelta. 
        function [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowByTrial, multiplierToSpikesPerSec, varargin)
            if isempty(spikeCell)
                rateCell = cell(size(spikeCell));
                timeCell = cell(size(spikeCell, 1), 0);
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
            [timeLabels, tbinsForHistcByTrial] = sf.binAlignmentMode.generateMultipleBinnedTimeVectors(...
                    tMinByTrial-tPadPre, tMaxByTrial+tPadPost, sf.binWidthMs); %#ok<ASGLU> % timeLabels is useful for debugging
                
            % timeCell contains time vector without padding bins
%             binOffsetStart = sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.binWidthMs);
%             binOffsetStop = sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.binWidthMs);
            
%             timeDeltaOffsetStart = sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.timeDelta);
%             timeDeltaOffsetStop = sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.timeDelta);
%             timeCell = BinAlignmentMode.generateMultipleBinnedTimeVectors(...
%                 tMinByTrial+binOffsetStart, tMaxByTrial+binOffsetStop, sf.binWidthMs, sf.binAlignmentMode);
%             timeCell = sf.binAlignmentMode.generateMultipleBinnedTimeVectors(...
%                 tMinByTrial+timeDeltaOffsetStart, tMaxByTrial+timeDeltaOffsetStop, sf.binWidthMs);
            timeCell = sf.binAlignmentMode.generateMultipleBinnedTimeVectors(...
                tMinByTrial, tMaxByTrial, sf.binWidthMs, sf.binWidthMs);
            
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
            if ~TrialDataUtilities.Data.isequaltol(sf.timeDelta, sf.binWidthMs)
                [rateCell, timeCell] = TrialDataUtilities.Data.resampleDataCellInTime(rateCell, timeCell, 'timeDelta', sf.timeDelta, ...
                    'timeReference', 0, 'binAlignmentMode', sf.binAlignmentMode, ...
                    'resampleMethod', sf.resampleMethod, 'uniformlySampled', true, ...
                    'tMinExcludingPadding', tMinByTrial, 'tMaxExcludingPadding', tMaxByTrial);
            end
        end
    end
    
end
