classdef SpikeFilter % < handle & matlab.mixin.Copyable
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    % Dependent properties are provided as a convenience, override the underlying methods
    properties(Dependent)
        preWindow
        postWindow
        padWindow % includes pre / post window plus binWidth time for the edge bins
        isCausal
    end

    methods(Abstract)
        % returns the [pre post] window in ms, outside of which the filter value is treated as 0
        windowMs = getFilterWindow(sf);
           
        % returns a function that maps time in ms to a filter amplitude at that time offset. MUST be vectorized
        % and return output the same size and class (typically single) as input
        fn = getFilterFn(sf);
    end
        
    methods        
        function t = get.preWindow(sf)
            t = sf.getPadWindow();
            t = t(1);
        end

        function t = get.postWindow(sf)
            t = sf.getPadWindow();
            t = t(2);
        end
        
        function w = get.padWindow(sf)
            w = sf.getPadWindow();
        end
                         
        function tf = get.isCausal(sf)
            tf = sf.getIsCausal();
        end
        function w = getPadWindow(sf)
            % pre and post window are already in ms to accommodate the filter
            % then we accommodate both
            w = sf.getFilterWindow();
        end

        function checkSettingsOkay(sf)
            % doesn't make sense to sample more finely than the
            % spikeBinWidth used
            
        end

        % spikeCell is nTrains x 1 cell array of time points
        %
        % tWindowPerTrial is nTrials x 2 matrix of tMin and tMax to grab
        % per trial, excluding any padding which should be accommodated for
        % when grabbing spikeCell
        %
        % multiplierToSpikesPerSec is the factor by which rates in the
        % spikeCell time units should be multiplied to get sec (e.g. 1000
        % for spikeCell in ms)
        % 
        % this subclass method should respect sf.timeDelta,
        % sf.binAlignmentMode, and sf.resampleMethod. We could do this work
        % here, but not doing so leaves flexibility to the subclass
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
            
            % normalization is critical to maintain filtered signal as a
            % spike rate
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
    
    methods(Access=protected) % Subclasses may wish to override these 
        % return a string describing this filter's particular parameters,
        % do not include the classname as this will be added automatically
        function str = subclassGetDescription(sf) %#ok<MANU>
            str = '';
        end
    end

    methods 
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
            yl(1) = min(0, yl(1));
            yl(2) = max(0, yl(2));
            
            ylim(yl);
            ax = AutoAxis.replace();
            
            ax.addAnchor(AutoAxis.AnchorInfo(h, AutoAxis.PositionType.Bottom, gca, AutoAxis.PositionType.Top, ax.axisPaddingTop));
            ax.update();
        end
        
        function str = getDescription(sf)
            % build description using subclass implemented method including
            % causality and classname
            subclassStr = sf.subclassGetDescription();
            classname = class(sf);
            if sf.isCausal()
                causalStr = 'Causal';
            else
                causalStr = 'Acausal';
            end
            if isempty(subclassStr)
                str = sprintf('%s %s', causalStr, classname);
            else
                str = sprintf('%s %s (%s)', causalStr, classname, subclassStr);
            end
        end
        
        function str = char(sf)
            str = sf.getDescription();
        end
    end

    methods
        function sf = SpikeFilter()
        end
       
        function [rateCell, timeCell] = filterSpikeTrainsWindowByTrial(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, varargin)
            % filters each trial individually, using a cell array to return
            % filtered rates and timeCell 
            sf.checkOkay();
            tWindowMat = [makecol(tMinByTrial), makecol(tMaxByTrial)];
            [rateCell, timeCell] = sf.subclassFilterSpikeTrains(spikeCell, tWindowMat, multiplierToSpikesPerSec, varargin{:});
        end
        
        function [rates, tvec] = filterSpikeTrainsWindowByTrialAsMatrix(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, varargin)
            % filters each trial individually and then embeds each filtered
            % trace in a nTrials x nTime x nUnits matrix, where missing samples are left as NaN
            % before and after each trial. tvec is the time vector that
            % indicates time along the columns.=
            p = inputParser;
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('interpolateMethod', 'linear', @ischar);
            p.addParameter('showProgress', true, @islogical);    
            p.addParameter('tMinByTrialExcludingPadding', [], @isvector);
            p.addParameter('tMaxByTrialExcludingPadding', [], @isvector);
            p.parse(varargin{:});
            
            % calls checkOkay
            % don't do resampling within since we'll do it all at once
            % below, by setting useTimeDelta false
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, 'useTimeDelta', false);
            
            % convert to matrix
            [rates, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(rateCell, timeCell, ...
                'assumeUniformSampling', true, ...
                'fixNonmonotonicTimes', false, ...
                'resampleMethod', sf.resampleMethod, ...
                'origDelta', sf.binWidthMs, ... % saves time if this is provided rather than inferred, requires 'useTimeDelta' false above
                'timeDelta', sf.timeDelta, ...
                'tMinExcludingPadding', p.Results.tMinByTrialExcludingPadding, ...
                'tMaxExcludingPadding', p.Results.tMaxByTrialExcludingPadding); 
        end
    end
    
    methods(Static)
        function sf = getDefaultFilter()
            sf = GaussianSpikeFilter();
        end
    end
end
