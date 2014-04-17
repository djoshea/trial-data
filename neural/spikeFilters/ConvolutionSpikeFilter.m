classdef ConvolutionSpikeFilter < SpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

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
        [filt, indZero] = getFilter(sf);
    end

    methods(Access=protected)
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPreWindow(sf)
            filtSize = length(sf.filter);
            t = filtSize - sf.indZero;
        end

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPostWindow(sf)
            t = sf.indZero - 1; 
        end

        % spikeCell is nTrains x 1 cell array of time points which will include 
        %   times in the preceding and postceding window
        function [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowByTrial, multiplierToSpikesPerSec)
            % build filter
            filt = sf.filter;
            filt = filt ./ sum(filt);
            
            nPre = sf.preWindow;
            nPost = sf.postWindow;

            tMinByTrial = floor(tWindowByTrial(:, 1));
            tMaxByTrial = floor(tWindowByTrial(:, 2));

            % filter via valid convolution, which automatically removes the padding
            nTrials = length(spikeCell);
            rateCell = cellvec(nTrials);
            timeCell = cellvec(nTrials);
            
            for i = 1:nTrials
                % get time vector including pre and post padding to accomodate filter
                tvecPad = tMinByTrial(i)-nPre:tMaxByTrial(i)+nPost+1;
                
                % timeCell contains time vector without padding
                timeCell{i} = (tMinByTrial(i):tMaxByTrial(i))';
                
                if ~isempty(spikeCell{i})
                    countsPad = histc(spikeCell{i}, tvecPad);
                    rateCell{i} = makecol(conv(countsPad(1:end-1), filt, 'valid') * multiplierToSpikesPerSec);
                else
                    rateCell{i} = zeros(size(timeCell{i}));
                end
            end
        end
    end

end
