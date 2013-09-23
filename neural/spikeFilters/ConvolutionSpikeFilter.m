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
        [filt indZero] = getFilter(sf);
    end

    methods(Access=protected)
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPreWindowMs(sf)
            filtSize = length(sf.filter);
            t = filtSize - sf.indZero;
        end

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        function t = getPostWindowMs(sf)
            t = sf.indZero - 1; 
        end

        % spikeCell is nTrains x 1 cell array of time points which will include 
        %   times in the preceding and postceding window
        function rates = subclassFilterSpikeTrains(sf, spikeCell, tWindow)
            % build filter
            filt = sf.filter;
            filt = filt ./ sum(filt);

            tMin = tWindow(1);
            tMax = tWindow(2);
            tvec = tMin:tMax-1;
            nTime = length(tvec);

            % get time vector including pre and post padding to accomodate filter
            nPre = sf.preWindowMs;
            nPost = sf.postWindowMs;
            tvecPad = tMin-nPre:tMax+nPost;
            nTimePad = length(tvecPad);
            
            % filter via valid convolution, which automatically removes the padding
            nTrains = length(spikeCell);
            rates = zeros(nTrains, nTime);
            for i = 1:nTrains
                if ~isempty(spikeCell{i})
                    countsPad = histc(spikeCell{i}, tvecPad);
                    rates(i, :) = conv(countsPad(1:end-1), filt, 'valid');
                end
            end

            % correct from spikes/ms to spikes/s
            rates = rates * 1000;
        end
    end

end
