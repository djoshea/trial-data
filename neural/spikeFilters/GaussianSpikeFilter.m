classdef GaussianSpikeFilter < ConvolutionSpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        % the std deviation of the gaussian around the peak
        sigma
        
        % the center of the Gaussian will be located this many ms into the
        % future
        delayPeak 
        
        % truncate the filter from looking further into the future at (ms in the future):
        truncateFuture

        % truncate the filter from looking further into the past (ms in the past) 
        truncatePast
    end

    methods
        function sf = GaussianSpikeFilter(varargin)
            p = inputParser;
            p.addParamValue('sigma', 20, @isscalar);
            % center the peak of the Gaussian impulse this many ms in the future
            p.addParamValue('delayPeak', 0, @isscalar);
            p.addParamValue('truncateFuture', Inf, @isscalar);
            p.addParamValue('truncatePast', Inf, @isscalar);
            p.parse(varargin{:});

            sf.sigma = p.Results.sigma;
            sf.delayPeak = p.Results.delayPeak;
            sf.truncateFuture = p.Results.truncateFuture;
            sf.truncatePast = p.Results.truncatePast;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if getFilterIndZero > 1
        function [filt, indZero] = getFilter(sf)
            sigmaMultiple = 3;
            % future is negative time
            % we care about 3 sigma in the future from the delayPeak 
            % unless we're truncating beyond a certain point in the future
            % regardless we must overlap with 0
            tMin = min(0, max(ceil(-sigmaMultiple*sf.sigma + 1/2) + sf.delayPeak, -sf.truncateFuture));
                
            % past is positive time
            % we care about 3 sigma in the past from the delayPeak
            % unless we're truncating beyond a certain point in the past
            % regardless we must overlap with 0
            tMax = max(0, min(floor(sigmaMultiple*sf.sigma) + sf.delayPeak, sf.truncatePast));

            % compute the gaussian
            t = tMin:tMax;
            filt = exp(-(t-sf.delayPeak).^2 / (2*sf.sigma^2));
            filt(t < -sf.truncateFuture | t > sf.truncatePast) = 0;

            indZero = find(t == 0);
        end
    end
    
    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms sigma, center at %.0f', sf.sigma, sf.delayPeak);
        end
    end

end
