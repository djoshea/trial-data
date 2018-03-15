classdef GaussianSpikeFilter < ConvolutionSpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        % the std deviation of the gaussian around the peak in ms
        sigma
        
        % the filter will extend sigma*halfWidthSigma into the past and
        % future to cover the Gaussian, unless truncated by the 
        halfWidthSigmas = 3;
        
        % the center of the Gaussian will be located this many ms into the
        % future
        delayPeak 
        
        % truncate the filter from looking further into the future at (ms
        % in the future), acausal taps
        truncateFuture

        % truncate the filter from looking further into the past (ms in the
        % past), causal taps
        truncatePast
    end

    methods
        function sf = GaussianSpikeFilter(varargin)
            p = inputParser;
            p.addParameter('sigma', 20, @isscalar); % in ms
            p.addParameter('halfWidthSigmas', 3, @isscalar);
            % center the peak of the Gaussian impulse this many ms in the future
            p.addParameter('delayPeak', 0, @isscalar);
            p.addParameter('truncateFuture', Inf, @isscalar);
            p.addParameter('truncatePast', Inf, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            sf = sf@ConvolutionSpikeFilter(p.Unmatched);
            sf.sigma = p.Results.sigma;
            sf.halfWidthSigmas = p.Results.halfWidthSigmas;
            sf.delayPeak = p.Results.delayPeak;
            sf.truncateFuture = p.Results.truncateFuture;
            sf.truncatePast = p.Results.truncatePast;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if getFilterIndZero > 1
        function [filt, indZero] = getFilter(sf)
            sigmaMultiple = sf.halfWidthSigmas;
            
            % we care about 3 sigma in the future from the delayPeak 
            % unless we're truncating beyond a certain point in the future
            % regardless we must overlap with 0
            tMin = min(0, max(ceil(-sigmaMultiple*sf.sigma) + sf.delayPeak, -sf.truncateFuture));

            % past is positive time (causal taps)
            % we care about 3 sigma in the past from the delayPeak
            % unless we're truncating beyond a certain point in the past
            % regardless we must overlap with 0
            tMax = max(0, min(floor(sigmaMultiple*sf.sigma) + sf.delayPeak, sf.truncatePast));

            % compute the gaussian
            t = TrialDataUtilities.Data.linspaceIntercept(tMin, sf.binWidthMs, tMax, 0);
            filt = exp(-(t-sf.delayPeak).^2 / (2*sf.sigma^2));
            
            % unnecessary unless our distribution doesn't have support over
            % zero, which is inadvisable
            filt(t < -sf.truncateFuture | t > sf.truncatePast) = 0;

            indZero = find(t == 0);
        end
    end
    
    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f sigma, center at %.0f', sf.sigma, sf.delayPeak);
        end
    end

end
