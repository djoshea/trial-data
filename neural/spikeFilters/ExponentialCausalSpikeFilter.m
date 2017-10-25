classdef ExponentialCausalSpikeFilter < ConvolutionSpikeFilter 
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        tauMs
        filterLength % length of filter in multiples of tau
    end

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms tau', sf.tauMs);
        end
    end
    
    methods
        function sf = ExponentialCausalSpikeFilter(varargin)
            p = inputParser;
            p.addParamValue('tauMs', 25, @isscalar);
            p.addParamValue('filterLengthTaus', 4, @isscalar);
            p.parse(varargin{:});

            sf.binAlignmentMode = BinAlignmentMode.Causal;
            sf.tauMs = p.Results.tauMs;
            sf.filterLength = p.Results.filterLengthTaus * sf.tauMs;
        end
               
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt indZero] = getFilter(sf)
            t = 0:sf.filterLength;
            filt = exp(-t/sf.tauMs);
            indZero = 1;
        end
    end

end
