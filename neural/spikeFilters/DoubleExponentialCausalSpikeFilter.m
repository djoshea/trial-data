classdef DoubleExponentialCausalSpikeFilter < ConvolutionSpikeFilter 
% based on Thompson KG, Hanes DP, Bichot NP, Schall JD (1996) Perceptual and mo- tor processing stages identified in the activity of macaque frontal eye field neurons during visual search. J Neurophysiol 76:4040?4055
% mimics a post-synaptic current

    properties
        tauGrowth = 1;
        tauDecay = 20;
        filterLength = 80;
    end

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms growth, %.0f ms decay', sf.tauGrowth, sf.tauDecay);
        end
    end
    
    methods
        function sf = DoubleExponentialCausalSpikeFilter(varargin)
            sf.binAlignmentMode = BinAlignmentMode.Causal;
        end

        function v = getBinAlignmentMode(~, ~)
           v = BinAlignmentMode.Causal;
        end
               
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            t = 0:sf.filterLength;
            filt = (1-exp(-t/sf.tauGrowth)).*exp(-t/sf.tauDecay);
            indZero = 1;
        end
    end

end
