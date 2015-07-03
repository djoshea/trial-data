classdef RectangularCausalSpikeFilter < ConvolutionSpikeFilter 

    properties
        widthMs
    end

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms rect', sf.widthMs);
        end
    end
    
    methods
        function sf = RectangularCausalSpikeFilter(varargin)
            p = inputParser;
            p.addOptional('widthMs', 25, @isscalar);
            p.parse(varargin{:});

            sf.widthMs = p.Results.widthMs;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            filt = onesvec(floor(sf.widthMs / sf.binWidthMs));
            indZero = 1; % this makes it causal
        end
    end

end
