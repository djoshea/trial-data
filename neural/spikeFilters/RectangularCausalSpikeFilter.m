classdef RectangularCausalSpikeFilter < ConvolutionSpikeFilter 

    properties
        widthMs % width of the rectangular filter. This may differ from the raw bin width which will be convolved over
    end

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms rect', sf.widthMs);
        end
    end
    
    methods
        function sf = RectangularCausalSpikeFilter(varargin)
            % params:
            %   widthMs : [default 25]
            p = inputParser;
            p.addOptional('widthMs', 25, @isscalar);
            p.parse(varargin{:});

            sf.widthMs = p.Results.widthMs;
            sf.binAlignmentMode = BinAlignmentMode.Causal;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            filt = onesvec(floor(sf.widthMs / sf.binWidthMs));
            indZero = 1; % this makes it causal
        end
        
        function v = getBinAlignmentMode(sf, v)
            v = BinAlignmentMode.Causal;
        end
    end

end
