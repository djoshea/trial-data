classdef NonOverlappingSpikeBinFilter < ConvolutionSpikeFilter 
% when using this class, be sure to set .timeDelta == .binWidthMs so that
% the time windows are sampled evenly. Otherwise an error will be thrown.

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms bins', sf.binWidthMs);
        end
    end
    
    methods
        function sf = NonOverlappingSpikeBinFilter(varargin)
            p = inputParser;
            p.addOptional('binWidthMs', 1, @isscalar);
            p.parse(varargin{:});

            sf.binWidthMs = p.Results.binWidthMs;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            filt = 1;
            indZero = 1;
        end
    end
    
    methods(Access=protected) % superclass overrides
        function checkTimeDeltaOkay(sf, timeDelta)
            assert(timeDelta == sf.binWidthMs, 'timeDelta must match binWidthMs in order to get non-overlapping bins');
        end
    end

end
