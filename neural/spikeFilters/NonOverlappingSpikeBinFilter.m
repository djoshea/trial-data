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
            % args: 'binWidthMs', #, binAlignmentMode, BinAlignmentMode.Acausal
            sf = sf@ConvolutionSpikeFilter(varargin{:});
            sf.timeDelta = sf.binWidthMs;
        end
        
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            filt = 1;
            indZero = 1;
        end
        
        function sf = postSetTimeDelta(sf)
            sf.binWidthMs = sf.timeDelta;
        end
        
        function checkOkay(sf) % superclass overrides
            assert(sf.timeDelta == sf.binWidthMs, 'timeDelta must match binWidthMs in order to get non-overlapping bins');
        end
    end
 
end
