classdef NonOverlappingSpikeBinFilter < ConvolutionSpikeFilter 
% when using this class, timeDelta will automatically return  binWidthMs so that
% the time windows are sampled evenly.

    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f ms bins', sf.binWidthMs);
        end
    end
    
    methods
        function sf = NonOverlappingSpikeBinFilter(varargin)
            % args: 'binWidthMs', #, binAlignmentMode, BinAlignmentMode.Acausal
            p = inputParser();
            p.addOptional('binWidthMs', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('timeDelta', [], @(x) isnumeric(x) && isscalar(x) || isempty(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            if isempty(p.Results.timeDelta)
                bin = p.Results.binWidthMs;
            else
                bin = p.Results.timeDelta;
            end
            
            sf = sf@ConvolutionSpikeFilter('binWidthMs', bin, 'timeDelta', bin, p.Unmatched);
        end
        
        function v = getTimeDelta(sf, v)
            v = sf.binWidthMs;
        end
            
        % filter used for convolution, as an impulse response which may 
        % have acausal elements if indZero > 1
        function [filt, indZero] = getFilter(sf)
            filt = 1;
            indZero = 1;
        end
        
        % keep binWidthMs == timeDelta
        function sf = postSetBinWidthMs(sf)
            % this would trigger a infinte loop if AbortSet were false,
            % so just in case
            if sf.timeDelta ~= sf.binWidthMs
                sf.timeDelta = sf.binWidthMs;
            end
        end
        
        function checkOkay(sf) % superclass overrides
            assert(sf.timeDelta == sf.binWidthMs, 'timeDelta must match binWidthMs in order to get non-overlapping bins');
        end
    end
 
end
