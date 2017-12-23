classdef NonOverlappingCausalSpikeBinFilter < NonOverlappingSpikeBinFilter 
% when using this class, be sure to set .timeDelta == .binWidthMs so that
% the time windows are sampled evenly. Otherwise an error will be thrown.

    methods
        function sf = NonOverlappingCausalSpikeBinFilter(varargin)
            sf = sf@NonOverlappingSpikeBinFilter(varargin{:});
            sf.binAlignmentMode = BinAlignmentMode.Causal;
        end
        
        function v = getBinAlignmentMode(sf, v)
           v = BinAlignmentMode.Causal;
        end
        
        function checkOkay(sf) % superclass overrides
            checkOkay@NonOverlappingSpikeBinFilter(sf);
            assert(sf.binAlignmentMode == BinAlignmentMode.Causal);
        end
    end
 
end
