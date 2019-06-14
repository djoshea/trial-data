classdef AnalogDataArbitrarySampling
   
    properties
        channelNames (:, 1) string
        
        data cell % nTrials x nAlign cell of nTime x 1 samples
        time cell % nTrials x nAlign cell of nTime x 1 vectors
        
        alignInfoSet
        conditionInfo
    end
    
    methods
        function d = AnalogDataArbitrarySampling(varargin)
            p = inputParser();
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            s = p.Unmatched;
            flds = fieldnames(s);
           
            for iF = 1:numel(flds)
                d.(flds{iF}) = s.(flds{iF});
            end
        end
    
    end
    
end