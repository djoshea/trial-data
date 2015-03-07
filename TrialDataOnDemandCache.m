classdef TrialDataOnDemandCache < handle & matlab.mixin.Copyable
    
    properties
        valid
        invalidCause
    end
    
    methods
        function flush(odc)
            odc.flushValid();
        end
        
        function flushValid(odc)
            odc.valid = [];
            odc.invalidCause = {};
        end
    end
    
end
