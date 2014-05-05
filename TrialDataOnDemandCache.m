classdef TrialDataOnDemandCache < handle & matlab.mixin.Copyable
    
    properties
        valid
    end
    
    methods
        function flush(odc)
            odc.valid
        end
    end
    
end
