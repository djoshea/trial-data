classdef TrialDataConditionAlignOnDemandCache < handle & matlab.mixin.Copyable
    
    properties
        alignSummary
    end
    
    methods
        function flush(odc)
            odc.alignSummary = [];
        end
    end
    
end