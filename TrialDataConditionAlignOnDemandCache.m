classdef TrialDataConditionAlignOnDemandCache < handle & matlab.mixin.Copyable
    
    properties
        eventCounts
        eventData
        alignSummarySet
    end
    
    methods
        function flush(odc)
            odc.alignSummarySet = [];
            odc.eventCounts = [];
            odc.eventData = [];
        end
    end
    
end
