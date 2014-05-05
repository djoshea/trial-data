classdef TrialDataConditionAlignOnDemandCache < TrialDataOnDemandCache
    
    properties
        eventCounts
        eventData
        alignSummarySet
    end
    
    methods
        function flush(odc)
            flush@TrialDataOnDemandCache(odc);
            odc.alignSummarySet = [];
            odc.eventCounts = [];
            odc.eventData = [];
        end
    end
    
end
