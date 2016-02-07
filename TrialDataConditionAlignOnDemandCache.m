classdef TrialDataConditionAlignOnDemandCache < TrialDataOnDemandCache
    
    properties
        eventCounts
        eventData
        alignSummarySet
    end
    
    methods
        function flush(odc)
            flush@TrialDataOnDemandCache(odc);
            odc.flushEventData();
        end
        
        function flushEventData(odc)
            odc.eventCounts = [];
            odc.eventData = [];
            odc.flushAlignSummarySet();
        end
        
        function flushAlignSummarySet(odc)
            odc.alignSummarySet = [];
        end
    end
end
