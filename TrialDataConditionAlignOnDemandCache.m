classdef TrialDataConditionAlignOnDemandCache < TrialDataOnDemandCache
    
    properties
        eventCounts
        eventData
        alignSummarySet
        
        conditionInfoRandomized 
        randomizedListsByCondition
        
    end
    
    methods
        function flush(odc)
            flush@TrialDataOnDemandCache(odc);
            odc.flushEventData();
        end
        
        function flushRandomized(odc)
            odc.conditionInfoRandomized = [];
            odc.randomizedListsByCondition = [];
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
