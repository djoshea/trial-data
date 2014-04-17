classdef AlignInfoOnDemandCache < handle & matlab.mixin.Copyable
    properties(SetAccess=?AlignInfo)
        timeInfo
        timeInfoValid
        
        computedValid
        
        markData
        markCounts
        
        markDataValid
        markCountsValid
        
        intervalStartData
        intervalStopData
        intervalCounts
        
        intervalStartDataValid
        intervalStopDataValid
        intervalCountsValid
    end
    
    methods
        function flush(odc)
            odc.timeInfo = [];
            odc.computedValid = [];

            odc.flushMarkData();
            odc.flushIntervalData();
            odc.flushManualValid();
        end
        
        function flushMarkData(odc)
            odc.markData = [];
            odc.markCounts = [];
            
            odc.markDataValid = [];
            odc.markCountsValid = [];
        end
        
        function flushIntervalData(odc)
            odc.intervalStartData = [];
            odc.intervalStopData = [];
            odc.intervalCounts = [];
            
            odc.intervalStartDataValid = [];
            odc.intervalStopDataValid = [];
            odc.intervalCountsValid = [];
        end
        
        function flushManualValid(odc)
            odc.timeInfoValid = [];
            
            odc.markDataValid = [];
            odc.markCountsValid = [];
            
            odc.intervalStartDataValid = [];
            odc.intervalStopDataValid = [];
            odc.intervalCountsValid = [];
        end
    end
end
