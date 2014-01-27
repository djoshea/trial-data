classdef AlignInfoOnDemandCache < handle & matlab.mixin.Copyable
    properties(SetAccess=?AlignInfo)
        timeInfo
        computedValid
        
        markData
        markCounts
        
        intervalStartData
        intervalStopData
        intervalCounts
    end
    
    methods
        function flush(odc)
            odc.timeInfo = [];
            odc.computedValid = [];

            odc.flushMarkData();
            odc.flushIntervalData();
        end
        
        function flushMarkData(odc)
            odc.markData = [];
            odc.markCounts = [];
        end
        
        function flushIntervalData(odc)
            odc.intervalStartData = [];
            odc.intervalStopData = [];
            odc.intervalCounts = [];
        end
    end
end
