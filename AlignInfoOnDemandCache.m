classdef AlignInfoOnDemandCache < handle & matlab.mixin.Copyable
    properties(SetAccess=?AlignInfo)
        timeInfo
        computedValid
        
        markData
        
        intervalStartData
        intervalStopData
    end
    
    methods
        function flush(odc)
            odc.timeInfo = [];
            odc.computedValid = [];
            odc.markData = [];
            odc.intervalStartData = [];
            odc.intervalStopData = [];
        end
        
        function flushMarkData(odc)
            odc.markData = [];
        end
        
        function flushIntervalData(odc)
            odc.intervalStartData = [];
            odc.intervalStopData = [];
        end
    end
end
