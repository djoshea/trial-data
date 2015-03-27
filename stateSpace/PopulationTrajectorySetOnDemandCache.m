classdef PopulationTrajectorySetOnDemandCache < handle & matlab.mixin.Copyable 

    properties
        basisNames
        basisUnits 
        
        basisValid
        basisInvalidCause
        
        alignSummaryData
        basisAlignSummaryLookup
        alignSummaryAggregated
        
        % single trial aligned data
        dataByTrial
        tMinForDataByTrial
        tMaxForDataByTrial
        tMinByTrial
        tMaxByTrial
        alignValidByTrial
        
        % time windows for trial averaged data
        tMinValidByAlignBasisCondition
        tMaxValidByAlignBasisCondition
        
        % trial-averaged data
        dataMean
        dataSem
        dataValid
        dataNTrials
        trialLists
        tMinForDataMean
        tMaxForDataMean
        
        dataMeanRandomized
        dataIntervalHigh
        dataIntervalLow
    end

    methods
        function flush(odc)
            odc.flushBasisInfo();
            odc.flushDataByTrial();
            odc.flushTrialAveragedData();
            odc.flushAlignSummaryData();
        end
        
        function flushBasisInfo(odc)
            odc.basisNames = {};
            odc.basisUnits = {};
            odc.basisValid = [];
        end
        
        function flushValid(odc)
            odc.basisValid = [];
            odc.basisInvalidCause = {};
        end
        
        function flushDataByTrial(odc)
            odc.dataByTrial = {};
            odc.tMinForDataByTrial = {};
            odc.tMaxForDataByTrial = {};
            odc.tMinByTrial = {};
            odc.tMaxByTrial = {};
            odc.alignValidByTrial = {};
        end
        
        function flushTrialAveragedData(odc)
            odc.dataMean = {};
            odc.dataSem = {};
            odc.dataValid = {};
            odc.dataNTrials = {};
            odc.trialLists = {};
            odc.tMinForDataMean = [];
            odc.tMaxForDataMean = [];
            odc.tMinValidByAlignBasisCondition = [];
            odc.tMaxValidByAlignBasisCondition = [];
            odc.flushRandomizedTrialAveragedData();
        end
        
        function flushAlignSummaryData(odc)
            odc.alignSummaryData = {};
            odc.basisAlignSummaryLookup = [];
            odc.alignSummaryAggregated = {};
        end

        function flushRandomizedTrialAveragedData(odc)
            odc.dataMeanRandomized = {};
            odc.dataIntervalHigh = {};
            odc.dataIntervalLow = {};
        end
        
    end

end
