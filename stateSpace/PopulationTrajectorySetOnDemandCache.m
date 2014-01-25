classdef PopulationTrajectorySetOnDemandCache < handle & matlab.mixin.Copyable 

    properties
        basisNames
        basisUnits
        
        alignSummaryData
        basisAlignSummaryLookup
        alignSummaryAggregated
        
        dataByTrial
        tMinForDataByTrial
        tMaxForDataByTrial
        tMinByTrial
        tMaxByTrial
        alignValidByTrial
        
        dataMean
        dataSem
        dataValid
        dataNTrials
        tMinForDataMean
        tMaxForDataMean
        tvecDataMean
        nTimeDataMean
        
        dataMeanRandomized
        dataIntervalHigh
        dataIntervalLow
    end

    methods
        function flush(odc)
            odc.basisNames = {};
            odc.basisUnits = {};
           
            odc.dataByTrial = {};
            odc.tMinForDataByTrial = {};
            odc.tMaxForDataByTrial = {};
            odc.tMinByTrial = {};
            odc.tMaxByTrial = {};
            odc.alignValidByTrial = {};
            
            odc.flushTrialAveragedData();
        end
        
        function flushTrialAveragedData(odc)
            odc.alignSummaryAggregated = {};
            odc.alignSummaryData = {};
            odc.basisAlignSummaryLookup = [];
            
            odc.dataMean = {};
            odc.dataSem = {};
            odc.dataValid = {};
            odc.dataNTrials = {};
            odc.tMinForDataMean = [];
            odc.tMaxForDataMean = [];
            odc.nTimeDataMean = [];
            odc.tvecDataMean = [];

            odc.flushRandomizedTrialAveragedData();
        end

        function flushRandomizedTrialAveragedData(odc)
            odc.dataMeanRandomized = {};
            odc.dataIntervalHigh = {};
            odc.dataIntervalLow = {};
        end
    end

end
