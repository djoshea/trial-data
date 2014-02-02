classdef PopulationTrajectorySetOnDemandCache < handle & matlab.mixin.Copyable 

    properties
        basisNames
        basisUnits
        
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
        
        % trial-averaged data
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
            odc.flushAlignSummaryData();
        end
        
        function flushTrialAveragedData(odc)
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
