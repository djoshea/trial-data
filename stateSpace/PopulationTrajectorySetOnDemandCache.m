classdef PopulationTrajectorySetOnDemandCache < handle & matlab.mixin.Copyable 

    properties
        basisNames
        basisUnits
        
        alignSummaryData
        basisAlignSummaryLookup
        
        dataByTrial
        tMinForDataByTrial
        tMaxForDataByTrial
        tMinByTrial
        tMaxByTrial
        alignValidByTrial
        
        dataMean
        dataErrorHigh
        dataErrorLow
        dataValid
        dataNTrials
        tMinForDataMean
        tMaxForDataMean
        tvecDataMean
        nTimeDataMean
        
        dataMeanRandomized
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
            odc.alignSummaryData = {};
            odc.basisAlignSummaryLookup = [];
            
            odc.dataMean = {};
            odc.dataErrorHigh = {};
            odc.dataErrorLow = {};
            odc.dataValid = {};
            odc.dataNTrials = {};
            odc.tMinForDataMean = [];
            odc.tMaxForDataMean = [];
            odc.nTimeDataMean = [];
        end
    end




end
