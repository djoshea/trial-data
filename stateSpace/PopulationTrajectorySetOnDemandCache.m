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
        
        dataByTrialCommonTimeGrouped
        
        % time windows for trial averaged data
        tMinValidByAlignBasisCondition
        tMaxValidByAlignBasisCondition
        
        % trial-averaged data
        dataMean
        dataSem
        dataValid
        dataNumTrialsRaw
        trialLists
        tMinForDataMean
        tMaxForDataMean

%         % noise estimates from scaled differences of trials
%         dataDifferenceOfTrialsScaledNoiseEstimate
%         
%         % data mean randomized via resampling and high/low intervals
%         dataIntervalHigh
%         dataIntervalLow
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
            odc.flushTrialAveragedData();
        end
        
        function flushTrialAveragedData(odc) % built by buildDataMean
            odc.dataByTrialCommonTimeGrouped = {}; % this is here as it depends on the grouping
            odc.dataMean = {};
            odc.dataSem = {};
            odc.tMinForDataMean = [];
            odc.tMaxForDataMean = [];
            odc.flushValid();
            odc.flushDataNumTrials();
            odc.flushTimeWindowsByAlignBasisCondition();
            odc.flushRandomizedTrialAveragedData();
%             odc.flushDifferenceOfTrialsNoiseEstimate();
            odc.flushAlignSummaryData();
        end
        
        function flushDataNumTrials(odc) % built by builddataNumTrials
            odc.dataValid = {};
            odc.dataNumTrialsRaw = {};
            odc.trialLists = {};
        end
        
        function flushTimeWindowsByAlignBasisCondition(odc) % built by buildTimeWindowsByAlignBasisCondition
            odc.tMinValidByAlignBasisCondition = [];
            odc.tMaxValidByAlignBasisCondition = [];
        end
        
        function flushAlignSummaryData(odc)
            odc.alignSummaryData = {};
            odc.basisAlignSummaryLookup = [];
            odc.alignSummaryAggregated = {};
        end

%         function flushRandomizedTrialAveragedData(odc)
%             odc.dataIntervalHigh = {};
%             odc.dataIntervalLow = {};
%         end
        
%         function flushDifferenceOfTrialsNoiseEstimate(odc)
%             odc.dataDifferenceOfTrialsScaledNoiseEstimate = [];
%         end
    end

end
