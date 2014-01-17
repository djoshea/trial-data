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
        dataSem
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
            odc.dataSem = {};
            odc.dataValid = {};
            odc.dataNTrials = {};
            odc.tMinForDataMean = [];
            odc.tMaxForDataMean = [];
            odc.nTimeDataMean = [];
        end

%         function selectAlongDimension(odc, dim, idx)
%             fields = properties(odc);
%             for i = 1:numel(fields)
%                 odc.(fields{i}) = TensorUtils.selectAlongDimension(odc.(fields{i}), dim, idx);
%             end
%         end
    end




end
