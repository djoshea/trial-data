classdef PopulationTrajectorySetOnDemandCache < handle & matlab.mixin.Copyable 

    properties
        basisNames
        basisUnits
        
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
            
            odc.flushTrialAveragedData
        end
        
        function flushTrialAveragedData(odc)
            odc.dataMean = {};
            odc.dataSem = {};
            odc.dataValid = {};
            odc.dataNTrials = {};
            odc.tMinForDataMean = {};
            odc.tMaxForDataMean = {};
        end

%         function selectAlongDimension(odc, dim, idx)
%             fields = properties(odc);
%             for i = 1:numel(fields)
%                 odc.(fields{i}) = TensorUtils.selectAlongDimension(odc.(fields{i}), dim, idx);
%             end
%         end
    end




end
