classdef PopulationTrajectorySetOnDemandCache < handle

    properties
        basisNames
        timeData
        data
        dataSem
        dataValid
        nTrialsData
    end

end

% extract data from trial data

for iBasis = 1:numel(trialDataSet)
     % figure out which channel to extract
     [rateMean, time, rateSem] = td.getSpikeRateFilteredMeanByGroup(unitStr{iBasis});


 end

end
