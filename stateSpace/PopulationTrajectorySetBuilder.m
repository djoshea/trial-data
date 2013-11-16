classdef PopulationTrajectorySetBuilder
    properties(SetAccess=protected)
        % if true, the raw sources used to build the data will be maintained within
        % this instance, which allows things like resampling and shuffling
        % to operate on the raw data sources
        hasTrialData

        % a cell array of unknown length which stores the union of all trial data 
        % instances used by all bases
        trialDataSet

        % a pointer list of size nBases x 1 which indicates which element of 
        % trialDataSources a particular basis derives its data 
        trialDataIndByBasis
    end
    
    properties(Dependent)
        trialDataByBasis % a dynamically computed cell array that returns the trial data associated with each basis
    end
end
