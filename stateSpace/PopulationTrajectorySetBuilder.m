classdef PopulationTrajectorySetBuilder
    methods(Static)
        function pset = fromAllUnitsInTrialData(td)
            units = td.listSpikeUnits();
            nUnits = numel(units);
            
            pset = PopulationTrajectorySet();
            
            tdca = TrialDataConditionAlign(td);
            pset.dataSources = {tdca};
            pset.basisDataSourceIdx = onesvec(nUnits);
            pset.basisDataSourceChannelNames = units;
            
            pset = pset.initialize();
        end
    end
end