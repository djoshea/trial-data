classdef PopulationTrajectorySetBuilder
    % these properties are used to construct PopulationTrajectorySets
    % manually by holding property values temporarily in the factory
    % instance. See PopulationTrajectorySet for documentation on the
    % meanings of these values
   
    % These properties hold temporary data for injecting into a PopulationTrajectorySet
    % when building 
    properties
        %% fSettings
        timeUnitName
        timeUnitsPerSecond
        timeDelta 
        spikeFilter
        minTrialsForTrialAveraging
        minFractionTrialsForTrialAveraging
        includeOnlyTrialsValidAllAlignments
        nRandomSamples
        randomSeed
        dataIntervalQuantileHigh
        dataIntervalQuantileLow

        %% fDescriptors
        alignDescriptorSet
        conditionDescriptor
        translationNormalization
        
        %% fBasisInfo
        basisNames
        basisUnits
        
        %% fDataSourceInfo
        dataSources
        basisDataSourceIdx
        basisDataSourceChannelNames
        
        %% fSingleTrial
        dataByTrial
        tMinForDataByTrial
        tMaxForDataByTrial
        alignValidByTrial
        tMinByTrial
        tMaxByTrial
        
        %% fTrialAvg
        tMinForDataMean
        tMaxForDataMean
        dataMean
        dataSem
        dataNTrials
        dataValid
        alignSummaryData
        basisAlignSummaryLookup
        
        %% fTrialAvgRandomized
        dataMeanRandomized
        dataIntervalHigh
        dataIntervalLow
    end
    
    % Lists of fields within PopulationTrajectorySet to simplify the
    % copying and assert ~isempty checks below.
    properties(Constant)
        fSettings = {'timeUnitName', 'timeUnitsPerSecond', 'timeDelta', 'spikeFilter', 'minTrialsForTrialAveraging', ...
            'minFractionTrialsForTrialAveraging', 'includeOnlyTrialsValidAllAlignments', ...
            'nRandomSamples', 'randomSeed', 'dataIntervalQuantileLow', 'dataIntervalQuantileHigh'};

        fDescriptors = {'alignDescriptorSet', 'conditionDescriptor', 'translationNormalization'};
        
        fBasisInfo = {'basisNames', 'basisUnits'};
        
        fDataSourceInfo = {'dataSources', 'basisDataSourceIdx', 'basisDataSourceChannelNames'};
        
        fSingleTrial = {'dataByTrial', 'tMinForDataByTrial', 'tMaxForDataByTrial', ...
            'alignValidByTrial', 'tMinByTrial', 'tMaxByTrial'};
        
        fTrialAvg = {'tMinForDataMean', 'tMaxForDataMean', 'dataMean', 'dataSem', ...
                'dataNTrials', 'dataValid', ...
                'alignSummaryData', 'basisAlignSummaryLookup'}
        
        fTrialAvgRandomized = {'dataIntervalHigh', 'dataIntervalLow', 'dataMeanRandomized'};
    end
        
    methods(Static)
        function pset = fromAllUnitsInTrialData(td)
            units = td.listSpikeUnits();
            nUnits = numel(units);
            
            pset = PopulationTrajectorySet();
            pset.datasetName = td.datasetName;
            
            tdca = TrialDataConditionAlign(td);
            pset.timeUnitName = tdca.timeUnitName;
            pset.timeUnitsPerSecond = tdca.timeUnitsPerSecond;
            pset.dataSources = {tdca};
            pset.basisDataSourceIdx = onesvec(nUnits);
            pset.basisDataSourceChannelNames = units;
            
            pset = pset.initialize();
        end
        
        function pset = fromMultipleTrialData(tdCell, channelNames)
            pset = PopulationTrajectorySet();
            %pset.datasetName = td.datasetName;
            
            nSources = numel(tdCell);
            for i = 1:nSources
                if ~isa(tdCell{i}, 'TrialDataConditionAlign')
                    tdCell{i} = TrialDataConditionAlign(tdCell{i});
                end
            end
            
            % check that all tdCell have same timeUnitsPerSecond
            assert(numel(unique(cellfun(@(td) td.timeUnitsPerSecond, tdCell))) == 1, ...
                'All data sources must have identical timeUnitsPerSecond');
            
            pset.timeUnitName = tdCell{1}.timeUnitName;
            pset.timeUnitsPerSecond = tdCell{1}.timeUnitsPerSecond;
            pset.dataSources = makecol(tdCell);
            pset.basisDataSourceIdx = makecol(1:nSources);
            pset.basisDataSourceChannelNames = makecol(channelNames);
            
            pset = pset.initialize();
        end
        
        function pset = fromAnalogChannelsInTrialData(td)
            chNames = td.listAnalogChannels();
            
            pset = PopulationTrajectorySet();
            pset.datasetName = td.datasetName;
            
            tdca = TrialDataConditionAlign(td);
            pset.timeUnitName = tdca.timeUnitName;
            pset.timeUnitsPerSecond = tdca.timeUnitsPerSecond;
            pset.dataSources = {tdca};
            pset.basisDataSourceIdx = onesvec(numel(chNames));
            pset.basisDataSourceChannelNames = chNames;
            
            pset = pset.initialize();
        end
    end
    
    methods
        function pset = buildManualWithSingleTrialData(bld)
            bld.assertNonEmpty(...
                PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                PopulationTrajectorySetBuilder.fDataSourceInfo, ...
                PopulationTrajectorySetBuilder.fSingleTrial, ...
                PopulationTrajectorySetBuilder.fTrialAvg);
            
            pset = bld.buildManualWithFields(...
                PopulationTrajectorySetBuilder.fSettings, ...
                PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                PopulationTrajectorySetBuilder.fDataSourceInfo, ...
                PopulationTrajectorySetBuilder.fSingleTrial, ...
                PopulationTrajectorySetBuilder.fTrialAvg, ...
                PopulationTrajectorySetBuilder.fTrialAvgRandomized);
        end
        
        function pset = buildManualWithTrialAveragedData(bld)
%             bld.assertNonEmpty(...
%                 PopulationTrajectorySetBuilder.fDescriptors, ...
%                 PopulationTrajectorySetBuilder.fBasisInfo, ...
%                 PopulationTrajectorySetBuilder.fTrialAvg);
            
            pset = bld.buildManualWithFields(...
                PopulationTrajectorySetBuilder.fSettings, ...
                PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                PopulationTrajectorySetBuilder.fTrialAvg, ...
                PopulationTrajectorySetBuilder.fTrialAvgRandomized);
        end
    end
    
    methods(Static)
        function bld = copyFromPopulationTrajectorySet(pset, fields)
            % copy values from population trajectory set, all fields by default
            bld = PopulationTrajectorySetBuilder();
            
            if nargin < 2
                fields = [...
                    PopulationTrajectorySetBuilder.fSettings, ...
                    PopulationTrajectorySetBuilder.fDescriptors, ...
                    PopulationTrajectorySetBuilder.fBasisInfo, ...
                    PopulationTrajectorySetBuilder.fDataSourceInfo, ...
                    PopulationTrajectorySetBuilder.fSingleTrial, ...
                    PopulationTrajectorySetBuilder.fTrialAvg ...
                    PopulationTrajectorySetBuilder.fTrialAvgRandomized
                    ]   ;
            end
            
            for iF = 1:length(fields)
                fld = fields{iF};
                bld.(fld) = pset.(fld);
            end
        end

        function bld = copySettingsDescriptorsFromPopulationTrajectorySet(pset)
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset, ...
                [ PopulationTrajectorySetBuilder.fSettings, ...
                  PopulationTrajectorySetBuilder.fDescriptors ]);
        end
        
        function psetManual = convertToManualWithSingleTrialData(pset)
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset);
            psetManual = bld.buildManualWithSingleTrialData();
        end
        
        function psetManual = convertToManualWithTrialAveragedData(pset)
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset);
            psetManual = bld.buildManualWithTrialAveragedData();
        end
    end
    
    methods(Access=protected)
        function assertNonEmpty(bld, varargin)
            % assert that all fields in field are non empty
            fields = [varargin{:}];
            isEmpty = cellfun(@(fld) isempty(bld.(fld)), fields);
            if any(isEmpty)
                error('Values must be specified for these fields: %s', ...
                    strjoin(fields(isEmpty), ', '));
            end
        end
        
        function pset = buildManualWithFields(bld, varargin)
            pset = PopulationTrajectorySet();
            pset.dataSourceManual = true;
            
            for iArg = 1:numel(varargin)
                fields = varargin{iArg};
                for iFld = 1:numel(fields)
                    fld = fields{iFld};
                    pset.(fld) = bld.(fld);
                end 
            end
            
            pset = pset.initialize();
        end
    end
end
