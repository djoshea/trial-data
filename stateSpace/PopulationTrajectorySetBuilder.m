classdef PopulationTrajectorySetBuilder
    % these properties are used to construct PopulationTrajectorySets
    % manually by holding property values temporarily in the factory
    % instance. See PopulationTrajectorySet for documentation on the
    % meanings of these values
    
    % These properties hold temporary data for injecting into a PopulationTrajectorySet
    % when building 
    properties
        %% fSettings
        dataUnits
        timeUnitName
        timeUnitsPerSecond
        timeDelta 
        spikeFilter
        minTrialsForTrialAveraging
        minFractionTrialsForTrialAveraging
        includeOnlyTrialsValidAllAlignments
        dataIntervalQuantileHigh
        dataIntervalQuantileLow

        %% fDescriptors
        alignDescriptorSet
        conditionDescriptor
        translationNormalization
        conditionDescriptorRandomized
        
        %% fBasisInfo
        basisNames
        basisUnits
        basisValid
        basisInvalidCause
        
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
        tMinValidByAlignBasisCondition
        tMaxValidByAlignBasisCondition
        tMinForDataMean
        tMaxForDataMean
        dataMean
        dataSem
        dataNTrials
        trialLists
        dataValid
        alignSummaryData
        basisAlignSummaryLookup
        dataDifferenceOfTrialsScaledNoiseEstimate
        
        %% fTrialAvgRandomized
        dataMeanRandomized
        dataSemRandomized
%         dataIntervalHigh
%         dataIntervalLow
    end
    
    % Lists of fields within PopulationTrajectorySet to simplify the
    % copying and assert ~isempty checks below.
    properties(Constant)
        fSettings = {'dataUnits', 'timeUnitName', 'timeUnitsPerSecond', 'timeDelta', 'spikeFilter', 'minTrialsForTrialAveraging', ...
            'minFractionTrialsForTrialAveraging', ...
            'dataIntervalQuantileLow', 'dataIntervalQuantileHigh'};

        fDescriptors = {'alignDescriptorSet', 'conditionDescriptor', 'translationNormalization', 'conditionDescriptorRandomized'};
        
        fBasisInfo = {'basisNames', 'basisUnits', 'basisValid', 'basisInvalidCause'};
        
        fDataSourceInfo = {'dataSources', 'basisDataSourceIdx', 'basisDataSourceChannelNames'};
        
        fSingleTrial = {'dataByTrial', 'tMinForDataByTrial', 'tMaxForDataByTrial', ...
            'alignValidByTrial', 'tMinByTrial', 'tMaxByTrial'};
        
        fTrialAvg = {'tMinValidByAlignBasisCondition', 'tMaxValidByAlignBasisCondition', ...
                'tMinForDataMean', 'tMaxForDataMean', 'dataMean', 'dataSem', ...
                'dataNTrials', 'dataValid', ...
                'alignSummaryData', 'basisAlignSummaryLookup', 'trialLists', ...
                'dataDifferenceOfTrialsScaledNoiseEstimate'}
        
        fTrialAvgRandomized = {'dataMeanRandomized', 'dataSemRandomized'};
        
        fCanBeEmptyExceptions = {'translationNormalization', 'conditionDescriptorRandomized'};
    end
        
    methods(Static)
        function pset = fromAllUnitsInTrialData(td)
            % each unit in TrialData becomes a basis
            units = td.listSpikeChannels();
            nUnits = numel(units);
            
            pset = PopulationTrajectorySet();
            pset.datasetName = td.datasetName;
            
            tdca = TrialDataConditionAlign(td);
            pset.timeUnitName = tdca.timeUnitName;
            pset.timeUnitsPerSecond = tdca.timeUnitsPerSecond;
            pset.dataSources = {tdca};
            pset.dataUnits = 'sp/s';
            pset.basisDataSourceIdx = onesvec(nUnits);
            pset.basisDataSourceChannelNames = units;
            
            pset = pset.setConditionDescriptor(td.conditionInfo);
            pset = pset.setAlignDescriptorSet(td.alignInfoSet);
            pset = pset.initialize();
        end
        
        function pset = fromAllUnitsInTrialDataSplitBySaveTag(td)
            % each unit from each save tag in TrialData becomes a basis
            td = TrialDataConditionAlign(td);
            units = td.listSpikeChannels();
            nUnits = numel(units);
            
            saveTags = td.getParamUnique('saveTag');
            nSaveTags = numel(saveTags);
            [~, saveTagIdxByTrial] = ismember(td.getParam('saveTag'), saveTags);
            
            % number of spikes in each save tag by each unit
            debug('Counting spikes in each save tag for each unit\n');
            spikeCountBySaveTag = zeros(nUnits, nSaveTags);
            for iU = 1:nUnits
                counts = accumarray(saveTagIdxByTrial, td.getSpikeCounts(units{iU}), [nSaveTags 1]); % accumarray is useful
                spikeCountBySaveTag(iU, :) = counts;
            end
            [keepUnits, keepSaveTags] = find(spikeCountBySaveTag > 0);
            
            pset = PopulationTrajectorySet();
            pset.dataUnits = 'sp/s';
            pset.datasetName = td.datasetName;
            pset.timeUnitName = td.timeUnitName;
            pset.timeUnitsPerSecond = td.timeUnitsPerSecond;
            
            % build separate TrialData for each save tag, then index the
            % data sources back into these by save tag
            debug('Splitting save tags\n');
            pset.dataSources = arrayfun(@(st) td.selectTrials(saveTagIdxByTrial == st), 1:numel(saveTags), 'UniformOutput', false);
            pset.basisDataSourceIdx = keepSaveTags;
            pset.basisDataSourceChannelNames = units(keepUnits);

            debug('Initializing PTS\n');
            pset = pset.initialize();
        end
        
        function pset = fromMultipleTrialData(tdCell, channelNames)
            % if only tdCell is provided, all spiking units in each td will
            % be used
            pset = PopulationTrajectorySet();
            %pset.datasetName = td.datasetName;
            
            nSources = numel(tdCell);
            iBasis = 1;
            for i = 1:nSources
                if ~isa(tdCell{i}, 'TrialDataConditionAlign')
                    tdCell{i} = TrialDataConditionAlign(tdCell{i});
                end
                units = tdCell{i}.listSpikeChannels();
                if nargin < 2
                    % use all spike channels in each file
                    for j = 1:numel(units)
                        basisDataSourceIdx(iBasis) = i; %#ok<AGROW,PROP>
                        basisDataSourceChannelNames{iBasis} = units{j}; %#ok<AGROW,PROP>
                        iBasis = iBasis + 1;
                    end
                else
                    % use specified channel names
                    basisDataSourceIdx(iBasis) = i; %#ok<PROP>
                    basisDataSourceChannelNames{iBasis} = channelNames{i}; %#ok<PROP> 
                    iBasis = iBasis + 1;
                end
            end
            
            % check that all tdCell have same timeUnitsPerSecond
            assert(numel(unique(cellfun(@(td) td.timeUnitsPerSecond, tdCell))) == 1, ...
                'All data sources must have identical timeUnitsPerSecond');
            
            pset.timeUnitName = tdCell{1}.timeUnitName;
            pset.timeUnitsPerSecond = tdCell{1}.timeUnitsPerSecond;
            pset.dataSources = makecol(tdCell);
            pset.basisDataSourceIdx = makecol(basisDataSourceIdx); %#ok<PROP>
            pset.basisDataSourceChannelNames = makecol(basisDataSourceChannelNames); %#ok<PROP>
            
            
            td = tdCell{1};
            pset = pset.setConditionDescriptor(td.conditionInfo);
            pset = pset.setAlignDescriptorSet(td.alignInfoSet);
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

        function bld = copyTrialAveragedOnlyFromPopulationTrajectorySet(pset, fields)
            % copy values from population trajectory set, all fields except
            % single trial data
            bld = PopulationTrajectorySetBuilder();
            
            if nargin < 2
                fields = [...
                    PopulationTrajectorySetBuilder.fSettings, ...
                    PopulationTrajectorySetBuilder.fDescriptors, ...
                    PopulationTrajectorySetBuilder.fBasisInfo, ...
                    PopulationTrajectorySetBuilder.fTrialAvg ...
                    PopulationTrajectorySetBuilder.fTrialAvgRandomized
                    ]   ;
            end
            
            for iF = 1:length(fields)
                fld = fields{iF};
                bld.(fld) = pset.(fld);
            end
        end
        
        function bld = copySettingsDescriptorsBasisInfoFromPopulationTrajectorySet(pset)
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset, ...
                [ PopulationTrajectorySetBuilder.fSettings, ...
                  PopulationTrajectorySetBuilder.fDescriptors, ...
                  PopulationTrajectorySetBuilder.fBasisInfo, ...
                  ]);
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
            fields = setdiff(fields, PopulationTrajectorySetBuilder.fCanBeEmptyExceptions);
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
