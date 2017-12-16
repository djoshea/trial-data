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
        spikeFilter
        minTrialsForTrialAveraging
        minFractionTrialsForTrialAveraging
        ignoreLeadingTrailingZeroSpikeTrials
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
        
        basisValidManual
        basisInvalidCauseManual
        
        basisValidTemporary
        basisInvalidCauseTemporary
               
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
        
        % fTrialAvgCrossValidation
        dataCachedSampledTrialsTensor
        dataCachedMeanExcludingSampledTrialsTensor
        dataCachedSampledTrialCounts
        
        %% fDiffTrialsNoise
        dataDifferenceOfTrialsScaledNoiseEstimate
        
        %% fTrialAvgRandomized
        dataMeanRandomized
        dataSemRandomized
        dataDifferenceOfTrialsScaledNoiseEstimateRandomized
%         dataIntervalHigh
%         dataIntervalLow
    end
    
    % Lists of fields within PopulationTrajectorySet to simplify the
    % copying and assert ~isempty checks below.
    properties(Constant)
        fSettings = {'dataUnits', 'timeUnitName', 'timeUnitsPerSecond', 'spikeFilter', 'minTrialsForTrialAveraging', ...
            'minFractionTrialsForTrialAveraging', 'ignoreLeadingTrailingZeroSpikeTrials', ...
            'dataIntervalQuantileLow', 'dataIntervalQuantileHigh'};

        fDescriptors = {'alignDescriptorSet', 'conditionDescriptor', 'translationNormalization', 'conditionDescriptorRandomized'};
        
        fBasisInfo = {'basisNames', 'basisUnits', ...
            'basisValidManual', 'basisInvalidCauseManual', ...
            'basisValidTemporary', 'basisInvalidCauseTemporary'};
        
        fDataSourceInfo = {'dataSources', 'basisDataSourceIdx', 'basisDataSourceChannelNames'};
        
        fSingleTrial = {'dataSources', 'dataByTrial', 'tMinForDataByTrial', 'tMaxForDataByTrial', ...
            'alignValidByTrial', 'tMinByTrial', 'tMaxByTrial'};
        
        fTrialAvg = {'tMinValidByAlignBasisCondition', 'tMaxValidByAlignBasisCondition', ...
                'tMinForDataMean', 'tMaxForDataMean', 'dataMean', 'dataSem', ...
                'dataNTrials', 'dataValid', ...
                'alignSummaryData', 'basisAlignSummaryLookup', 'trialLists'};
            
        fTrialAvgCrossValidation = {'dataCachedSampledTrialsTensor', 'dataCachedMeanExcludingSampledTrialsTensor', ...
                'dataCachedSampledTrialCounts'};
            
        fDiffTrialsNoise = {'dataDifferenceOfTrialsScaledNoiseEstimate'};
        
        fTrialAvgRandomized = {'dataMeanRandomized', 'dataSemRandomized', ...
            'dataDifferenceOfTrialsScaledNoiseEstimateRandomized'};
        
        fCanBeEmptyExceptions = {'translationNormalization', 'conditionDescriptorRandomized', ...
            'basisInvalidCauseTemporary', 'basisValidTemporary', ...
            'basisValidManual', 'basisInvalidCauseManual'};
    end
        
    methods(Static)
        function pset = fromAllUnitsInTrialData(tdca)
            if ~isa(tdca, 'TrialDataConditionAlign')
                tdca = TrialDataConditionAlign(tdca);
            end
            
            % each unit in TrialData becomes a basis
            units = tdca.listSpikeChannels();
            nUnits = numel(units);
            
            pset = PopulationTrajectorySet();
            pset.datasetName = tdca.datasetName;
            
            pset.timeUnitName = tdca.timeUnitName;
            pset.timeUnitsPerSecond = tdca.timeUnitsPerSecond;
            pset.dataSources = {tdca};
            pset.dataUnits = 'sp/s';
            pset.basisDataSourceIdx = onesvec(nUnits);
            pset.basisDataSourceChannelNames = units;
            
            ci = tdca.conditionInfo.fixAllValueLists();
            pset = pset.setConditionDescriptor(ci);
            pset = pset.setAlignDescriptorSet(tdca.alignInfoSet);
            pset = pset.initialize();
        end
        
        function pset = fromAllUnitsInTrialDataSplitBySaveTag(td)
            % each unit from each save tag in TrialData becomes a basis
            if ~isa(td, 'TrialDataConditionAlign')
                td = TrialDataConditionAlign(td);
            end
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
            pset.dataSources = arrayfun(@(st) td.selectTrials(saveTagIdxByTrial == st), (1:numel(saveTags))', 'UniformOutput', false);
            pset.basisDataSourceIdx = keepSaveTags;
            pset.basisDataSourceChannelNames = units(keepUnits);

            debug('Initializing PTS\n');
            pset = pset.setConditionDescriptor(td.conditionInfo);
            pset = pset.setAlignDescriptorSet(td.alignInfoSet);
            pset = pset.initialize();
        end
        
        function pset = fromMultipleTrialData(tdCell, varargin)
            p = inputParser();
            p.addOptional('channelNames', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.parse(varargin{:});
            
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
                if isempty(p.Results.channelNames)
                    % use all spike channels in each file
                    for j = 1:numel(units)
                        basisDataSourceIdx(iBasis) = i; %#ok<AGROW>
                        basisDataSourceChannelNames{iBasis} = units{j}; %#ok<AGROW>
                        iBasis = iBasis + 1;
                    end
                else
                    % use specified channel names
                    basisDataSourceIdx(iBasis) = i;
                    basisDataSourceChannelNames{iBasis} = p.Results.channelNames{i}; 
                    iBasis = iBasis + 1;
                end
            end
            
            % check that all tdCell have same timeUnitsPerSecond
            assert(numel(unique(cellfun(@(td) td.timeUnitsPerSecond, tdCell))) == 1, ...
                'All data sources must have identical timeUnitsPerSecond');
            
            pset.timeUnitName = tdCell{1}.timeUnitName;
            pset.timeUnitsPerSecond = tdCell{1}.timeUnitsPerSecond;
            pset.dataSources = makecol(tdCell);
            pset.basisDataSourceIdx = makecol(basisDataSourceIdx); 
            pset.basisDataSourceChannelNames = makecol(basisDataSourceChannelNames); 
            if ~isempty(p.Results.spikeFilter)
                pset.spikeFilter = p.Results.spikeFilter;
            end
            
            tdca = tdCell{1};
            pset = pset.setConditionDescriptor(tdca.conditionInfo);
            pset = pset.setAlignDescriptorSet(tdca.alignInfoSet);
            pset = pset.initialize();
        end
        
        function pset = fromAnalogChannelsInTrialData(td, chNames, varargin)
            if nargin < 2
                chNames = td.listAnalogChannels();
            end
            p = inputParser();
            p.addParameter('timeDelta', td.getAnalogTimeDelta(chNames), @isscalar);
            p.parse(varargin{:});
            
            pset = PopulationTrajectorySet();
            pset.datasetName = td.datasetName;
            
            if ~isa(td, 'TrialDataConditionAlign')
                td = TrialDataConditionAlign(td);
            end
            pset.timeUnitName = td.timeUnitName;
            pset.timeUnitsPerSecond = td.timeUnitsPerSecond;
            pset.dataSources = {td};
            pset.basisDataSourceIdx = onesvec(numel(chNames));
            pset.basisDataSourceChannelNames = chNames;
            
            % don't want spiking filter to add padding
            pset.spikeFilter = NonOverlappingSpikeBinFilter('timeDelta', p.Results.timeDelta);
            
            pset = pset.setConditionDescriptor(td.conditionInfo);
            pset = pset.setAlignDescriptorSet(td.alignInfoSet);
            pset = pset.initialize();
        end
    end
    
    methods(Static)
        function [toCopy, toCheckNonEmpty] = listFieldsSingleTrial(varargin)
            p = inputParser();
            p.addParameter('includeDataSourceInfo', false, @islogical); % used when construction psets with dataSourceManual == true
            p.parse(varargin{:});
            
            toCopy = [PopulationTrajectorySetBuilder.fSettings, ...
                PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                PopulationTrajectorySetBuilder.fDataSourceInfo, ...
                PopulationTrajectorySetBuilder.fSingleTrial, ...
                PopulationTrajectorySetBuilder.fTrialAvg, ...
                PopulationTrajectorySetBuilder.fTrialAvgRandomized, ...
                PopulationTrajectorySetBuilder.fTrialAvgCrossValidation];
            
            toCheckNonEmpty = [PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                ... %PopulationTrajectorySetBuilder.fDataSourceInfo, ...
                PopulationTrajectorySetBuilder.fSingleTrial, ...
                PopulationTrajectorySetBuilder.fTrialAvg];
            
            if p.Results.includeDataSourceInfo
                toCopy = [toCopy, PopulationTrajectorySetBuilder.fDataSourceInfo];
                toCheckNonEmpty = [toCheckNonEmpty, PopulationTrajectorySetBuilder.fDataSourceInfo];
            end
        end
        
        function [toCopy, toCheckNonEmpty] = listFieldsTrialAverage(varargin)
            p = inputParser();
            p.addParameter('includeDiffTrialsNoise', true, @islogical); % this can be slow so we make it optional
            p.addParameter('includeDataSourceInfo', false, @islogical); % used when construction psets with dataSourceManual == true
            p.parse(varargin{:});
            
            toCopy = [PopulationTrajectorySetBuilder.fSettings, ...
                PopulationTrajectorySetBuilder.fDescriptors, ...
                PopulationTrajectorySetBuilder.fBasisInfo, ...
                PopulationTrajectorySetBuilder.fTrialAvg, ...
                PopulationTrajectorySetBuilder.fTrialAvgRandomized, ...
                PopulationTrajectorySetBuilder.fTrialAvgCrossValidation];
            
            toCheckNonEmpty = {};
            
            if p.Results.includeDiffTrialsNoise
                toCopy = [toCopy, PopulationTrajectorySetBuilder.fDiffTrialsNoise];
            end
            
            if p.Results.includeDataSourceInfo
                toCopy = [toCopy, PopulationTrajectorySetBuilder.fDataSourceInfo];
                toCheckNonEmpty = [toCheckNonEmpty, PopulationTrajectorySetBuilder.fDataSourceInfo];
            end
        end
    end
    
    methods(Static)
        function bld = copyFromPopulationTrajectorySet(pset, fields)
            % copy values from population trajectory set, all fields by default
            bld = PopulationTrajectorySetBuilder();
            
            if nargin < 2
                fields = PopulationTrajectorySetBuilder.listFieldsSingleTrial();
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
                fields = PopulationTrajectorySetBuilder.listFieldsTrialAverage();
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
        
        function psetManual = convertToManualWithSingleTrialData(pset, varargin)
            fields = PopulationTrajectorySetBuilder.listFieldsSingleTrial(varargin{:});
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset, fields);
            psetManual = bld.buildManualWithSingleTrialData(varargin{:});
        end
        
        function psetManual = convertToManualWithTrialAveragedData(pset, varargin)
            % params include 'includeDiffTrialsNoise'
            fields = PopulationTrajectorySetBuilder.listFieldsTrialAverage(varargin{:});
            bld = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset, fields);
            psetManual = bld.buildManualWithTrialAveragedData(varargin{:});
        end
        
        function pset = concatenatePopulationTrajectorySets(psetCell)
            if numel(psetCell) == 1
                pset = psetCell{1};
                return;
            end
            
            bld = PopulationTrajectorySetBuilder.copySettingsDescriptorsFromPopulationTrajectorySet(psetCell{1});
            
            % determine simultaneity (when the data source is the same
            % across each
            simultaneous = all(cellfun(@(pset) pset.simultaneous, psetCell));
            if simultaneous
                commonDataSource = psetCell{1}.dataSources{1};
                for iP = 2:numel(psetCell)
                    simultaneous = simultaneous && isequal(commonDataSource, psetCell{iP}.dataSources{1});
                end
            end
            
            manual = cellfun(@(pset) pset.dataSourceManual, psetCell);
            if any(manual)
                error('Not yet implemented for manual psets');
            end
 
            % then copy the data
            fields = PopulationTrajectorySetBuilder.fBasisInfo;   
            getFn = @(fld) cellfun(@(pset) makecol(pset.(fld)), psetCell, 'UniformOutput', false);
            for iF = 1:length(fields)
                fld = fields{iF};
                vals = getFn(fld);
                bld.(fld) = cat(1, vals{:});
            end
            
            nBasesTotal = sum(cellfun(@(pset) pset.nBases, psetCell));
            
            if simultaneous
                bld.dataSources = pset.dataSource;
                bld.basisDataSourceIdx = ones(nBasesTotal, 1);
            else
                vals = getFn('dataSources');
                bld.dataSources = cat(1, vals{:});
                vals = getFn('basisDataSourceIdx');
            
                % offset the data source idx according to the concatenation
                nSources = cellfun(@(pset) pset.nDataSources, psetCell);
                [vals, whichPset] = TensorUtils.catWhich(1, vals{:});
                accum = nSources(1);
                for iP = 2:numel(psetCell)
                    vals(whichPset == iP) = vals(whichPset == iP) + accum;
                    accum = accum + nSources(iP);
                end
                bld.basisDataSourceIdx = vals;
            end
            
            vals = getFn('basisDataSourceChannelNames');
            bld.basisDataSourceChannelNames = cat(1, vals{:});
            
            if simultaneous
                pset = bld.buildAutoWithSingleTrialData();
            else
                pset = bld.buildAutoWithTrialAveragedData();
            end
        end

    end
    
    methods
        function pset = buildAutoWithSingleTrialData(bld)
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.listFieldsSingleTrial('includeDataSourceInfo', true);
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildAutoWithFields(toCopy);
        end
        
        function pset = buildAutoWithTrialAveragedData(bld, varargin)
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.listFieldsTrialAverage(...
                'includeDiffTrialsNoise', true, 'includeDataSourceInfo', true);
            
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildAutoWithFields(toCopy);
        end
            
        function pset = buildManualWithSingleTrialData(bld)
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.listFieldsSingleTrial();
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildManualWithFields(toCopy);
        end
        
        function pset = buildManualWithTrialAveragedData(bld, varargin)
            p = inputParser();
            p.addParameter('includeDiffTrialsNoise', true, @islogical); % this can be slow so we make it optional
            p.parse(varargin{:});
            
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.listFieldsTrialAverage(...
                'includeDiffTrialsNoise', p.Results.includeDiffTrialsNoise);
            
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildManualWithFields(toCopy);
        end
    end
    
    methods(Access=protected)
        function assertNonEmpty(bld, varargin)
            % assert that all fields in field are non empty
            fields = cat(1, varargin{:});
            fields = setdiff(fields, PopulationTrajectorySetBuilder.fCanBeEmptyExceptions);
            isEmpty = cellfun(@(fld) isempty(bld.(fld)), fields);
            if any(isEmpty)
                error('Values must be specified for these fields: %s', ...
                    strjoin(fields(isEmpty), ', '));
            end
        end
        
        function pset = buildManualWithFields(bld, fields)
            pset = PopulationTrajectorySet();
            pset.dataSourceManual = true;
            
            for iFld = 1:numel(fields)
                fld = fields{iFld};
                pset.(fld) = bld.(fld);
            end 
            
            pset = pset.initialize();
            pset.dataSourceManual = true;
        end
        
        function pset = buildAutoWithFields(bld, fields)
            pset = PopulationTrajectorySet();
            pset.dataSourceManual = false;
            
            for iFld = 1:numel(fields)
                fld = fields{iFld};
                pset.(fld) = bld.(fld);
            end 
            
            pset = pset.initialize();
            pset.dataSourceManual = false;
            pset = pset.applyAlignDescriptorSet();
            pset = pset.applyConditionDescriptor();
        end
    end
end
