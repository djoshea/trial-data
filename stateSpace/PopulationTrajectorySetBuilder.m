classdef PopulationTrajectorySetBuilder
    % these properties are used to construct PopulationTrajectorySets
    % manually by holding property values temporarily in the factory
    % instance. See PopulationTrajectorySet for documentation on the
    % meanings of these values
    
    % These properties hold temporary data for injecting into a PopulationTrajectorySet
    % when building 
    properties
        data
        propMeta
    end
   
    % initial construction methods
    methods(Static)
        function pset = fromAllUnitsInTrialData(tdca, varargin)
            p = inputParser();
            p.addParameter('spikeFilter', [], @(x) isempty(x) || isa(x, 'SpikeFilter'));
            p.parse(varargin{:});
            
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
            if ~isempty(p.Results.spikeFilter)
                pset.spikeFilter = p.Results.spikeFilter;
            end
            
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
        
        function pset = fromDataMeanTensor(dataMean_NbyCbyT, varargin)
            p = inputParser();
            
            p.addOptional('time', [], @(x) isempty(x) || iscell(x) || isvector(x));
            p.addParameter('basisNames', {}, @iscellstr);
            p.addParameter('zeroEventNames', {}, @iscellstr);
            p.addParameter('conditionsSize', [], @(x) isempty(x) || isvector(x));
            p.addParameter('axisNames', {}, @iscellstr);
            p.addParameter('valuesAlongAxes', {}, @iscell);
            p.addParameter('dataSem', [], @(x) isempty(x) || iscell(x) || isnumeric(x));
            p.addParameter('timeUnitName', 'ms', @ischar);
            p.addParameter('timeUnitsPerSecond', 1000, @isscalar);
            p.parse(varargin{:});
            
            if ~iscell(dataMean_NbyCbyT)
                dataMean_NbyCbyT = {dataMean_NbyCbyT};
            end
            if isempty(p.Results.dataSem)
                dataSem = cellfun(@(x) nan(size(x)), dataMean_NbyCbyT, 'UniformOutput', false);
            else
                dataSem = p.Results.dataSem;
                if ~iscell(dataSem)
                    dataSem = {dataSem};
                end
            end
            T = cellfun(@(d) size(d, 3), dataMean_NbyCbyT);
            C = size(dataMean_NbyCbyT{1}, 2);
            N = size(dataMean_NbyCbyT{1}, 1);
            if isempty(p.Results.conditionsSize)
                conditionsSize = C;
            else
                conditionsSize = p.Results.conditionsSize;
                assert(prod(conditionsSize) == C);
            end
            cd = ConditionDescriptor.createManualWithSize(conditionsSize, 'axisNames', p.Results.axisNames, 'valuesAlongAxes', p.Results.valuesAlongAxes);
                
            if isempty(p.Results.time)
                time = (1:T)';
            else
                time = p.Results.time;
            end
            if ~iscell(time)
                time = {time};
            end
            timeDelta = time{1}(2) - time{1}(1);
            tMinDataMean = cellfun(@min, time);
            tMaxDataMean = cellfun(@max, time);
            A = numel(dataMean_NbyCbyT);
            
            if isempty(p.Results.zeroEventNames)
                zeroEventNames = arrayfun(@(i) sprintf('Align%d', i), (1:A)', 'UniformOutput', false);
            end
            alignDescriptorSet = cell(A, 1);
            for iA = 1:A
                ev = zeroEventNames{iA};
                alignDescriptorSet{iA} = AlignDescriptor().zero(ev).start(ev, tMinDataMean(iA)).stop(ev, tMaxDataMean(iA));
            end
            
            [tMinValidByAlignBasisCondition, tMaxValidByAlignBasisCondition] = deal(nan(A, N, C));
            dataValid = false(A, N, C);
            for iA = 1:A
                tMinValidByAlignBasisCondition(iA, :, :) = TensorUtils.findNAlongDim(~isnan(dataMean_NbyCbyT{1}), 3, 1, 'first');
                tMaxValidByAlignBasisCondition(iA, :, :) = TensorUtils.findNAlongDim(~isnan(dataMean_NbyCbyT{1}), 3, 1, 'last');
                dataValid(iA, :, :) = any(~isnan(dataMean_NbyCbyT{iA}), 3);
            end
            
            if isempty(p.Results.basisNames)
                basisNames = arrayfun(@(i) sprintf('basis%d', i), (1:N)', 'UniformOutput', false);
            else
                basisNames = p.Results.basisNames;
            end
            
            pset = PopulationTrajectorySet();
            pset.dataSourceManual = true;
            pset.alignDescriptorSet = alignDescriptorSet;
            pset.basisNames = basisNames;
            pset.timeUnitName = p.Results.timeUnitName;
            pset.timeUnitsPerSecond = p.Results.timeUnitsPerSecond;
            pset.spikeFilter = NonOverlappingSpikeBinFilter('timeDelta', timeDelta);
            pset.conditionDescriptor = cd;
            pset.tMinForDataMean = tMinDataMean;
            pset.tMaxForDataMean = tMaxDataMean;
            pset.dataMean = dataMean_NbyCbyT;
            pset.dataSem = dataSem;
            pset.dataValid = dataValid;
            pset.dataNTrials = double(dataValid);
            pset.tMinValidByAlignBasisCondition = tMinValidByAlignBasisCondition;
            pset.tMaxValidByAlignBasisCondition = tMaxValidByAlignBasisCondition;
            
            alignSummaryData = cell(1, A);
            for iA = 1:A
                alignSummaryData{iA} = AlignSummary.buildEmptyFromConditionAlignDescriptor(cd, pset.alignDescriptorSet{iA}, tMinDataMean(iA), tMaxDataMean(iA));
            end
            pset.alignSummaryData = alignSummaryData;
            pset.basisAlignSummaryLookup = ones(N, 1);
            pset = pset.initialize();
        end
    end
    
    % methods to list the properties that are involved in specific contexts
    methods(Static)
        function propMetaStruct = filterPropMetaNonEmptyRequired(propMetaStruct)
            
            propMetaStruct = PropertyShapeMeta.filterMetaStruct(propMetaStruct, @(meta, name) meta.getAttrWithDefault('emptyOkay', false), true);
        end
        
        function [toCopy, toCheckNonEmpty] = getPropMetaInGroups(groups)
            toCopy = PopulationTrajectorySet.getPropMetaInGroup(groups);
            toCheckNonEmpty = PopulationTrajectorySetBuilder.filterPropMetaNonEmptyRequired(toCopy);
        end
            
        function [toCopy, toCheckNonEmpty] = getPropMetaSingleTrial(varargin)
            p = inputParser();
            p.addParameter('includeDataSourceInfo', true, @islogical); % used when construction psets with dataSourceManual == true
            p.parse(varargin{:});
            
            groups = {'settings', 'descriptors', 'basisInfo', 'singleTrial', ...
                 'trialAverage', 'randomized', 'noise'}
            if p.Results.includeDataSourceInfo
                groups{end+1} = 'dataSourceInfo';
            end
            
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.getPropMetaInGroups(groups);
        end
        
        function [toCopy, toCheckNonEmpty] = getPropMetaTrialAverage(varargin)
            p = inputParser();
            p.addParameter('includeDiffTrialsNoise', true, @islogical); % this can be slow so we make it optional
%             p.addParameter('includeDataSourceInfo', true, @islogical); % used when construction psets with dataSourceManual == true
            p.parse(varargin{:});
            
            groups = {'settings', 'descriptors', 'basisInfo', ...
                 'trialAverage', 'randomized'}
%             if p.Results.includeDataSourceInfo
%                 groups{end+1} = 'dataSourceInfo';
%             end
            if p.Results.includeDiffTrialsNoise
                groups{end+1} = 'noise';
            end
            
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.getPropMetaInGroups(groups);
        end
    end
    
    % methods to copy fields from pset into a builder
    methods(Static)
        function bld = copySingleTrialFromPopulationTrajectorySet(pset)
            % copy values from population trajectory set, all fields by default
            bld = PopulationTrajectorySetBuilder();
            propMeta = PopulationTrajectorySetBuilder.getPropMetaSingleTrial();
            bld.data = pset.getValuesFromPropMeta(propMeta);
        end

        function bld = copyTrialAveragedOnlyFromPopulationTrajectorySet(pset)
            % copy values from population trajectory set, all fields except
            % single trial data
            bld = PopulationTrajectorySetBuilder();
            propMeta = PopulationTrajectorySetBuilder.getPropMetaTrialAverage();
            bld.data = pset.getValuesFromPropMeta(propMeta);
        end
        
        function bld = copySettingsDescriptorsBasisInfoFromPopulationTrajectorySet(pset)
            bld = PopulationTrajectorySetBuilder();
            propMeta = PopulationTrajectorySetBuilder.getPropMetaInGroups({'settings', 'descriptors', 'basisInfo'});
            bld.data = pset.getValuesFromPropMeta(propMeta);
        end
        
        function bld = copySettingsDescriptorsFromPopulationTrajectorySet(pset)
            bld = PopulationTrajectorySetBuilder();
            propMeta = PopulationTrajectorySetBuilder.getPropMetaInGroups({'settings', 'descriptors'});
            bld.data = pset.getValuesFromPropMeta(propMeta);
        end
    end
    
    methods(Static)
        function psetManual = convertToAutoWithManualDataByTrial(pset, varargin)
            bld = PopulationTrajectorySetBuilder.copySingleTrialFromPopulationTrajectorySet(pset);
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
                bld.data.(fld) = cat(1, vals{:});
            end
            
            nBasesTotal = sum(cellfun(@(pset) pset.nBases, psetCell));
            
            if simultaneous
                bld.data.dataSources = pset.dataSource;
                bld.data.basisDataSourceIdx = ones(nBasesTotal, 1);
            else
                vals = getFn('dataSources');
                bld.data.dataSources = cat(1, vals{:});
                vals = getFn('basisDataSourceIdx');
            
                % offset the data source idx according to the concatenation
                nSources = cellfun(@(pset) pset.nDataSources, psetCell);
                [vals, whichPset] = TensorUtils.catWhich(1, vals{:});
                accum = nSources(1);
                for iP = 2:numel(psetCell)
                    vals(whichPset == iP) = vals(whichPset == iP) + accum;
                    accum = accum + nSources(iP);
                end
                bld.data.basisDataSourceIdx = vals;
            end
            
            vals = getFn('basisDataSourceChannelNames');
            bld.data.basisDataSourceChannelNames = cat(1, vals{:});
            
            if simultaneous
                pset = bld.buildAutoWithSingleTrialData();
            else
                pset = bld.buildAutoWithTrialAveragedData();
            end
        end
    end
    
    % methods to build a pset from the fields that were stored
    methods
        function pset = buildWithFields(bld, propMetaToCopy, varargin)
            p = inputParser();
            p.addParameter('dataInfoManual', false, @islogical);
            p.addParameter('dataByTrialManual', false, @islogical);
            p.addParameter('dataMeanManual', false, @islogical);
            p.parse(varargin{:});
            
            pset = PopulationTrajectorySet();
            
            % first set the properties that determine what is manual and what is automatically computed
            % these settings will dictate where the subsequent properties are stored
            pset.dataInfoManual = p.Results.dataInfoManual;
            pset.dataByTrialManual = p.Results.dataByTrialManual;
            pset.dataMeanManual = p.Results.dataMeanManual;
            
            flds = fieldnames(propMetaToCopy);
            for f = 1:numel(flds)
                if isfield(bld.data, flds{f})
                    pset.(flds{f}) = bld.data.(flds{f});
                end
            end

            pset = pset.initialize();
        end
        
        function pset = buildAuto(bld, varargin)
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.getPropMetaSingleTrial(varargin{:});
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildWithFields(toCopy);
        end
        
        function pset = buildManualDataInfo(bld, varargin)
            [toCopy, toCheckNonEmpty] = PopulationTrajectorySetBuilder.getPropMetaSingleTrial(varargin{:});
            bld.assertNonEmpty(toCheckNonEmpty);
            pset = bld.buildWithFields(toCopy);
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
%             fields = setdiff(fields, PopulationTrajectorySetBuilder.fCanBeEmptyExceptions);
            isEmpty = cellfun(@(fld) isempty(bld.data.(fld)), fields);
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
                if isfield(bld.data, fld)
                    pset.(fld) = bld.data.(fld);
                elseif pset.propCanBeSetCurrently(fld)
                    % warn since this property should theoretically be set 
                    warning('PTS Builder missing property value for %s', fld);
                end
            end 
            
            pset = pset.initialize();
        end
        
        function pset = buildAutoWithFields(bld, fields)
            pset = PopulationTrajectorySet();
            pset.dataSourceManual = false;
            
            for iFld = 1:numel(fields)
                fld = fields{iFld};
                pset.(fld) = bld.data.(fld);
            end 
            
            pset = pset.initialize();
            pset.dataSourceManual = false;
            pset = pset.applyAlignDescriptorSet();
            pset = pset.applyConditionDescriptor();
        end
    end
end
