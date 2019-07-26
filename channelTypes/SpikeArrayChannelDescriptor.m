classdef SpikeArrayChannelDescriptor < ChannelDescriptor
    properties(Dependent)
        hasWaveforms
        hasSortQualityEachTrial
        hasBlankingRegions
        nChannels
        nElectrodes % number of unique electrodes
        array string

        subChannelNames
    end

    properties
        % nChannels x  1 arrays corresponding to each channel
        electrodes
        units

        waveformsField = '';
        waveformsUnits = '';
        waveformsScaleFromLims = [];
        waveformsScaleToLims = [];
        waveformsOriginalDataClass = '';
        waveformsTime = []; % common time vector to be shared for ALL waveforms for this channel

        waveformsNumChannels = 1;
        waveformsInfo % used for spatial coordinates or offsets of the recorded waveforms
       
        sortQualityEachTrialField = '';

        blankingRegionsField = ''; % refers to a field which conveys times where spikes from this channel are to be considered "unobserved"

        sortQuality = NaN; % numeric scalar metric of sort quality
        sortMethod = '';

        spikeThreshold = NaN; % set if known, otherwise this will be estimated from the waveforms
        
        isColumnOfSharedMatrix = false; % this field shares a data field with other channels
        primaryDataFieldColumnIndex = 1; % which column am I?   
    end

    methods(Access=protected)
        function cd = SpikeArrayChannelDescriptor(name, electrodes, units)
            cd = cd@ChannelDescriptor(name);
            cd.electrodes = makecol(electrodes);
            cd.units = makecol(units);
            
            assert(~any(isnan(cd.electrodes) & isnan(cd.units)));
            cd = cd.initialize();
        end
    end

    methods % built in channel descriptor functions
        function cd = initialize(cd)
            cd.fieldIds = {'spikes'};
            cd.dataFields = {cd.name};
            cd.elementTypeByField = cd.CELL;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};

            if cd.hasWaveforms
                cd.fieldIds{end+1} = 'waveforms';
                cd.dataFields{end+1} = cd.waveformsField;
                cd.elementTypeByField(end+1) = cd.CELL;
                cd.originalDataClassByField{end+1} = cd.waveformsOriginalDataClass;
                cd.unitsByField{end+1} = cd.waveformsUnits;
            end

            if cd.hasSortQualityEachTrial
                cd.fieldIds{end+1} = 'sortQuality';
                cd.dataFields{end+1} = cd.sortQualityEachTrialField;
                cd.elementTypeByField(end+1) = cd.CELL;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end

            if cd.hasBlankingRegions
                cd.fieldIds{end+1} = 'blankingRegions';
                cd.dataFields{end+1} = cd.blankingRegionsField;
                cd.elementTypeByField(end+1) = cd.CELL;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end
            
            cd.catAlongFirstDimByField = false(cd.nFields, 1);
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = SpikeArrayChannelImpl(cd);
        end

        % used by trial data when it needs to change field names
        function name = suggestFieldName(cd, fieldIdx)
            suggest = {cd.name};
            if cd.hasWaveforms
                suggest{end+1} = sprintf('%s_waveforms', cd.name);
            end
            if cd.hasSortQualityEachTrial
                suggest{end+1} = sprintf('%s_sortQualityByTrial', cd.name);
            end
            if cd.hasBlankingRegions
                suggest{end+1} = sprintf('%s_blankingRegions', cd.name);
            end

            if fieldIdx <= numel(suggest)
                name = suggest{fieldIdx};
            else
                name = sprintf('%s_f%d', fieldIdx);
            end
        end

        function cd = addWaveformsField(cd, waveField, varargin)
            p = inputParser;
            p.addParameter('time', [], @isvector);
            p.addParameter('units', 'uV', @ischar);
            p.addParameter('scaleFromLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('scaleToLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('dataClass', '', @ischar);
            p.parse(varargin{:});

            if nargin < 2 || isempty(waveField)
                waveField = sprintf('%s_waveforms', cd.name);
            end

            cd.waveformsField = waveField;
            cd.waveformsTime = p.Results.time;
            cd.waveformsUnits = p.Results.units;
            cd.waveformsOriginalDataClass = p.Results.dataClass;
            cd.waveformsScaleFromLims = p.Results.scaleFromLims;
            cd.waveformsScaleToLims = p.Results.scaleToLims;
            cd = cd.initialize();
        end

        function cd = removeWaveformsField(cd)
            cd.waveformsField = '';
            cd.waveformsUnits = '';
            cd.waveformsScaleFromLims = [];
            cd.waveformsScaleToLims = [];
            cd.waveformsOriginalDataClass = '';
            cd.waveformsTime = [];

            cd = cd.initialize();
        end

        function cd = addSortQualityEachTrialField(cd, field)
            cd.warnIfNoArgOut(nargout);
            if nargin < 2 || isempty(field)
                field = sprintf('%s_sortQualityByTrial', cd.name);
            end
            cd.sortQualityEachTrialField = field;
            cd = cd.initialize();
        end

        function cd = addBlankingRegionsField(cd, field)
            cd.warnIfNoArgOut(nargout);
            if nargin < 2 || isempty(field)
                field = sprintf('%s_blankingRegions', cd.name);
            end
            cd.blankingRegionsField = field;
            cd = cd.initialize();
        end

        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (cd.scaleFromLims -> cd.scaleToLims)
            data = convertDataCellOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelDescriptor.scaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
        end

        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            data = convertDataSingleOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelDescriptor.scaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
        end

        function data = convertAccessDataCellToMemory(cd, fieldIdx, data)
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelDescriptor.unscaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
            data = convertAccessDataCellToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end

        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelDescriptor.unscaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
            data = convertAccessDataSingleToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end

        function waveData = scaleWaveforms(cd, waveData)
            if iscell(waveData)
                waveData = cd.convertDataCellOnAccess(2, waveData);
            else
                waveData = cd.convertDataSingleOnAccess(2, waveData);
            end
        end

        function waveData = unscaleWaveforms(cd, waveData)
            if iscell(waveData)
                waveData = cd.convertAccessDataCellToMemory(2, waveData);
            else
                waveData = cd.convertAccessDataSingleToMemory(2, waveData);
            end
        end

        function [cd, dataFieldRenameStruct] = rename(cd, newName)
            cd.warnIfNoArgOut(nargout);
            % also rename _waveforms field if it matches
            oldName = cd.name;
            [cd, dataFieldRenameStruct] = rename@ChannelDescriptor(cd, newName);
            if ~isempty(cd.waveformsField)
                oldWave = cd.waveformsField;
                if strcmp(oldWave, sprintf('%s_waveforms', oldName))
                    newWave = sprintf('%s_waveforms', newName);
                    dataFieldRenameStruct.(oldWave) = newWave;
                    cd.waveformsField = newWave;
                    cd = cd.initialize();
                end
            end
        end

        function array = get.array(cd)
            array = cd.name;
        end

        function tf = get.hasWaveforms(cd)
            tf = ~isempty(cd.waveformsField);
        end

        function tf = get.hasSortQualityEachTrial(cd)
            tf = ~isempty(cd.sortQualityEachTrialField);
        end

        function tf = get.hasBlankingRegions(cd)
            tf = ~isempty(cd.blankingRegionsField);
        end
        
        function n = get.nElectrodes(cd)
            n = numel(unique(cd.electrodes));
        end

        function type = getType(~)
            type = 'spike';
        end

        function str = describe(cd)
            str = sprintf('Array %s (%d channels)', cd.name, cd.nChannels);
        end

        function n = get.nChannels(cd)
            n = numel(cd.electrodes);
        end

        function subChannelNames = get.subChannelNames(cd)
            subChannelNames = cd.getSubChannelNames();
        end

        function cd = inferAttributesFromData(cd, varargin)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            assert(numel(varargin) == 1, 'Spike Channel descriptor takes exactly 1 data cell');

            cd.originalDataClassByField = {ChannelDescriptor.getCellElementClass(varargin{1})};
            cd.elementTypeByField = cd.VECTOR;
        end
    end

    methods
        function cd = buildSubChannelDescriptor(cd, nameOrIdx) 
            if isnumeric(nameOrIdx)
                index = nameOrIdx;
            else
                [tf, index] = ismember(nameOrIdx, cd.subChannelNames);
                assert(all(tf), 'Channel not found');
            end
            cd = cd.buildIndividualSubChannel(index);
        end
        
        function cds = buildIndividualSubChannel(cd, uidx)
            names = SpikeChannelDescriptor.generateNameListFromArrayElectrodeUnit(...
                cd.array, cd.electrodes(uidx), cd.units(uidx));

            if isscalar(cd.spikeThreshold)
                thresh = repmat(cd.spikeThreshold, numel(names));
            else
                thresh = cd.spikeThreshold(uidx);
            end
            
            args = {...
                'waveformsField', cd.waveformsField, ...
                'waveformsTime', cd.waveformsTime, ...
                'waveformsUnits', cd.waveformsUnits, ...
                'waveformsScaleFromLims', cd.waveformsScaleFromLims, ...
                'waveformsScaleToLims', cd.waveformsScaleToLims, ...
                'waveformsOriginalDataClass', cd.waveformsOriginalDataClass, ...
                'waveformsNumChannels', cd.waveformsNumChannels, ...
                'waveformsInfo', cd.waveformsInfo, ...
                'sortQualityEachTrialField', cd.sortQualityEachTrialField, ...
                'sortQuality', cd.sortQuality, ...
                'blankingRegionsField', cd.blankingRegionsField, ...
                'isColumnOfArray', true};
                
            for i = numel(names):-1:1
                cds(i) = SpikeChannelDescriptor.build(names{i}, ...
                    'primaryDataFieldColumnIndex', uidx(i), ...
                    'spikeThreshold', thresh(i), args{:});
            end
        end
        
        function cds = buildIndividualSubChannelByElectrodeUnit(cd, elec, unit)
            uidx = nan(numel(elec), 1);
            for iE = 1:numel(elec)
                thisIdx = find(cd.electrodes == elec(iE) & cd.units == unit(iE));
                if isempty(thisIdx)
                    error('Array %s does not have electrode %d, unit %d', cd.array, elec(iE), unit(iE));
                end
                if numel(thisIdx) > 1
                    error('Array %s has multiple channels with electrode %d, unit %d', cd.array, elec(iE), unit(iE));
                end
                
                uidx(iE) = thisIdx;
            end
            
            cds = cd.buildIndividualSubChannel(uidx);
        end
        
        function ch_idx = findSubChannelByElectrodeUnit(cd, elec, unit)
            ch_idx = nan(size(elec));
            for iE = 1:numel(elec)
                if isnan(elec(iE))
                    thisIdx = find(isnan(cd.electrodes) & cd.units == unit(iE));
                else
                    thisIdx = find(cd.electrodes == elec(iE) & cd.units == unit(iE));
                end
                if ~isempty(thisIdx)
                    ch_idx(iE) = thisIdx;
                end
            end
        end

        function [names, chidx] = listNamedSubChannels(cd)
            names = SpikeChannelDescriptor.generateNameListFromArrayElectrodeUnit(cd.array, cd.electrodes, cd.units);
            chidx = (1:numel(names))';
        end

        function [tf, idx] = hasSubChannel(cd, names)
            if ischar(names), names = {names}; end
            N = numel(names);
            tf = false(N, 1);
            idx = nan(N, 1);
            for iN = 1:N
                name = names{iN};
                if contains(name, '(')
                    % paren indexed channel
                    [array, chidx] = ChannelDescriptor.parseIndexedChannelName(name); %#ok<*PROPLC>
                    tf(iN) = strcmp(array, cd.name) && TrialDataUtilities.Data.indexInRange(chidx, cd.nChannels);
                else
                    % arrayElec_Unit
                    [array, elec, unit] = SpikeChannelDescriptor.parseArrayElectrodeUnit(name);
                    chidx = cd.findSubChannelByElectrodeUnit(elec, unit);
                    tf(iN) = strcmp(array, cd.name) && TrialDataUtilities.Data.indexInRange(chidx, cd.nChannels);
                end
            end
        end
    end

    methods(Static)
        function cd = build(name, electrodes, units)
            cd = SpikeArrayChannelDescriptor(name, electrodes, units);
        end
        
        function cd = buildFromSubChNames(subChNames)
            [array, elec, unit] = arrayfun(@(subname) SpikeChannelDescriptor.parseArrayElectrodeUnit(subname), string(subChNames));
            name = array(1);
            assert(all(array == name), 'Array names must be identical');
            
            cd = SpikeArrayChannelDescriptor.build(name, elec, unit);
        end
        
        function cls = getSubChannelClass()
            cls = 'SpikeChannelDescriptor';
        end
    end
end
