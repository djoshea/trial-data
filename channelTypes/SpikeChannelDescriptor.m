classdef SpikeChannelDescriptor < ChannelDescriptor
    properties(Dependent)
        hasWaveforms
        hasSortQualityEachTrial
        hasBlankingRegions
        array string
        electrode
        unit
    end

    properties
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

        spikeThreshold = []; % set if known, otherwise this will be estimated from the waveforms
        
        isColumnOfArray = false; % part of spike channel array descriptor
        primaryDataFieldColumnIndex = 1; % which column am I? Should be 1 if not part of array
    end
    
%     properties(SetAccess=protected)
%         arrayManual = '';
%         electrodeManual = [];
%         unitManual = [];
%     end

    methods(Access=protected)
        function cd = SpikeChannelDescriptor(name)
            cd = cd@ChannelDescriptor(name);

            cd = cd.initialize();
        end
    end

    methods
        function cd = initialize(cd)
            cd.dataFields = {cd.name};
            cd.fieldIds = {'spikes'};
            cd.elementTypeByField = cd.VECTOR;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};

            if cd.hasWaveforms
                cd.fieldIds{end+1} = 'waveforms';
                cd.dataFields{end+1} = cd.waveformsField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = cd.waveformsOriginalDataClass;
                cd.unitsByField{end+1} = cd.waveformsUnits;
            end

            if cd.hasSortQualityEachTrial
                cd.fieldIds{end+1} = 'sortQuality';
                cd.dataFields{end+1} = cd.sortQualityEachTrialField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end

            if cd.hasBlankingRegions
                cd.fieldIds{end+1} = 'blankingRegions';
                cd.dataFields{end+1} = cd.blankingRegionsField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end
            
            cd.catAlongFirstDimByField = false(cd.nFields, 1);
            
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = SpikeChannelImpl(cd);
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

%         function cd = setArrayElectrodeUnit(cd, array, electrode, unit)
%             cd.warnIfNoArgOut(nargout);
%             assert(ischar('array'));
%             assert(isnumeric(electrode));
%             assert(isnumeric(unit));
%             cd.electrodeManual = electrode;
%             cd.arrayManual = array;
%             cd.unitManual = unit;
%         end

%         function cd = clearManualArrayElectrode(cd)
%             cd.warnIfNoArgOut(nargout);
%             cd.electrodeManual = [];
%             cd.arrayManual = [];
%         end

        function cd = addWaveformsField(cd, waveField, varargin)
            p = inputParser;
            p.addParameter('time', [], @isvector);
            p.addParameter('units', 'uV', @ischar);
            p.addParameter('scaleFromLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('scaleToLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('dataClass', '', @ischar);
            p.addParameter('waveformsNumChannels', 1, @isscalar);
            p.addParameter('waveformsInfo', [], @(x) true); % used for spatial coordinates or offsets of the recorded waveforms
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
            cd.waveformsNumChannels = p.Results.waveformsNumChannels;
            cd.waveformsInfo = p.Results.waveformsInfo;
            cd = cd.initialize();
        end

        function cd = removeWaveformsField(cd)
            cd.waveformsField = '';
            cd.waveformsUnits = '';
            cd.waveformsScaleFromLims = [];
            cd.waveformsScaleToLims = [];
            cd.waveformsOriginalDataClass = '';
            cd.waveformsTime = [];
            cd.waveformsNumChannels = 0;
            cd.waveformsInfo = [];

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
%             % waveData is either matrix or cell of matrices
%             dtype = cd.accessClassByField{2};
%             if isempty(cd.waveformsScaleFromLims) || isempty(cd.waveformsScaleToLims)
%                 scaleFn = @(x) cast(x, dtype);
%             else
%                 fromDelta = cd.waveformsScaleFromLims(2) - cd.waveformsScaleFromLims(1);
%                 toDelta = cd.waveformsScaleToLims(2) - cd.waveformsScaleToLims(1);
%
%                 % convert to dtype and scale
%                 scaleFn = @(x) (cast(x, dtype) - cd.waveformsScaleFromLims(1)) / fromDelta * toDelta + cd.waveformsScaleToLims(1);
%             end
%             if iscell(waveData)
%                 waveData = cellfun(scaleFn, waveData, 'UniformOutput', false);
%             else
%                 waveData = scaleFn(waveData);
%             end
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

        function name = getNameWithUpdatedArray(cd, array)
            name = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(array, cd.electrode, cd.unit);
        end

        function name = getNameWithUpdatedElectrode(cd, electrode)
            name = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(cd.array, electrode, cd.unit);
        end

        function name = getNameWithUpdatedUnit(cd, unit)
            name = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(cd.array, cd.electrode, unit);
        end

        function array = get.array(cd)
            array = SpikeChannelDescriptor.parseArrayElectrodeUnit(cd.name);
        end

        function elec = get.electrode(cd)
            [~, elec] = SpikeChannelDescriptor.parseArrayElectrodeUnit(cd.name);
        end

        function unit = get.unit(cd)
            [~, ~, unit] = SpikeChannelDescriptor.parseArrayElectrodeUnit(cd.name);
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

        function type = getType(~)
            type = 'spike';
        end

        function str = describe(cd)
            str = sprintf('Unit %s', cd.name);
        end

        function cd = inferAttributesFromData(cd, varargin)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            assert(numel(varargin) == 1, 'Spike Channel descriptor takes exactly 1 data cell');

            cd.originalDataClassByField = {ChannelDescriptor.getCellElementClass(varargin{1})};
            cd.elementTypeByField = cd.VECTOR;
        end
    end

    methods(Static)
        function cd = build(name, varargin)
            p = inputParser();
            p.addParameter('waveformsField', '', @ischar);
            p.addParameter('waveformsTime', [], @(x) isempty(x) || isvector(x));
            p.addParameter('waveformsUnits', 'uV', @ischar);
            p.addParameter('waveformsScaleFromLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('waveformsScaleToLims', [], @(x) isvector(x) || isempty(x));
            p.addParameter('waveformsOriginalDataClass', '', @ischar);
            p.addParameter('waveformsNumChannels', 1, @isscalar);
            p.addParameter('waveformsInfo', [], @(x) true); % used for spatial coordinates or offsets of the recorded waveforms
            p.addParameter('sortQualityEachTrialField', '', @ischar);
            p.addParameter('sortQuality', NaN, @isscalar)
            p.addParameter('blankingRegionsField', '', @ischar);
            p.addParameter('spikeThreshold', NaN, @isscalar);
            p.addParameter('isColumnOfArray', false, @islogical);
            p.addParameter('primaryDataFieldColumnIndex', 1, @isscalar);
            p.parse(varargin{:});
            
            cd = SpikeChannelDescriptor(name);
            if ~isempty(p.Results.waveformsField)
                cd = cd.addWaveformsField(p.Results.waveformsField, 'time', p.Results.waveformsTime, ...
                    'units', p.Results.waveformsUnits, ...
                    'scaleFromLims', p.Results.waveformsScaleFromLims, ...
                    'scaleToLims', p.Results.waveformsScaleToLims, ...
                    'dataClass', p.Results.waveformsOriginalDataClass, ...
                    'waveformsNumChannels', p.Results.waveformsNumChannels, ...
                    'waveformsInfo', p.Results.waveformsInfo);
            end
            
            if ~isempty(p.Results.sortQualityEachTrialField)
                cd = cd.addSortQualityEachTrialField(p.Results.sortQualityEachTrialField);
            end
            cd.sortQuality = p.Results.sortQuality;
            
            if ~isempty(p.Results.blankingRegionsField)
                cd = cd.addBlankingRegionsField(p.Results.blankingRegionsField);
            end
            
            cd.spikeThreshold = p.Results.spikeThreshold;
            
            cd.isColumnOfArray = p.Results.isColumnOfArray;
            cd.primaryDataFieldColumnIndex = p.Results.primaryDataFieldColumnIndex;
        end

%         function cd = buildFromUnitName(name)
%             % attempt to parse the unit name into array electrode# and
%             % unit#
%             cd = SpikeChannelDescriptor(name);
%         end

        function cd = buildFromArrayElectrodeDotUnit(electrodeDotUnit, varargin)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(electrodeDotUnit);
            cd = SpikeChannelDescriptor(name, varargin{:});
        end

        function fld = convertUnitNameToChannelName(unitStr)
            fld = ['unit', strrep(unitStr, '.', '_')];
        end

        function cd = buildFromArrayElectrodeUnit(array, electrode, unit, varargin)
            p  = inputParser();
            p.addParameter('maxElectrode', [], @(x) isempty(x) || isscalar(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            name = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(array, electrode, unit, p.Results.maxElectrodes);
            cd = SpikeChannelDescriptor.buildFromUnitName(name, p.Unmatched);
        end

        function [array, electrode, unit] = parseArrayElectrodeUnit(unitName)
            tokens = regexp(unitName, '^(?<array>[A-Za-z0-9_]*?)(?<electrode>\d+)[_+](?<unit>\d+)$', 'names', 'once');
            if isempty(tokens)
                % try with just array and unit
                tokens = regexp(unitName, '^(?<array>[A-Za-z0-9_]*?)(?<unit>\d+)$', 'names', 'once');

                if isempty(tokens)
                    array = '';
                    electrode = NaN;
                    unit = NaN;
                else
                    array = tokens.array;
                    electrode = NaN;
                    unit =str2double(tokens.unit);
                end
            else
                array = tokens.array;
                electrode =str2double(tokens.electrode);
                unit = str2double(tokens.unit);
            end
        end

        function name = generateNameFromArrayElectrodeUnit(array, electrode, unit, maxElectrode)
            if nargin < 4 || isempty(maxElectrode)
                if isnan(electrode)
                    maxElectrode = max(99, unit);
                else
                    maxElectrode = max(99, electrode);
                end
            end
            nZeros = floor(log10(double(maxElectrode))) + 1;
            
            if isnan(electrode)
                name = sprintf("%s%0*d", array, nZeros, unit);
            else
                name = sprintf("%s%0*d_%d", array, nZeros, electrode, unit);
            end
        end
        
        function names = generateNameListFromArrayElectrodeUnit(arrays, electrodes, units, maxElectrode)
            if nargin < 4
                if ~all(isnan(electrodes))
                    maxElectrode = max(99, max(electrodes));
                else
                    maxElectrode = max(99, max(units)); % no electrode numbers
                end
            end
            
            arrays = string(arrays);
            N = max([numel(arrays) numel(electrodes) numel(units)]);
            if numel(arrays) == 1
                arrays = repmat(arrays, N, 1);
            end
            if isscalar(electrodes)
                electrodes = repmat(electrodes, N, 1);
            end
            if isscalar(units)
                units = repmat(units, N, 1);
            end
               
            names = arrayfun(@(a,e,u) SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(a, e, u, maxElectrode), ...
                makecol(arrays), makecol(electrodes), makecol(units));
        end
        
        function cls = getSubChannelClass()
            cls = '';
        end
    end

end
