classdef SpikeChannelDescriptor < ChannelDescriptor
    properties(Dependent)
        hasWaveforms
        hasSortQualityEachTrial
        hasBlankingRegions
        array
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
        
        sortQualityEachTrialField = '';
        
        blankingRegionsField = ''; % refers to a field which conveys times where spikes from this channel are to be considered "unobserved"
        
        sortQuality = NaN; % numeric scalar metric of sort quality
        sortMethod = '';
        
        spikeThreshold = []; % set if known, otherwise this will be estimated from the waveforms
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
            cd.elementTypeByField = cd.VECTOR;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};
            
            if cd.hasWaveforms
                cd.dataFields{end+1} = cd.waveformsField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = cd.waveformsOriginalDataClass;
                cd.unitsByField{end+1} = cd.waveformsUnits;
            end
            
            if cd.hasSortQualityEachTrial
                cd.dataFields{end+1} = cd.sortQualityEachTrialField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end
            
            if cd.hasBlankingRegions
                cd.dataFields{end+1} = cd.blankingRegionsField;
                cd.elementTypeByField(end+1) = cd.NUMERIC;
                cd.originalDataClassByField{end+1} = 'double';
                cd.unitsByField{end+1} = '';
            end
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
            p.addParamValue('time', [], @isvector);
            p.addParamValue('units', 'uV', @ischar);
            p.addParamValue('scaleFromLims', [], @isvector);
            p.addParamValue('scaleToLims', [], @isvector);
            p.addParamValue('dataClass', '', @ischar);
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
        
        
        function waveData = scaleWaveforms(cd, waveData)
            % waveData is either matrix or cell of matrices
            dtype = cd.accessClassByField{2};
            if isempty(cd.waveformsScaleFromLims) || isempty(cd.waveformsScaleToLims)
                scaleFn = @(x) cast(x, dtype);
            else
                fromDelta = cd.waveformsScaleFromLims(2) - cd.waveformsScaleFromLims(1);
                toDelta = cd.waveformsScaleToLims(2) - cd.waveformsScaleToLims(1);
                
                % convert to dtype and scale
                scaleFn = @(x) (cast(x, dtype) - cd.waveformsScaleFromLims(1)) / fromDelta * toDelta + cd.waveformsScaleToLims(1);
            end
            if iscell(waveData)
                waveData = cellfun(scaleFn, waveData, 'UniformOutput', false);
            else
                waveData = scaleFn(waveData);
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
        function cd = build(name)
            cd = SpikeChannelDescriptor(name);
        end
        
        function cd = buildFromUnitName(name)
            % attempt to parse the unit name into array electrode# and
            % unit#
            cd = SpikeChannelDescriptor(name);
        end
        
        function cd = buildFromArrayElectrodeDotUnit(electrodeDotUnit)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(electrodeDotUnit);
            cd = SpikeChannelDescriptor(name);
        end
        
        function fld = convertUnitNameToChannelName(unitStr)
            fld = ['unit', strrep(unitStr, '.', '_')];
        end
        
        function cd = buildFromArrayElectrodeUnit(array, electrode, unit)
            name = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(array, electrode, unit);
            cd = SpikeChannelDescriptor.buildFromUnitName(name);
        end
        
        function [array, electrode, unit] = parseArrayElectrodeUnit(unitName)
            tokens = regexp(unitName, '(?<array>[A-Za-z_]*)(?<electrode>\d+)[_+](?<unit>\d+)', 'names', 'once');
            if isempty(tokens)
                array = '';
                electrode = NaN;
                unit = NaN;
            else
                array = tokens.array;
                electrode =str2double(tokens.electrode);
                unit = str2double(tokens.unit);
            end
        end
        
        function name = generateNameFromArrayElectrodeUnit(array, electrode, unit)
            name = sprintf('%s%02d_%d', array, electrode, unit);
        end
    end
    
end
