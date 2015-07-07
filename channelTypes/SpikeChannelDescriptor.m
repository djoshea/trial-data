classdef SpikeChannelDescriptor < ChannelDescriptor
    properties
        waveformsTime = []; % common time vector to be shared for ALL waveforms for this channel
        quality = NaN;
        sortMode = NaN;
    end
    
    properties(SetAccess=protected)
        waveformsField = '';
    end
    
    properties(Dependent)
        hasWaveforms
        unitStr
        unit
        electrode
    end
    
    properties(Constant)
        SORT_THRESHOLD = 1;
        SORT_ONLINE = 2;
        SORT_MANUAL = 3;
    end

    methods(Access=protected)
        function cd = SpikeChannelDescriptor(name)
            cd = cd@ChannelDescriptor(name); 
            cd.dataFields = {cd.name};
            cd.elementTypeByField = cd.VECTOR;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};
        end
    end
    
    methods
        function cd = addWaveformsField(cd, waveField, varargin)
            p = inputParser;
            p.addParamValue('time', [], @isvector);
            p.parse(varargin{:});
            
            if nargin < 2 || isempty(waveField)
                waveField = sprintf('%s_waveforms', cd.name);
            end
            
            cd.warnIfNoArgOut(nargout);
            cd.dataFields = {cd.name; waveField};
            cd.waveformsField = waveField;
            cd.elementTypeByField = [cd.VECTOR; cd.NUMERIC];
            cd.originalDataClassByField = {'', ''};
            cd.unitsByField = {'', 'mV'};
            cd.waveformsTime = p.Results.time;
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
                    cd.dataFields{2} = newWave;
                end
            end
        end 
        
        function u = get.unitStr(cd)
            if isempty(cd.name)
                u = '';
            else
                u = SpikeChannelDescriptor.convertChannelNameToUnitName(cd.name);
            end
        end
        
        function e = get.electrode(cd)
            [e, ~] = SpikeChannelDescriptor.convertChannelNameToElectrodeUnit(cd.name);
        end
        
        function tf = get.hasWaveforms(cd)
            tf = ~isempty(cd.waveformsField);
        end
        
        function u = get.unit(cd)
            [~, u] = SpikeChannelDescriptor.convertChannelNameToElectrodeUnit(cd.name);
        end
        
        function type = getType(~)
            type = 'spike';
        end

        function str = describe(cd)
            str = sprintf('Unit %s', cd.name);  
        end

        function dataFields = getDataFields(cd)
            dataFields = {cd.name};
            if ~isempty(cd.waveformsField)
                dataFields{end+1} = cd.waveformsField;
            end
        end

        function cd = inferAttributesFromData(cd, varargin)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');
            
            assert(numel(varargin) == 1, 'Spike Channel descriptor takes exactly 1 data cell');
            
            cd.originalDataClassByField = {ChannelDescriptor.getCellElementClass(varargin{1})};
            cd.elementTypeByField = cd.VECTOR;
        end
    end
    
    methods(Static)
        function cd = buildFromUnitName(name)
            cd = SpikeChannelDescriptor(name);
        end

        function cd = buildFromUnitStr(unitName)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitName);
            cd = SpikeChannelDescriptor(name);
        end
            
        function fld = convertUnitNameToChannelName(unitStr)
            fld = ['unit', strrep(unitStr, '.', '_')];
        end
        
        function [unitStr, valid] = convertChannelNameToUnitName(ch)
            info = regexp(ch, 'unit(?<channel>\d+)_(?<unit>\d+)', 'names', 'once');
            if isempty(info)
                unitStr = '';
                valid = false;
                %error('Could not parse channel name %s as unit', ch);
            else
                if isempty(info.unit)
                    unitStr = sprintf('%s', info.channel);
                else
                    unitStr = sprintf('%s.%s', info.channel, info.unit);
                end
                valid = true;
            end
        end
        
        function [electrode, unit] = convertUnitNameToElectrodeUnit(unitName)
            tokens = regexp(unitName, '(?<electrode>\d+)\.(?<unit>\d+)', 'names');
            electrode =str2double(tokens.electrode);
            unit = str2double(tokens.unit);
        end
        
                
        function [electrode, unit] = convertChannelNameToElectrodeUnit(ch)
            unitName = SpikeChannelDescriptor.convertChannelNameToUnitName(ch);
            [electrode, unit] = SpikeChannelDescriptor.convertUnitNameToElectrodeUnit(unitName);
        end
    end

end
