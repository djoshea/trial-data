classdef SpikeChannelDescriptor < ChannelDescriptor
    properties
        waveformsField
        
        unit
        electrode
        
        quality
    end

    methods
        function cd = SpikeChannelDescriptor(name)
            cd = cd@ChannelDescriptor(name); 
        end
        
        function type = getType(cdesc)
            type = 'spike';
        end

        function str = describe(cdesc)
            str = sprintf('Spike %s', cdesc.name);  
        end

        function dataFields = getDataFields(cd)
            dataFields = {cd.name};
            if ~isempty(cd.waveformsField)
                dataFields{end+1} = cd.waveformsField;
            end
        end
        
        function cd = inferAttributesFromData(cd, dataCell)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            cd.dfd = NumericVectorField();
            cd.storageDataClass = 'double';
        end
    end
    
    methods(Static)
        function fld = convertUnitNameToChannelName(unitStr)
            fld = ['unit', strrep(unitStr, '.', '_')];
        end
        
        function fld = convertChannelNameToUnitName(ch)
            info = regexp(ch, 'unit(?<channel>\d+)_(?<unit>\d+)', 'names', 'once');
            if isempty(info)
                error('Could not parse channel name %s as unit', ch);
            else
                if isempty(info.unit)
                    fld = sprintf('%s', info.channel);
                else
                    fld = sprintf('%s.%s', info.channel, info.unit);
                end
            end
        end
        
        function cd = buildChannelDotUnit(unitStr)
            name = SpikeChannelDescriptor.convertUnitNameToChannelName(unitStr);
            cd = SpikeChannelDescriptor(name);
            [ch.electrode, cd.unit] = parseUnitName(unitStr);
        end
    end

end
