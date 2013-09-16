classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        timeField
    end

    methods
        % default to name_time if not manually specified
        function field = get.timeField(cdesc)
            if isempty(cdesc.timeField)
                if isempty(cdesc.name)
                    field = '';
                else
                    field = sprintf('%s_time', cdesc.name);
                end
            else
                field = cdesc.timeField;
            end
        end

        function type = getType(cdesc)
            type = 'analog';
        end

        function fields = getDataFields(cdesc)
            fields = {cdesc.name, cdesc.timeField};
        end

        function str = describe(cdesc)
            str = sprintf('Analog (%s)', cdesc.name, cdesc.units);  
        end

        function cd = AnalogChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
    end

end
