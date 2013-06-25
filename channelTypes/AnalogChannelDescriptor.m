classdef AnalogChannelDescriptor < ChannelDescriptor
    methods
        function type = getType(cdesc)
            type = 'analog';
        end

        function str = describe(cdesc)
            str = sprintf('Analog (%s)', cdesc.name, cdesc.units);  
        end

        function dataFields = getExtraDataFields(cdesc)
            dataFields = {'time'};
        end

        function cd = AnalogChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
    end

end
