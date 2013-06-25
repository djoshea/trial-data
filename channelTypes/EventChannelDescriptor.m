classdef EventChannelDescriptor < ChannelDescriptor
    methods
        function type = getType(cdesc)
            type = 'event';
        end

        function str = describe(cdesc)
            str = sprintf('Event');  
        end

        function dataFields = getExtraDataFields(cdesc)
            dataFields = {'tags'};
        end

        function cd = EventChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
    end

end
