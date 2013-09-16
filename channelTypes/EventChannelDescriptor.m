classdef EventChannelDescriptor < ChannelDescriptor
    properties
        tagFields = {}
    end

    methods
        function type = getType(cdesc)
            type = 'event';
        end

        function str = describe(cdesc)
            str = sprintf('Event');  
        end

        function dataFields = getDataFields(cdesc)
            dataFields = {cdesc.name; cdesc.tagFields{:}};
        end

        function cd = EventChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
        
        function cd = inferAttributesFromData(cd, dataCell)
            cd = inferAttributesFromData@ChannelDescriptor(cd, dataCell);
        end
    end
    
    methods(Static)
        function cd = buildScalarEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
            cd.units = timeUnits;
        end
    end

end
