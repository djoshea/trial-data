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
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            cd.dfd = NumericVectorField();
            cd.storageDataClass = 'double';
        end
    end
    
    methods(Static)
        function cd = buildSingleEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
            cd.dfd = ScalarField();
            if nargin > 1
                cd.units = timeUnits;
            end
        end
        
        function cd = buildMultipleEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
            cd.dfd = NumericVectorField();
            if nargin > 1
                cd.units = timeUnits;
            end
        end
    end

end
