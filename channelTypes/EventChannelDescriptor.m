classdef EventChannelDescriptor < ChannelDescriptor
    methods
        function cd = EventChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.VECTOR;
        end
        
        function type = getType(~)
            type = 'event';
        end

        function str = describe(cd)
            str = sprintf('%s', cd.name);  
        end

        function cd = inferAttributesFromData(cd, dataCell)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');
            cd.originalDataClassByField = {ChannelDescriptor.getCellElementClass(dataCell)};
        end
    end
    
    methods(Static)
        function cd = buildSingleEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
            cd.elementTypeByField = cd.SCALAR;
            if nargin > 1
                cd.unitsByField = {timeUnits};
            end
        end
        
        function cd = buildMultipleEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
            cd.elementTypeByField = cd.VECTOR;
            if nargin > 1
                cd.unitsByField = {timeUnits};
            end
        end
    end

end
