classdef EventChannelDescriptor < ChannelDescriptor
    methods(Access=protected)
        function cd = EventChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.VECTOR;
        end
    end
    
    methods    
        function type = getType(~)
            type = 'event';
        end

        function str = describe(cd)
            str = sprintf('Event %s', cd.name);  
        end

        function cd = inferAttributesFromData(cd, dataCell)
            cd.warnIfNoArgOut(nargout);
            cd.originalDataClassByField = {ChannelDescriptor.getCellElementClass(dataCell)};
        end
    end
    
    methods(Static)
        function cd = buildSingleEvent(name, timeUnits)
            cd = EventChannelDescriptor(name);
%             cd.elementTypeByField = cd.SCALAR;
%           Switching this to be vector always for consistency
            cd.elementTypeByField = cd.VECTOR;
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
