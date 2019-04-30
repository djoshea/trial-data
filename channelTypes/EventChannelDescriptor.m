classdef EventChannelDescriptor < ChannelDescriptor
    
    properties
        color % default color used when used in a mark or interval 
    end
    
    methods(Access=protected)
        function cd = EventChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.VECTOR;
            cd.fieldIds = {'times'};
            cd = cd.initialize();
        end
    end
    
    methods    
        function cd = initialize(cd)
            % check and repair internal consistency
            cd.warnIfNoArgOut(nargout);
            if isempty(cd.fieldIds)
                cd.fieldIds = {'times'};
            end
            
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = EventChannelImpl(cd);
        end
        
        
        function type = getType(~)
            type = 'event';
        end

        function str = describe(cd)
            str = sprintf('Event %s', cd.name);  
        end

        function cd = inferAttributesFromData(cd, dataCell)
            cd.warnIfNoArgOut(nargout);
            cd.originalDataClassByField = {ChannelImpl.getCellElementClass(dataCell)};
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
        
        function buildEventWithTagFields(name, timeUnits, tagFields)
            error('not yet implemented');
        end
        
        function cls = getSubChannelClass()
            cls = '';
        end
    end

end
