classdef AnalogChannelDescriptor < ChannelDescriptor
    properties(Dependent)
        timeField
    end
    
    methods
        function cd = AnalogChannelDescriptor(name, timeField)
            cd = cd@ChannelDescriptor(name);
            if nargin < 2
                timeField = sprintf('%s_time', cd.name);
            end
            
            cd.dataFields = {name, timeField};
            cd.originalDataClassByField = {'double', 'double'};
            cd.elementTypeByField = [cd.VECTOR, cd.VECTOR];
        end
        
        function f = get.timeField(cd)
            if cd.nFields < 2
                f = '';
            else
                f = cd.dataFields{2};
            end
        end
        
        function type = getType(~)
            type = 'analog';
        end

        function fields = getDataFields(cd)
            fields = {cd.name, cd.timeField};
        end

        function str = describe(cd)
            if isempty(cd.unitsPrimary)
                str = cd.name;
            else
                str = sprintf('%s (%s)', cd.name, cd.unitsPrimary); 
            end
        end

        function cd = inferAttributesFromData(cd, dataCell, timeCell)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');
            dataClass = ChannelDescriptor.getCellElementClass(dataCell);
            timeClass = ChannelDescriptor.getCellElementClass(timeCell);
            cd.originalDataClassByField = {dataClass, timeClass};
            if strcmp(dataClass, 'cell')
                cd.elementTypeByField = [cd.CELL, cd.VECTOR];
            else
                cd.elementTypeByField = [cd.VECTOR, cd.VECTOR]; 
            end
        end
    end
    
     methods(Static)
        function cd = buildVectorAnalog(name, timeField, units, timeUnits)
            cd = AnalogChannelDescriptor(name, timeField);
            if nargin > 3
                cd.unitsByField = {units, timeUnits};
            elseif nargin > 2
                cd.unitsByField = {units, ''};
            end
        end
    end

end
