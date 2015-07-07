classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
    end
    
    properties(Dependent)
        timeField
    end
    
    methods(Access=protected)
        function cd = AnalogChannelDescriptor(name, timeField)
            cd = cd@ChannelDescriptor(name);
            if nargin < 2
                timeField = sprintf('%s_time', cd.name);
            end
            
            cd.dataFields = {name, timeField};
            cd.originalDataClassByField = {'double', 'double'};
            cd.elementTypeByField = [cd.VECTOR, cd.VECTOR];
        end
    end
        
    methods
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
            cd.warnIfNoArgOut(nargout);
            dataClass = ChannelDescriptor.getCellElementClass(dataCell);
            timeClass = ChannelDescriptor.getCellElementClass(timeCell);
            cd.originalDataClassByField = {dataClass, timeClass};
            if strcmp(dataClass, 'cell')
                cd.elementTypeByField = [cd.CELL, cd.VECTOR];
            else
                cd.elementTypeByField = [cd.VECTOR, cd.VECTOR]; 
            end
        end
        
        function [cd, dataFieldRenameStruct] = rename(cd, newName)
            cd.warnIfNoArgOut(nargout);
            % also rename _time field if it matches
            oldName = cd.name;
            [cd, dataFieldRenameStruct] = rename@ChannelDescriptor(cd, newName);
            oldTimeField = sprintf('%s_time', oldName);
            if strcmp(cd.dataFields{2}, oldTimeField)
                newTimeField = sprintf('%s_time', newName);
                dataFieldRenameStruct.(cd.dataField{2}) = newTimeField;
                cd.dataFields{2} = newTimeField;
            end
        end 
        
        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (cd.scaleFromLims -> cd.scaleToLims)
            data = convertDataCellOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                scaleFromLow = cd.scaleFromLims(2);
                scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
                scaleToLow = cd.scaleToLims(2);
                scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
                data = cellfun(@(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow, data, 'UniformOutput', false);
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
        
        function cd = buildVectorAnalogFromValues(name, timeField, units, timeUnits, dataCell, timeCell)
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits);
            cd = cd.inferAttributesFromData(dataCell, timeCell);
        end
    end

end
