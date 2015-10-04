classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
        
        isColumnOfSharedMatrix = false; % this field shares a data field with other channels
        primaryDataFieldColumnIndex = 1; % which column am I?    
    end
    
    properties(Dependent)
        % set this to make the data field be different from .name (the default).
        % this is only allowed for shared matrix analog signals
        primaryDataField
        
        timeField
    end
    
    properties(Hidden)
        primaryDataFieldManual 
    end
    
    methods
        function f = get.primaryDataField(cd)
            if ~cd.isColumnOfSharedMatrix || isempty(cd.primaryDataFieldManual)
                f = cd.name;
            else
                f = cd.primaryDataFieldManual;
            end
        end
        
        function cd = set.primaryDataField(cd, f)
            assert(cd.isColumnOfSharedMatrix, 'Primary data field cannot be set unless isColumnOfSharedMatrix is set to true.');
            assert(ischar(f) && isvector(f), 'Field name must be string');
            cd.primaryDataFieldManual = f;
            cd = cd.initialize();
        end
        
        function cd = separateFromColumnOfSharedMatrix(cd)
            % transform this cd so that it's not a column of a shared
            % matrix anymore
            cd.warnIfNoArgOut(nargout);
            cd.isColumnOfSharedMatrix = false;
            cd.primaryDataFieldColumnIndex = 1;
            cd.primaryDataFieldManual = '';
            cd = cd.initialize();
        end
        
        function cd = withNoScaling(cd)
            cd.warnIfNoArgOut(nargout);
            cd.scaleFromLims = [];
            cd.scaleToLims = [];
        end
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
        function cd = initialize(cd)
            cd.warnIfNoArgOut(nargout);
            cd.dataFields = {cd.primaryDataField, cd.timeField};
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
            fields = {cd.primaryDataField, cd.timeField};
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
        function cd = buildVectorAnalog(name, timeField, units, timeUnits, varargin)
            p = inputParser();
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'AnalogChannelDescriptor')); % used by subclasses
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);
            p.addParameter('dataClass', 'double', @ischar);
            p.addParameter('timeClass', 'double', @ischar);
            p.parse(varargin{:});
            
            if isempty(p.Results.channelDescriptor)
                cd = AnalogChannelDescriptor(name, timeField);
            else
                cd = p.Results.channelDescriptor;
            end
            
            if nargin > 3
                cd.unitsByField = {units, timeUnits};
            elseif nargin > 2
                cd.unitsByField = {units, ''};
            end
            
            % set the scale on access limits
            cd.scaleFromLims = p.Results.scaleFromLims;
            cd.scaleToLims = p.Results.scaleToLims;
            
            % set the data and time class appropriately
            cd.originalDataClassByField = {p.Results.dataClass, p.Results.timeClass};
            cd = cd.initialize();
        end
        
        % dataClass and timeClass are the in-memory representation of the
        % data values and time values (e.g. 'double', 'int32', etc)
        function cd = buildSharedMatrixColumnAnalog(name, dataFieldName, dataFieldColumnIndex, ...
                timeField, units, timeUnits, varargin)
            assert(ischar(dataFieldName));
            assert(isscalar(dataFieldColumnIndex));
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits, varargin{:});
            cd.isColumnOfSharedMatrix = true;
            cd.primaryDataField = dataFieldName;
            cd.primaryDataFieldColumnIndex = dataFieldColumnIndex;
            cd = cd.initialize();
        end
        
        function cd = buildVectorAnalogFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            p = inputParser();
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);
            p.parse(varargin{:});
            
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits);
            cd = cd.inferAttributesFromData(dataCell, timeCell);
            
            cd.scaleFromLims = p.Results.scaleFromLims;
            cd.scaleToLims = p.Results.scaleToLims;
            cd = cd.initialize();
        end
    end

end
