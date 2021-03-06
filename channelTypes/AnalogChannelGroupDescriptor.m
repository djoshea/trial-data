classdef AnalogChannelGroupDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
    end
    
    properties(Dependent)
        timeField
        hasScaling
    end

    methods
        function tf = get.hasScaling(cd)
            tf = ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims) && ~isequal(cd.scaleFromLims, cd.scaleToLims);
        end
        
        function cd = withNoScaling(cd)
            cd.warnIfNoArgOut(nargout);
            cd.originalDataClassByField{1} = 'single';
            cd.scaleFromLims = [];
            cd.scaleToLims = [];
        end
    end
    
    methods(Access=protected)
        function cd = AnalogChannelGroupDescriptor(name, timeField)
            cd = cd@ChannelDescriptor(name);
            if nargin < 2
                timeField = sprintf('%s_time', cd.name);
            end
            
            cd.dataFields = {name, timeField};
            cd.originalDataClassByField = {'double', 'double'};
            cd.elementTypeByField = [cd.NUMERIC, cd.VECTOR];
        end
    end
        
    methods
        function cd = initialize(cd)
            cd.warnIfNoArgOut(nargout);
            cd.dataFields = {cd.name, cd.timeField};
        end
        
        % used by trial data when it needs to change field names
        function name = suggestFieldName(cd, fieldIdx)
            if fieldIdx == 1
                name = cd.name;
            elseif fieldIdx == 2
                name = sprintf('%s_time', cd.name);
            else                
                name = sprintf('%s_f%d', fieldIdx);
            end
        end
        
        function f = get.timeField(cd)
            if cd.nFields < 2
                f = '';
            else
                f = cd.dataFields{2};
            end
        end
        
        function type = getType(~)
            type = 'analogGroup';
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
            dataCell = {tdi.trials.(group.signalNames{1})};
                    
            cd.originalDataClassByField = {dataClass, timeClass};
            if strcmp(dataClass, 'cell')
                cd.elementTypeByField = [cd.CELL, cd.VECTOR];
            else
                cd.elementTypeByField = [cd.NUMERIC, cd.VECTOR]; 
            end
        end
        
        function [cd, dataFieldRenameStruct] = rename(cd, newName, renameTime)
            if nargin < 3
                renameTime = true;
            end
            cd.warnIfNoArgOut(nargout);
            
            % also rename _time field if it matches
            oldName = cd.name;
            [cd, dataFieldRenameStruct] = rename@ChannelDescriptor(cd, newName);
            
            if renameTime
                oldTimeField = sprintf('%s_time', oldName);
                if strcmp(cd.dataFields{2}, oldTimeField)
                    newTimeField = sprintf('%s_time', newName);
                    dataFieldRenameStruct.(cd.dataFields{2}) = newTimeField;
                    cd.dataFields{2} = newTimeField;
                end
            end
        end 
        
        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (cd.scaleFromLims -> cd.scaleToLims)
            data = convertDataCellOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims) && ~isempty(data)
                scaleFromLow = cd.scaleFromLims(1);
                scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
                scaleToLow = cd.scaleToLims(1);
                scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
                scaleFn = @(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
                if iscell(data)
                    data = cellfun(scaleFn, data, 'UniformOutput', false);
                else
                    data = scaleFn(data);
                end
            end
        end
        
        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            data = convertDataSingleOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                scaleFromLow = cd.scaleFromLims(2);
                scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
                scaleToLow = cd.scaleToLims(2);
                scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
                data = (data-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
            end
        end
        
        function data = convertAccessDataCellToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                scaleFromLow = cd.scaleFromLims(2);
                scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
                scaleToLow = cd.scaleToLims(2);
                scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
                data = cellfun(@(d) (d-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow, data, 'UniformOutput', false);
            end
            data = convertAccessDataCellToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end
        
        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                scaleFromLow = cd.scaleFromLims(2);
                scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
                scaleToLow = cd.scaleToLims(2);
                scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
                data = (data-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow;
            end
            data = convertAccessDataSingleToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end
        
        function cdIndividual = buildIndividualSubChannel(cd, name, index, units)
            if nargin < 4
                units = cd.unitsByField{1};
            end
            cdIndividual = AnalogChannelDescriptor.buildSharedMatrixColumnAnalog(name, cd.name, index, cd.timeField, units, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2});
        end
    end
    
     methods(Static)
        function cd = buildAnalogGroup(name, timeField, units, timeUnits, varargin)
            p = inputParser();
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'AnalogChannelGroupDescriptor')); % used by subclasses
            p.addParameter('scaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('scaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('dataClass', 'double', @ischar);
            p.addParameter('timeClass', 'double', @ischar);
            p.parse(varargin{:});
            
            if isempty(p.Results.channelDescriptor)
                cd = AnalogChannelGroupDescriptor(name, timeField);
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
        
        function cd = buildAnalogGroupFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            p = inputParser();
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'AnalogChannelDescriptor')); % used by subclasses
            p.parse(varargin{:});
            
            if isempty(p.Results.channelDescriptor)
                cd = AnalogChannelGroupDescriptor.buildAnalogGroup(name, timeField, units, timeUnits);
            else
                cd = p.Results.channelDescriptor;
            end
            cd = cd.inferAttributesFromData(dataCell, timeCell);
            
            if nargin > 3
                cd.unitsByField = {units, timeUnits};
            elseif nargin > 2
                cd.unitsByField = {units, ''};
            end
            
            cd.scaleFromLims = p.Results.scaleFromLims;
            cd.scaleToLims = p.Results.scaleToLims;
            
            cd = cd.initialize();
        end
    end

end
