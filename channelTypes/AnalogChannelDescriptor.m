classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
        
        isColumnOfSharedMatrix = false; % this field shares a data field with other channels
        primaryDataFieldColumnIndex = 1; % which column am I?
        
        % if this channel just a virtually transformed version of a
        % different channel, what is the other channel's name?
        transformChannelNames string = [];
        
        % function handle to apply to other channel's data
        % function should be prepared to take arbitrary arguments using
        % varargin
        % but the suggested signature is:
        % outData = fn(inData, requestedColIdx, varargin)
        transformFn = []; 
        
        % 'simple': out_one_trial = fn(in_one_trial), no slicing, one
        % transformChannel
        % 'manual': channel does all work and accepts all data at once
        transformFnMode = '';
    end
    
    properties(Dependent)
        % set this to make the data field be different from .name (the default).
        % this is only allowed for shared matrix analog signals
        primaryDataField
        dataClass
        timeClass
        
        timeField
        timeUnits
        
        hasScaling
        
        isTransform
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
        
        function tf = get.isTransform(cd)
            tf = ~isempty(cd.transformChannelNames);
        end
        
        function cd = updateGroup(cd, newGroupDescriptor)
            cd.warnIfNoArgOut(nargout);
            assert(isa(newGroupDescriptor, 'AnalogChannelGroupDescriptor'));
            assert(cd.isColumnOfSharedMatrix)
            cd.primaryDataField = newGroupDescriptor.name;
            cd.timeField = newGroupDescriptor.timeField;
        end  
        
        function grp = getGroupName(cd, ~)
            if cd.isColumnOfSharedMatrix
                grp = cd.primaryDataField;
            else
                grp = '';
            end
        end
        
        function cd = set.primaryDataField(cd, f)
            assert(cd.isColumnOfSharedMatrix, 'Primary data field cannot be set unless isColumnOfSharedMatrix is set to true.');
            assert(ischar(f) && isvector(f), 'Field name must be string');
            cd.primaryDataFieldManual = f;
            cd = cd.initialize();
        end
        
        function cd = set.timeField(cd, f)
            cd.dataFields{2} = f;
        end
        
        function c = get.dataClass(cd)
            c = cd.originalDataClassByField{1};
        end
        
        function c = get.timeClass(cd)
            c = cd.originalDataClassByField{2};
        end
        
        function tf = get.hasScaling(cd)
            tf = ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims) && ~isequal(cd.scaleFromLims, cd.scaleToLims);
        end
        
        function cd = separateFromColumnOfSharedMatrix(cd, newTimeField)
            % transform this cd so that it's not a column of a shared
            % matrix anymore
            cd.warnIfNoArgOut(nargout);
            cd.isColumnOfSharedMatrix = false;
            cd.elementTypeByField(1) = cd.VECTOR; 
            cd.primaryDataFieldColumnIndex = 1;
            cd.primaryDataFieldManual = '';
            
            if nargin > 1
                cd.dataFields{2} = newTimeField;
            end
            cd = cd.initialize();
        end
        
        function cd = withNoScaling(cd)
            cd.warnIfNoArgOut(nargout);
            cd.originalDataClassByField{1} = 'single';
            cd.scaleFromLims = [];
            cd.scaleToLims = [];
        end
        
        function cdGroup = buildGroupChannelDescriptor(cd, varargin)
            cdGroup = AnalogChannelGroupDescriptor.buildAnalogGroup(cd.primaryDataField, cd.timeField, ...
                cd.unitsByField{1}, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2}, varargin{:});
            if isempty(cdGroup.sampleSize)
                cdGroup.sampleSize = cd.primaryDataFieldColumnIndex;
            else
                cdGroup.sampleSize = max(cdGroup.sampleSize, cd.primaryDataFieldColumnIndex);
            end
               
            cdGroup = cdGroup.setSubChannelInfo(cd.primaryDataFieldColumnIndex, cd.name, cd.unitsPrimary);
        end
    end
    
    methods(Access=protected)
        function cd = AnalogChannelDescriptor(name, timeField, varargin)
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
        
        function f = get.timeUnits(cd)
            if cd.nFields < 2
                f = '';
            else
                f = cd.unitsByField{2};
            end
        end
        
        function type = getType(~)
            type = 'analog';
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
                dataFieldRenameStruct.(cd.dataFields{2}) = newTimeField;
                cd.dataFields{2} = newTimeField;
            end
        end 
        
        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (cd.scaleFromLims -> cd.scaleToLims)
            data = convertDataCellOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelDescriptor.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = cellfun(@(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow, data, 'UniformOutput', false);
            end
        end
        
        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            data = convertDataSingleOnAccess@ChannelDescriptor(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelDescriptor.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = (data-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
            end
        end
        
        function data = convertAccessDataCellToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelDescriptor.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = cellfun(@(d) (d-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow, data, 'UniformOutput', false);
            end
            data = convertAccessDataCellToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end
        
        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelDescriptor.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = (data-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow;
            end
            data = convertAccessDataSingleToMemory@ChannelDescriptor(cd, fieldIdx, data);
        end
        
        function [dataCell, timeCell] = computeTransformDataRaw(cd, td, varargin)
            % if shared column, use the slice arg to subselect the
            % appropriate column
            if cd.isColumnOfSharedMatrix
                sliceArgs = {cd.primaryDataFieldColumnIndex};
            else
                sliceArgs = {};
            end
            
            [dataCell, timeCell] = AnalogChannelGroupDescriptor.doComputeTransformData(cd, td, 'slice', sliceArgs, varargin{:});
        end
    end
    
     methods(Static)
        function cd = buildVectorAnalog(name, timeField, units, timeUnits, varargin)
            p = inputParser();
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'AnalogChannelDescriptor')); % used by subclasses
            p.addParameter('scaleFromLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('scaleToLims', [], @(x) isempty(x) || isvector(x));
            p.addParameter('dataClass', 'double', @ischar);
            p.addParameter('timeClass', 'double', @ischar);
            p.addParameter('transformChannelNames', {}, @(x) iscellstr(x) || isstring(x));
            p.addParameter('transformFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('transformFnMode', '', @ischar);
            
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
            
            cd.transformChannelNames = string(p.Results.transformChannelNames);
            cd.transformFn = p.Results.transformFn;
            cd.transformFnMode = p.Results.transformFnMode;
            
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
            cd.elementTypeByField(1) = cd.NUMERIC; % otherwise we'll throw errors when validating the shared matrix data as a whole
            cd.primaryDataField = dataFieldName;
            cd.primaryDataFieldColumnIndex = dataFieldColumnIndex;
            cd = cd.initialize();
        end
        
        function cd = buildVectorAnalogFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            p = inputParser();
            p.addParameter('scaleFromLims', [], @isvector);
            p.addParameter('scaleToLims', [], @isvector);
            p.addParameter('channelDescriptor', [], @(x) isa(x, 'AnalogChannelDescriptor')); % used by subclasses
            p.parse(varargin{:});
            
            if isempty(p.Results.channelDescriptor)
                cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits);
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
        
        function cd = buildTransformAnalogChannel(name, cdList, transformFn, varargin)
            cdo = cdList{1};
            transformChannelNames = string(cellfun(@(cd) cd.name, cdList, 'UniformOutput', false));
            
            p = inputParser();
            p.addParameter('units', cdo.unitsPrimary, @ischar);
            p.addParameter('scaleFromLims', cdo.scaleFromLims, @(x) isempty(x) || isvector(x));
            p.addParameter('scaleToLims', cdo.scaleToLims, @(x) isempty(x) || isvector(x));
            p.addParameter('dataClass', cdo.originalDataClassByField{1}, @ischar);
            p.addParameter('timeClass', cdo.originalDataClassByField{2}, @ischar);
            p.addParameter('transformFnMode', '', @ischar);
            
            p.parse(varargin{:});
            
            if isempty(p.Results.transformFnMode) 
                % determine automatically whether the transform fn is
                % simple
                if nargin(transformFn) == 1
                    transformFnMode = 'simple';
                else
                    transformFnMode = 'default';
                end
            else
                transformFnMode = p.Results.transformFnMode;
            end
            
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, cdo.timeField, p.Results.units, cdo.timeUnits, ...
                'transformChannelNames', transformChannelNames, 'transformFn', transformFn, ...
                'transformFnMode', transformFnMode, rmfield(p.Results, {'units', 'transformFnMode'}));            
        end
        
        function tf = testChannelsShareTimeField(cdCell)
            assert(iscell(cdCell));
            
            timeField = cellfun(@(cd) cd.timeField, cdCell, 'UniformOutput', false);
            tf = numel(unique(timeField)) == 1;
        end
        
%         function [cdCell, data] = sharedMatrixCheckConvertDataAndUpdateClass(cd, data)
%             % equivalent of
%             % checkConvertDataAndUpdateMemoryClassToMakeCompatible(cd,
%             % fieldIdx, data) but operates more effectively on groups of
%             % shared data
%         end

        function cls = getSubChannelClass()
            cls = '';
        end
    end

end
