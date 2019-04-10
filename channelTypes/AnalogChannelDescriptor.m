classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
        
        isColumnOfSharedMatrix = false; % this field shares a data field with other channels (whose name is stored in primaryDataFieldManual
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
        
        % if specified, this channel will be uniform (isUniform == true) and 
        % must be constrained to have / will be trusted to have a fixed sampling rate
        % this requires some upfront processing to ensure the samples are even
        % bue saves time during access and resampling
        timeDelta = []
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
        
        isUniform % has manually specified timeDelta
    end
    
    properties(Hidden)
        primaryDataFieldManual 
    end
    
    methods
        function f = get.primaryDataField(cd)
            f = cd.getPrimaryDataField();
        end
        
        function f = getPrimaryDataField(cd) % so that the get.primaryDataField can be overridden if necessary
            if cd.isColumnOfSharedMatrix
                f = cd.primaryDataFieldManual;
            elseif cd.isTransform
                f = cd.transformChannelNames{1};
            else
                f = cd.name;
            end
        end
        
        function tf = get.isTransform(cd)
            tf = ~isempty(cd.transformChannelNames);
        end
        
        function tf = get.isUniform(cd)
            tf = ~isempty(cd.timeDelta);
        end
        
        function cd = updateGroup(cd, newGroupDescriptor)
            cd.warnIfNoArgOut(nargout);
            assert(isa(newGroupDescriptor, 'AnalogChannelGroupDescriptor'));
            assert(cd.isColumnOfSharedMatrix)
            cd.primaryDataFieldManual = newGroupDescriptor.name;
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
            cd = cd.initialize();
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
            cd = cd.initialize();
        end
    end
        
    methods
        function cd = initialize(cd)
            cd.warnIfNoArgOut(nargout);
            
            cd.dataFields = {cd.primaryDataField, cd.timeField};
            cd.fieldIds = {'data', 'time'};
            
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = AnalogChannelImpl(cd);
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
            dataClass = ChannelDescriptor.getCellElementClass(dataCell); %#ok<*PROPLC>
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
