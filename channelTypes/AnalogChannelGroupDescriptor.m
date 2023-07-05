classdef AnalogChannelGroupDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
        
        subChannelNames (:, 1) string % nChannels x 1 cellstr of sub channel names
        subChannelUnits (:, 1) string % nChannels x 1 cellstr of sub channel units
        
        sampleSize % size of a single sample
        
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
        timeField
        timeUnits
        hasScaling
        
        nChannels
        
        % is this channel just a transformed version of a different channel?
        isTransform
        
        isUniform
    end

    methods
        function tf = get.hasScaling(cd)
            tf = ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims) ...
                && ~isequal(cd.scaleFromLims, cd.scaleToLims);
        end
        
        function cd = withNoScaling(cd)
            cd.warnIfNoArgOut(nargout);
            cd.originalDataClassByField{1} = 'single';
            cd.scaleFromLims = [];
            cd.scaleToLims = [];
        end
        
        function n = get.nChannels(cd)
            if isempty(cd.sampleSize)
                n = NaN;
            else
                n = prod(cd.sampleSize);
            end
        end
        
        function v = get.subChannelNames(cd)
            subChannelNames = cd.getSubChannelNames(cd.subChannelNames); %#ok<*PROP>
            if isempty(subChannelNames)
                if ~isnan(cd.nChannels)
                    v = repmat("", cd.nChannels, 1);
                else
                    v = string([]);
                end
            else
                if isnan(cd.nChannels)
                    v = string([]);
                else
                    v = subChannelNames(1:cd.nChannels);
                end
            end
        end
        
        function names = getSubChannelNames(cd, namesStored) %#ok<INUSL>
            % allows sub classes to override
            names = namesStored;
        end
        
        function vals = getMissingValueByField(cd)
            vals = getMissingValueByField@ChannelDescriptor(cd);
            accClasses = string(cd.accessClassByField);
            vals{1} = nan(0, cd.nChannels, accClasses{1});
        end
        
        function v = get.subChannelUnits(cd)
            if isempty(cd.subChannelUnits)
                if ~isnan(cd.nChannels)
                    v = repmat("", cd.nChannels, 1);
                else
                    v = string([]);
                end
            else
                if isnan(cd.nChannels)
                    v = string([]);
                else
                    v = cd.subChannelUnits(1:cd.nChannels);
                end
            end
        end
        
        function tf = get.isUniform(cd)
            tf = ~isempty(cd.timeDelta);
        end
    end
    
    methods(Access=protected)
        function cd = AnalogChannelGroupDescriptor(name, timeField, sampleSize)
            cd = cd@ChannelDescriptor(name);
            if nargin < 2 || isempty(timeField)
                timeField = sprintf('%s_time', cd.name);
            end
            
            cd.dataFields = {char(name), char(timeField)};
            cd.originalDataClassByField = {'double', 'double'};
            cd.elementTypeByField = [cd.NUMERIC, cd.VECTOR];
            cd.sampleSize = sampleSize;
        end
    end
        
    methods
        function cd = initialize(cd)
            cd.warnIfNoArgOut(nargout);
            
            if isempty(cd.dataFields)
                cd.dataFields = {cd.primaryDataField, cd.timeField};
            end
            
            if isempty(cd.fieldIds)
                cd.fieldIds = {'data', 'time'};
            end
            
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = AnalogChannelGroupImpl(cd);
        end

        % used by trial data when it needs to change field names
        function name = suggestFieldName(cd, fieldIdx)
            fieldIdx = cd.lookupFieldId(fieldIdx);
            if fieldIdx == 1
                name = cd.name;
            elseif fieldIdx == 2
                name = sprintf('%s_time', cd.name);
            else                
                name = suggestFieldName@ChannelDescriptor(cd, fieldIdx);
            end
        end
        
        function f = get.timeField(cd)
            if cd.nFields < 2
                f = '';
            else
                f = cd.dataFields{2};
            end
        end
        
        function u = get.timeUnits(cd)
            if cd.nFields < 2
                u = '';
            else
                u = cd.unitsByField{2};
            end
        end
        
        function type = getType(~)
            type = 'analogGroup';
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
            cd.sampleSize = ChannelDescriptor.getCellElementSize(dataCell);
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
        
        function cdIndividual = buildIndividualSubChannel(cd, index)
            units = cd.subChannelUnits{index};
            name = cd.subChannelNames{index};
            if isempty(name)
                name = AnalogChannelGroupDescriptor.generateDefaultSubChannelName(cd, index);
            end
            cdIndividual = AnalogChannelDescriptor.buildSharedMatrixColumnAnalog(name, cd.name, index, cd.timeField, units, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2}, ...
                'transformChannelNames', cd.transformChannelNames, 'transformFn', cd.transformFn, ...
                'transformFnMode', cd.transformFnMode);
        end
        
        function cd = buildSubChannelDescriptor(cd, nameOrIdx) 
            if isnumeric(nameOrIdx)
                index = nameOrIdx;
            else
                [tf, index] = ismember(nameOrIdx, cd.subChannelNames);
                assert(all(tf), 'Channel not found');
            end
            cd = cd.buildIndividualSubChannel(index);
        end
        
        function sz = getSampleSize(cd, td)
            if isempty(cd.sampleSize)
                fld = cd.dataFieldPrimary;
                for iT = 1:td.nTrials
                    if ~isempty(td.data(iT).(fld))
                        sz = size(td.data(iT).(fld));
                        sz = sz(2:end);
                        return
                    end
                end
            else
                sz = makerow(cd.sampleSize);
            end
        end  
        
        
        function [tf, idx] = hasSubChannel(cd, name)
            [tf, idx] = ismember(string(name), cd.subChannelNames);
        end
        
        function [names, chidx] = listNamedSubChannels(cd)
            names = cd.subChannelNames;
            mask = strlength(names) > 0;
            names = string(names);
            names = names(mask);
            chidx = find(mask);
        end
        
        function cd = setSubChannelInfo(cd, chidx, names, units)
            cd.warnIfNoArgOut(nargout);
            assert(~isnan(cd.nChannels), 'Number of channels not manually set for AnalogChannelGroupDescriptor');
            assert(all(TrialDataUtilities.Data.indexInRange(chidx, cd.nChannels)));
            cd.subChannelNames(chidx) = string(names);
            if nargin > 3 && ~isempty(units)
                cd.subChannelUnits(chidx) = string(units);
            end
        end
        
        function cd = filterSubChannels(cd, colIdx, newSampleSize)
            % generate an equivalent analog channel group with reduced size
            cd.warnIfNoArgOut(nargout);
            
            if nargin < 3 || isempty(newSampleSize)
                if islogical(colIdx)
                    newSampleSize = nnz(colIdx);
                else
                    newSampleSize = numel(colIdx);
                end
            end
            
            names = cd.subChannelNames;
            units = cd.subChannelUnits;
            cd.sampleSize = newSampleSize;
            cd.subChannelNames = names(colIdx);
            cd.subChannelUnits = units(colIdx); 
        end
    end
    
    % transforms of other groups
    methods 
        function tf = get.isTransform(cd)
            tf = ~isempty(cd.transformChannelNames);
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
            p.addParameter('sampleSize', [], @isvector);
            p.addParameter('subChannelNames', [], @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('subChannelUnits', [], @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('transformChannelNames', {}, @(x) iscellstr(x) || isstring(x));
            p.addParameter('transformFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('transformFnMode', '', @ischar);
            p.parse(varargin{:});
            
            if isempty(p.Results.channelDescriptor)
                cd = AnalogChannelGroupDescriptor(name, timeField, p.Results.sampleSize);
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
            if ~isempty(cd.transformChannelNames)
                assert(~isempty(p.Results.sampleSize), 'Sample size required for transform groups');
            end
            cd.transformFn = p.Results.transformFn;
            cd.transformFnMode = p.Results.transformFnMode;
            cd.sampleSize = p.Results.sampleSize;
            
            subChannelNames = string(p.Results.subChannelNames);
            subChannelUnits = string(p.Results.subChannelUnits);
            if isempty(p.Results.sampleSize) && ~isempty(subChannelNames)
                cd.sampleSize = numel(subChannelNames);
            else
                assert(isempty(subChannelNames) || prod(cd.sampleSize) == numel(subChannelNames), 'Sub channel names must match prod(sampleSize)');
                assert(isempty(subChannelUnits) || prod(cd.sampleSize) == numel(subChannelUnits), 'Sub channel names must match prod(sampleSize)');
            end
            % these fields should only be accessed after sample size is correctly set because they are autopopulated in get.subChannel* to be sampleSize long
            cd.subChannelNames = subChannelNames;
            cd.subChannelUnits = subChannelUnits;
            
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
        
        function cd = buildTransformAnalogGroup(name, cdList, transformFn, outputSize, varargin)
            cdo = cdList(1);
            transformChannelNames = arrayfun(@(cd) string(cd.name), cdList);
            
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
            
            cd = AnalogChannelGroupDescriptor.buildAnalogGroup(name, cdo.timeField, p.Results.units, cdo.timeUnits, ...
                'transformChannelNames', transformChannelNames, 'transformFn', transformFn, ...
                'transformFnMode', transformFnMode, ...
                'sampleSize', outputSize, rmfield(p.Results, {'units', 'transformFnMode'}));            
        end
        
        function name = generateDefaultSubChannelName(cd, index)
            nCh = cd.nChannels;
            nPad = ceil(log10(nCh));
            name = sprintf('%s_%0*d', cd.name, nPad, index);
        end
        
        function cls = getSubChannelClass()
            cls = 'AnalogChannelDescriptor';
        end
    end

end
