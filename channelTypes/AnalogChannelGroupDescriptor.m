classdef AnalogChannelGroupDescriptor < ChannelDescriptor
    properties
        scaleFromLims
        scaleToLims
        
        subChannelNames string % nChannels x 1 cellstr of sub channel names
        subChannelUnits string % nChannels x 1 cellstr of sub channel units
        
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
    end
    
    properties(Dependent)
        timeField
        timeUnits
        hasScaling
        
        nChannels
        
        % is this channel just a transformed version of a different channel?
        isTransform
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
    end
    
    methods(Access=protected)
        function cd = AnalogChannelGroupDescriptor(name, timeField, sampleSize)
            cd = cd@ChannelDescriptor(name);
            if nargin < 2 || isempty(timeField)
                timeField = sprintf('%s_time', cd.name);
            end
            
            cd.dataFields = {name, timeField};
            cd.originalDataClassByField = {'double', 'double'};
            cd.elementTypeByField = [cd.NUMERIC, cd.VECTOR];
            cd.sampleSize = sampleSize;
        end
    end
        
    methods
        function cd = initialize(cd)
            cd.warnIfNoArgOut(nargout);
            if cd.isTransform
                primaryField = cd.transformChannelNames{1};
            else
                primaryField = cd.name;
            end
                
            cd.dataFields = {primaryField, cd.timeField};
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
        
        function [dataCell, timeCell] = computeTransformDataRaw(cd, td, varargin)
            [dataCell, timeCell] = AnalogChannelGroupDescriptor.doComputeTransformData(cd, td, varargin{:});
        end
             
        function [tf, idx] = hasSubChannel(cd, name)
            [tf, idx] = ismember(string(name), cd.subChannelNames);
        end
        
        function [names, chidx] = listNamedSubChannels(cd)
            names = cd.subChannelNames;
            mask = ~cellfun(@isempty, cellstr(names));
            names = string(names);
            names = names(mask);
            chidx = find(mask);
        end
        
        function cd = setSubChannelInfo(cd, chidx, names, units)
            cd.warnIfNoArgOut(nargout);
            assert(~isnan(cd.nChannels), 'Number of channels not manually set for AnalogChannelGroupDescriptor');
            assert(all(TrialDataUtilities.Data.indexInRange(chidx, cd.nChannels)));
            cd.subChannelNames(chidx) = string(names);
            if ~isempty(units)
                cd.subChannelUnits(chidx) = string(units);
            end
        end
    end
    
    % transforms of other groups
    methods 
        function tf = get.isTransform(cd)
            tf = ~isempty(cd.transformChannelNames);
        end
    end
    
    methods(Static)
        function [dataCell, timeCell] = doComputeTransformData(cd, td, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical);
            p.addParameter('slice', {}, @(x) true); % this is used to index specifically into each sample
            p.addParameter('timeCell', [], @(x) true);
            p.addParameter('sort', false, @islogical);
            p.parse(varargin{:});
            
            % get raw data from trial data
            [dataCell, timeCell] = td.getAnalogChannelGroupMulti(cd.transformChannelNames, 'raw', true, ...
                'applyScaling', p.Results.applyScaling, 'sort', p.Results.sort);
            
            % use transform function
            switch cd.transformFnMode
                case 'simple'
                    % simple function, call it on one trial at a time and we'll
                    % handle the slicing
                    args = p.Results.slice;
                    if ~iscell(args)
                        args = {args};
                    end

                    prog = ProgressBar(numel(dataCell), 'Computing transform analog channel on the fly');
                    for iT = 1:numel(dataCell)
                        if ~isempty(dataCell{iT})
                            prog.update(iT);
                            dataCell{iT} = cd.transformFn(dataCell{iT});
                            if ~isempty(args)
                                dataCell{iT} = dataCell{iT}(:, args{:}); % take a slice through the data
                            end
                        end
                    end
                    prog.finish();
                    timeCell = timeCell(:, 1);
                    
                case 'default'
                    % more advanced function, can handle all data at once and slicing
                    if ~isempty(p.Results.slice)
                        args = p.Results.slice;
                        if ~iscell(args)
                            args = {args};
                        end
                    else
                        args = {};
                    end
                    
                   [dataCell, timeCell] = cd.transformFn(dataCell, timeCell, 'slice', args, 'scalingApplied', p.Results.applyScaling);
            end
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
