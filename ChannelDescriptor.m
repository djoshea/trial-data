classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    % Use the factory builder methods in subclasses rather than constructing
    % directly

    properties
        description = ''; % extended description

        meta % anything you'd like
    end

    properties(SetAccess={?TrialDataInterface,?TrialData,?ChannelDescriptor})
        name = ''; % short name, must be valid field name
    end

    properties(SetAccess={?TrialDataInterface,?TrialData,?ChannelDescriptor})
        groupName = ''; % name of group to which this channel belongs

        special = false; % whether this channel is a "special" identifier channel used by TrialData

        required = true; % whether data for this channel is required or not
    end

    % set by factory builder methods or inferAttributesFromData
    % EACH OF THE FOLLOWING PROPERTIES MUST HAVE THE SAME SIZE (nFields x 1)
    properties(SetAccess=protected)
        % one of the element type Constants defined below
        elementTypeByField = [];

        % name of each field in the data struct
        dataFields = {}

        % string describing units of each field
        unitsByField = {}

        % original class() of each field, used to maintain the same type in memory .memoryClassByFIeld
        originalDataClassByField = {};

        catAlongFirstDimByField = {}; % manual override specifying that values may be concatenated along dim 1, this will set collectAsCell to be false
    end

    properties(Constant, Hidden) % element type constants
        UNKNOWN = 0;
        BOOLEAN = 1; % scalar logical
        SCALAR = 2; % scalar non-logical
        VECTOR = 3;
        NUMERIC = 4;
        STRING = 5;
        DATENUM = 6;
        CELL = 7;
    end

    properties(Dependent)
        nFields

        % nFields x 1 arrays or cells
        collectAsCellByField
        missingValueByField
        isNumericScalarByField % true for logical as well
        isBooleanByField
        isStringByField
        isVectorByField
        isScalarByField
        isNumericScalarByField

        % indicates whether a given field is shareable between multiple
        % channels. If a field is marked as shareable, it will be copied
        % when this channel's data is updated, so that the other channel's
        % are not affected. See getIsShareableByField for implementation
        isShareableByField

        % data class of each field as accessed
        accessClassByField % cell array specifying data class to convert each

        % data class of each field as stored in .data
        memoryClassByField

        % persistent storage data class of each field (or empty if unknown / mixed)
        storageClassByField

        % these are shortcuts to the first element of the properties above
        % since the first data field is considered the primary data field
        unitsPrimary
        dataFieldPrimary
    end

    methods(Abstract)
        % return a string with this channels short type description string
        type = getType(cdesc);

        % return a string with a short description of the channel excluding
        % channel name
        str = describe(cdesc);

        % update internal properties so as to handle specific data types in
        % vararargin{1}, varargin{2}, ... varargin{nDataFields}
        cd = inferAttributesFromData(cd, varargin);
    end

    methods(Access=protected) % Constructor
        function cd = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});

            cd.name = p.Results.name;
        end
    end

    methods
        function cd = initialize(cd)
            % gives the class a chance to reinitialize itself after the
            % user makes changes
            cd.warnIfNoArgOut(nargout);
        end

        % used by trial data when it needs to change field names
        function name = suggestFieldName(cd, fieldIdx)
            if fieldIdx == 1
                name = cd.name;
            else
                name = sprintf('%s_f%d', fieldIdx);
            end
        end

        function data = convertDataToCorrectVectorFormat(cd, iF, data)
            % if collectAsCellByField(iF) is true, convert to cell
            % if false, convert to numeric vector with appropriate missing
            % values inserted to equalize length

            if cd.collectAsCellByField(iF)
                if ~iscell(data)
                    data = num2cell(data);
                end
            else
                % convert to numeric vector, check for non-scalar values,
                % and fill missing values
                if iscell(data)
                    if ismember(cd.elementTypeByField(iF), [cd.VECTOR, cd.NUMERIC])
                        % numeric data cat along first dim
                        nonEmpty = ~cellfun(@isempty, data);
                        if cd.elementTypeByField(iF) == cd.VECTOR
                            data = cellfun(@makerow, data, 'UniformOutput', false);
                        end
                        try
                            mat = cell2mat(data(nonEmpty));
                        catch
                            error('Numeric data that is to be concatenated along first dim has uneven sizing');
                        end

                        data = TensorUtils.inflateMaskedTensor(mat, 1, nonEmpty, cd.missingValueByField{iF});
                    elseif ismember(cd.elementTypeByField(iF), cd.STRING)
                        if iscell(data)
                            % TODO address this
                            emptyMask = cellfun(@(x) isempty(x) || (isscalar(x) && ismissing(x)), data);
                            data(emptyMask) = {''};
                        end
                        assert(iscellstr(data) || isstring(data), 'String data field must be cellstr or string');
                        data = string(data);
                        assert(isempty(data) || isvector(data));

                        % replace empty values with missing?
                        emptyMask = arrayfun(@(x) x == "", data);
                        data(emptyMask) = string(missing());
                    else
                        missingVal = cd.missingValueByField{iF};
                        assert(isscalar(missingVal));
                        nVals = cellfun(@numel, data);
                        if any(nVals > 1)
                            throwError('Data must contain scalar values for each trial');
                        end
                        [data{nVals==0}] = deal(missingVal);
                        nel = numel(data);
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, class(missingVal));
                        data = ChannelDescriptor.cellCast(data, newClass);
                        data = cat(1, data{:});
                        assert(isempty(data) || isvector(data) && numel(data) == nel);
                    end
                end
            end

            function throwError(varargin)
                error(['Error in channel %s, field %d: ' varargin{1}], cd.name, iF, varargin{2:end});
            end
        end

        function [cd, data] = checkConvertDataAndUpdateMemoryClassToMakeCompatible(cd, fieldIdx, data)
            % this function does a couple of things. It takes a cell or
            % vector of data destined for a specific field of this channel
            % descriptor. It first checks whether this data is at all
            % acceptable for this field type (i.e. BOOLEAN, SCALAR, VECTOR,
            % above). It will throw an error if not.
            %
            % If the data are compatible, it will then check the class of
            % data against my .memoryClassByField{fieldIdx}. If it is
            % possible to convert data to memClass, data will be converted.
            % If not, then memClass will be updated to reflect the change.
            % For example, if memClass is single and data(i) = uint16(1),
            % data(i) will be converted to single(1). However, if memClass
            % is uint16 and data(i) = double(1.5), then memClass will be changed to
            % double. This will not change any other data set on this field
            % in the TrialData instance, but this preexisting data will be
            % cast into the access class on access anyway, so it's not
            % necessary to worry about it now.if not we
            % change the class to match the new format.

            % if meant to be collected as a vector, do that first as it
            % simplifies checking the types
            iF = fieldIdx;

            % convert to cell or vector depending on collectAsCellByField
            data = cd.convertDataToCorrectVectorFormat(iF, data);

            memClass = cd.memoryClassByField{iF};
            switch cd.elementTypeByField(iF)
                case cd.BOOLEAN
                    data(isnan(data)) = false;
                    convertedData = logical(data);
                    if any(convertedData ~= data)
                        throwError('Data must be logical or convertible to logical vector');
                    end
                    newClass = memClass;

                case {cd.SCALAR, cd.DATENUM}
                    newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
                    data = ChannelDescriptor.icast(data, newClass);

                case cd.VECTOR
                    if cd.collectAsCellByField(iF) && iscell(data)
                        okay = cellfun(@(x) isempty(x) || isvector(x), data);
                        if ~all(okay)
                            throwError('Data cell contents must be vectors or empty');
                        end
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);
                        data = ChannelDescriptor.cellCast(data, newClass);
                    else
                        newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
%                         if strcmp(newClass, 'logical')
%                             data(isnan(data)) = false;
%                         end
                        data = ChannelDescriptor.icast(data, newClass);
                    end

                    if cd.collectAsCellByField(iF) && iscell(data)
                        % we'll be cat'ing them along dim 1.
                        data = cellfun(@makecol, data, 'UniformOutput', false);
                    end

                case cd.NUMERIC
                    if cd.collectAsCellByField(iF) && iscell(data)
                        okay = cellfun(@(x) isempty(x) || (isnumeric(x) || islogical(x)), data);
                        if ~all(okay)
                            throwError('Data cell contents must be numeric');
                        end
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);
                        data = ChannelDescriptor.cellCast(data, newClass);
                    else
                        newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
%                         if strcmp(newClass, 'logical')
%                             data(isnan(data)) = false;
%                         end
                        data = ChannelDescriptor.icast(data, newClass);
                    end

                case cd.STRING
                    okay = cellfun(@(x) isempty(x) || (ischar(x) && isvector(x)), data);
                    if ~all(okay)
                        throwError('Data cell contents must be strings');
                    end
                    data = cellfun(@makerow, data, 'UniformOutput', false);
                    newClass = 'char';

                case cd.CELL
                    okay = cellfun(@(x) isempty(x) || iscell(x), data);
                    if ~all(okay)
                        throwError('Data cell contents must be cells');
                    end
                    % this used to be 'cell', but it wasn't correct for spike
                    % array waveforms. might need to be fixed if this isn't
                    % sufficient
                    newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);

                otherwise
                    throwError('Unknown element type')
            end

            cd.originalDataClassByField{iF} = newClass;

            function throwError(varargin)
                error(['Error in channel %s, field %d: ' varargin{1}], cd.name, fieldIdx, varargin{2:end});
            end
        end

        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = ChannelDescriptor.icast(data, accClass);
            end
        end

        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % when data is accessed by TrialData, this is one additional
            % chance to cast or adjust the data stored in .data. This
            % default implementation casts from the memory class to the
            % access class on demand
            %memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};

            data = ChannelDescriptor.cellCast(data, accClass);
            % data = cellfun(@(a) cast(a, accClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            if ~cd.collectAsCellByField(fieldIdx)
                % for vector types, we can makerow the contents to ensure
                % that cell2mat works as intended
                if cd.isVectorByField(fieldIdx)
                    data = cellfun(@makerow, data, 'UniformOutput', false);
                end

                % because of the way missingValue works, empty trials will
                % be NaN for invalid trials when setting channel data
                % trials missing values will have a single Nan
                nanMask = cellfun(@(x) isscalar(x) && ismissing(x), data);
                if any(~nanMask)
                    data = cat(1, data{~nanMask});
                else
                    data = data{1}([], :, :, :, :, :, :);
%                     data = zeros(0, 1);
                end
                data = TensorUtils.inflateMaskedTensor(data, 1, ~nanMask);
            end
        end

        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = ChannelDescriptor.icast(data, memClass);
            end
        end

        function data = convertAccessDataCellToMemory(cd, fieldIdx, data)
            % when data is accessed by TrialData, this is one additional
            % chance to cast or adjust the data stored in .data. This
            % default implementation casts from the access class to the
            % memory class on demand
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = cellfun(@(a) ChannelDescriptor.icast(a, memClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            end
        end

        function c = getAccessClassByField(cd)
            % implements get.accessClassByField
            memClass = cd.getMemoryClassByField();
            c = cell(1, cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.DATENUM, cd.CELL}
                        if ismember(memClass{iF}, {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single'})
                            c{iF} = 'single';
                        elseif ismember(memClass{iF}, {'logical'})
                            c{iF} = 'logical';
                        elseif ismember(memClass{iF}, {'categorical'})
                            c{iF} = 'categorical';
                        elseif ismember(memClass{iF}, {'string'})
                            c{iF} = 'string';
                        else
                            c{iF} = 'double';
                        end
                    case cd.STRING
                        c{iF} = 'string';
                    otherwise
                        c{iF} = '';
                end
            end
        end

        function c = getMemoryClassByField(cd)
            % implements get.memoryClassByField, differs in that numeric
            % types default to their original data class
            c = cell(1, cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.DATENUM}
                        if ~isempty(cd.originalDataClassByField{iF})
                            c{iF} = cd.originalDataClassByField{iF};
                        else
                            c{iF} = 'double';
                        end
                    case cd.STRING
                        c{iF} = 'string';
                    otherwise
                        c{iF} = cd.originalDataClassByField{iF};
                end
            end
        end

        function c = getStorageClassByField(cd)
            c = cell(1, cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.CELL}
                        if ~isempty(cd.originalDataClassByField{iF})
                            c{iF} = cd.originalDataClassByField{iF};
                        else
                            c{iF} = 'double';
                        end
                    case cd.DATENUM
                        c{iF} = 'double';
                    case cd.STRING
                        c{iF} = 'string';
%                     case cd.CELL
%                         c{iF} = 'cell';
                    otherwise
                        c{iF} = '';
                end
            end
        end

        function str = getAxisLabelPrimary(cd)
            if isempty(cd.unitsByField{1})
                str = sprintf('%s', cd.name);
            else
                str = sprintf('%s (%s)', cd.name, cd.unitsPrimary);
            end
        end

        function [cd, dataFieldRenameStruct] = rename(cd, newName)
            % renames the channel and any fields within
            % dataFieldRenameStruct indicates which .data fields will need
            % to be renamed, as dataFieldRenameStruct.oldName = newName
            cd.warnIfNoArgOut(nargout);
            oldName = cd.name;
            cd.name = newName;

            dataFieldRenameStruct = struct();
            if strcmp(cd.dataFields{1}, oldName)
                dataFieldRenameStruct.(oldName) = newName;
                cd.dataFields{1} = newName;
            end
            cd = cd.initialize();
        end

        function cd = renameDataField(cd, iF, newName)
            cd.warnIfNoArgOut(nargout);
            assert(iF ~= 1, 'Should not rename primary field using this method');
            cd.dataFields{iF} = newName;
            cd = cd.initialize();
        end

        function cd = setPrimaryUnits(cd, unitName)
            assert(ischar(unitName));
            cd.unitsByField{1} = unitName;
        end

        % look at the data struct and check that the correct data fields
        % are present in the R struct
        function [ok, missing] = checkData(cd, data)
            fields = cd.dataFields;
            if ~all(isfield(data, fields))
                ok = false;
                missing = fields(~isfield(data, fields));
            else
                ok = true;
                missing = {};
            end
        end

        function data = addMissingFields(cd, data)
            % manually add the fields needed by this channel to data
            % struct, filling them with the appropriate missing value
            for iF = 1:cd.nFields
                fld = cd.dataFields{iF};
                if ~isfield(data, fld)
                    % insert field with appropriate missing values
                    [data(:).(fld)] = deal(cd.missingValueByField{iF});
                end
            end
        end

        function data = correctMissingValueInData(cd, iF, data)
            missingValue = cd.missingValueByField{iF};

            if iscell(data)
                replace = cellfun(@(x) isempty(x) || (isscalar(x) && ismissing(x)), values);
                [data(replace)] = deal(missingValue);
            else
                replace = isnan(values);
                data(replace) = missingValue;
            end
        end

        function [tf, idx] = hasSubChannel(~, name)
            % overwritten by subclasses
            if ischar(name)
                tf = false;
                idx = NaN;
            else
                tf = false(numel(name), 1);
                idx = nan(numel(name), 1);
            end
        end

        function cd = buildSubChannelDescriptor(~, nameOrIdx) %#ok<INUSD>
            error('Not implemented');
        end

        function [names, chidx] = listNamedSubChannels(~)
            names = string([]);
            chidx = [];
        end
    end

    methods(Sealed)
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect', ...
                    class(obj));
            end
        end
    end

    % Dependent property definitions, some refer out to getProperty methods
    % of the same name to give subclasses an opportunity to redefine their
    % implementation
    methods
        function grp = get.groupName(cd)
            grp = cd.getGroupName(cd.groupName);
        end

        % allows for sub classes to overwrite
        function groupName = getGroupName(cd, groupName)

        end
        function n = get.nFields(cd)
            n = numel(cd.dataFields);
        end

        function tf = getIsShareableByField(cd)
            tf = true(1, cd.nFields);
            tf(1) = false;
        end

        function vals = get.missingValueByField(cd)
            missingVals = {false, NaN, [], [], string(missing), NaN, []};
            vals = missingVals(cd.elementTypeByField);
            accClasses = string(cd.accessClassByField); % was memory class, changed to access class so that int types converted to single get filled as NaN
            for iF = 1:numel(vals)
                if ~cd.collectAsCellByField(iF) && ismember(cd.elementTypeByField, [cd.VECTOR, cd.NUMERIC])
                    % replace [] with NaN since this will be collected as
                    % matrix
                    vals{iF} = nan(1, accClasses{iF});

                elseif cd.elementTypeByField(iF) == cd.SCALAR
                    if strcmp(accClasses{iF}, 'logical')
                        vals{iF} = false;
                    elseif ismember(accClasses(iF), ["string", "char"])
                        vals{iF} = string(missing);
                    else
                        vals{iF} = cast(vals{iF}, accClasses{iF});
                    end
                    vals{iF} = ChannelDescriptor.icast(vals{iF}, accClasses{iF});

                end
            end
        end

        function c = get.accessClassByField(cd)
            % defer to method to make it overrideable in subclasses
            c = cd.getAccessClassByField();
        end

        function c = get.memoryClassByField(cd)
            % defer to method to make it overrideable in subclasses
            c = cd.getMemoryClassByField();
        end

        function c = get.storageClassByField(cd)
            c = cd.getStorageClassByField();
        end

        function tf = get.collectAsCellByField(cd)
            tf = ~cd.isScalarByField;
            if ~isempty(cd.catAlongFirstDimByField)
                tf = tf & makecol(~cd.catAlongFirstDimByField);
            end
        end

        function tf = get.isStringByField(cd)
            tf = cd.elementTypeByField == cd.STRING;
        end

        function tf = get.isScalarByField(cd) % returns true for boolean and scalar including string
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM, cd.STRING]);
        end

        function tf = get.isNumericScalarByField(cd)
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM]) & ...
                ismember(cd.accessClassByField, {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64', 'double', 'single'});
        end

        function tf = get.isVectorByField(cd)
            tf = ismember(cd.elementTypeByField, cd.VECTOR);
        end

        function tf = get.isBooleanByField(cd)
            tf = cd.elementTypeByField == cd.BOOLEAN;
        end

        function tf = get.isNumericScalarByField(cd)
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM]);
        end

        function tf = get.isShareableByField(cd)
            % defer to overrideable function
            tf = cd.getIsShareableByField();
        end

        function u = get.unitsPrimary(cd)
            if isempty(cd.unitsByField)
                u = '';
            else
                u = cd.unitsByField{1};
            end
        end

        function f = get.dataFieldPrimary(cd)
            if isempty(cd.dataFields)
                f = '';
            else
                f = cd.dataFields{1};
            end
        end
    end

    methods(Static) % Utility methods
        function data = icast(data, newClass)
            if ~isa(data, newClass)
                if strcmp(newClass, 'logical')
                    data(isnan(data)) = false;
                    data = logical(data);
                elseif strcmp(newClass, 'categorical')
                    data = categorical(data);
                else
                    data = cast(data, newClass);
                end
            end
        end

        function data = scaleData(data, scaleFromLims, scaleToLims)
            if isempty(scaleFromLims) || isempty(scaleToLims)
                return;
            end
            scaleFromLow = scaleFromLims(2);
            scaleFromRange = scaleFromLims(2) - scaleFromLims(1);
            scaleToLow = scaleToLims(2);
            scaleToRange = scaleToLims(2) - scaleToLims(1);
            if iscell(data)
                data = cellfun(@(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow, data, 'UniformOutput', false);
            else
                data = (data-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
            end
        end

        function data = unscaleData(data, scaleFromLims, scaleToLims)
            data = ChannelDescriptor.scaleData(data, scaleToLims, scaleFromLims);
        end

        function cls = getCellElementClass(dataCell)
            if isempty(dataCell)
                cls = 'double';
            elseif ~iscell(dataCell)
                cls = class(dataCell);
            else
                nonEmpty = ~cellfun(@isempty, dataCell);
                if ~any(nonEmpty)
                    cls = 'double'; % assume double if no values found
                else
                    first = find(nonEmpty, 1);
                    cls = class(dataCell{first});
                end
            end
        end

        function sz = getCellElementSize(dataCell)
            if isempty(dataCell)
                sz = [];
            elseif ~iscell(dataCell)
                sz = size(dataCell);
            else
                nonEmpty = ~cellfun(@isempty, dataCell);
                if ~any(nonEmpty)
                    sz = [];
                else
                    first = find(nonEmpty, 1);
                    sz = size(dataCell{first});
                end
            end
        end

        function data = cellCast(data, newClass)
            for i = 1:numel(data)
                data{i} = ChannelDescriptor.icast(data{i}, newClass);
            end
        end
        
        function [ch, idx] = parseIndexedChannelName(name)
            % if it is a name like analogGroup(5), make a sub channel
            % channelDescriptor for analogGroup referring to column 5
            tokens = regexp(name, '(?<ch>[\w_]+)\((?<idx>\d+)\)', 'names');
            if isempty(tokens)
                ch = name;
                idx = NaN;
            else
                ch = tokens.ch;
                idx = str2double(tokens.idx);
            end
        end

        function cls = getSubChannelClass()
            cls = 'ChannelDescriptor';
        end
    end
end
