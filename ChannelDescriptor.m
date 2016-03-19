classdef ChannelDescriptor < matlab.mixin.Heterogeneous
% Use the factory builder methods in subclasses rather than constructing
% directly

    properties
        description = ''; % extended description

        meta % anything you'd like
    end

    properties(SetAccess={?TrialDataInterface,?TrialData})
        name = ''; % short name, must be valid field name
    end

    properties(SetAccess={?TrialDataInterface,?TrialData})
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
    end

    properties(Constant, Hidden) % element type constants
        UNKNOWN = 0;
        BOOLEAN = 1;
        SCALAR = 2;
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
        isBooleanByField
        isStringByField
        isVectorByField
        isScalarByField

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
            p.addOptional('name', '', @ischar);
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
        
        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = cast(data, accClass);
            end
        end

        function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % when data is accessed by TrialData, this is one additional
            % chance to cast or adjust the data stored in .data. This
            % default implementation casts from the memory class to the
            % access class on demand
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = cellfun(@(a) cast(a, accClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            end
        end
        
        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = cast(data, memClass);
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
                data = cellfun(@(a) cast(a, memClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            end
        end

        function c = getAccessClassByField(cd)
            % implements get.accessClassByField
            memClass = cd.getMemoryClassByField();
            c = cellvec(cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.DATENUM}
                        if ismember(memClass{iF}, {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single'})
                            c{iF} = 'single';
                        else
                            c{iF} = 'double';
                        end
                    case cd.STRING
                        c{iF} = 'char';
                    otherwise
                        c{iF} = '';
                end
            end
        end

        function c = getMemoryClassByField(cd)
            % implements get.memoryClassByField, differs in that numeric
            % types default to their original data class
            c = cellvec(cd.nFields);
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
                        c{iF} = 'char';
                    otherwise
                        c{iF} = cd.originalDataClassByField{iF};
                end
            end
        end

        function c = getStorageClassByField(cd)
            c = cellvec(cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC}
                        if ~isempty(cd.originalDataClassByField{iF})
                            c{iF} = cd.originalDataClassByField{iF};
                        else
                            c{iF} = 'double';
                        end
                    case cd.DATENUM
                        c{iF} = 'double';
                    case cd.STRING
                        c{iF} = 'char';
                    case cd.CELL
                        c{iF} = 'cell';
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

        % do any replacement of missing values, etc.
        % also adjust the data class inside channelDescriptor to match the
        % data
        function [data, cd] = repairData(cd, data)
            % replace empty values with appropriate missing value
            for iF = 1:cd.nFields
                fld = cd.dataFields{iF};
                missingValue = cd.missingValueByField{iF};

                % replace empty values, with the appropriate missing value
                % so that [] ends up as '' or NaN where appropriate
                values = {data.(fld)};
                replace = cellfun(@isempty, values);
                [data(replace).(fld)] = deal(missingValue);

                if cd.elementTypeByField(iF) == cd.BOOLEAN
                    % manually convert to logical
                    for i = 1:numel(data)
                        if ~isnan(data(i).(fld))
                            data(i).(fld) = logical(data(i).(fld));
                        else
                            data(i).(fld) = false;
                        end
                    end
                end

                if cd.elementTypeByField(iF) == cd.VECTOR
                   % manually convert to column
                    for i = 1:numel(data)
                        data(i).(fld) = makecol(data(i).(fld));
                    end
                end

                % check the data types match the specified memory data
                % types
                maskEmpty = arrayfun(@(t) isempty(t.(fld)), data);
                if all(maskEmpty), continue, end;
                % get all unique data classes ignoring empty (which are
                % typically created as [] and tend to be double)
                actualClasses = unique(arrayfun(@(t) class(t.(fld)), data(~maskEmpty), 'UniformOutput', false));
                if ~ismember(cd.originalDataClassByField{iF}, actualClasses) && numel(actualClasses) == 1
                    % all the actual data is a different class than its
                    % specified in channelDescriptor and the classes are
                    % uniform. Presumably the wrong class is specified in
                    % the channelDescriptor, so we should change it
                    cd.originalDataClassByField{iF} = actualClasses{1};
                end

                % now convert everything to the memory class, which will
                % match the originalDataClassByField for numeric types
                for i = 1:numel(data)
                    memClass = cd.memoryClassByField{iF};
                    if ~strcmp(class(data(i).(fld)), memClass) && ...
                      ismember(memClass, {'double', 'single', 'logical', 'char', 'int8', 'uint8', 'uint16', 'int16', 'uint32', 'int32', 'unit64', 'int64'})
                        data(i).(fld) = cast(data(i).(fld), memClass);
                    end
                end
            end
        end

        function data = convertDataToMemoryClass(cd, data)
            % change the storage class type of the data to the in-memory
            % class
            for iF = 1:cd.nFields
                from = cd.storageClassByField{iF};
                to = cd.memoryClassByField{iF};
                if ~isempty(from) && ~isempty(to) && ~strcmp(from, to)
                    data = structConvertFieldValues(data, to, cd.dataFields{iF});
                end
            end
        end

        function data = convertDataToStorageClass(cd, data)
            for iF = 1:cd.nFields
                to = cd.storageClassByField{iF};
                from = cd.memoryClassByField{iF};
                if ~isempty(from) && ~isempty(to) && ~strcmp(from, to)
                    data = structConvertFieldValues(data, to, cd.dataFields{iF});
                end
            end
        end

        function vec = vectorWithMissingValue(cd, n, fieldnum)
            if nargin < 3
                fieldnum = 1;
            end
            missing = cd.missingValueByField{fieldnum};
            if cd.collectAsCellByField(fieldnum)
                vec = repmat({missing}, n, 1);
            else
                vec = repmat(missing, n, 1);
            end
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
        function n = get.nFields(cd)
            n = numel(cd.dataFields);
        end

        function tf = getIsShareableByField(cd)
            tf = true(cd.nFields, 1);
            tf(1) = false;
        end

        function vals = get.missingValueByField(cd)
            missingVals = {false, NaN, [], [], '', {}, []};
            vals = missingVals(cd.elementTypeByField);
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
        end

        function tf = get.isStringByField(cd)
            tf = cd.elementTypeByField == cd.STRING;
        end

        function tf = get.isScalarByField(cd) % returns true for boolean and scalar
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM]);
        end

        function tf = get.isBooleanByField(cd)
            tf = cd.elementTypeByField == cd.BOOLEAN;
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

        function cls = getCellElementClass(dataCell)
            if isempty(dataCell)
                cls = 'double';
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
    end
end
