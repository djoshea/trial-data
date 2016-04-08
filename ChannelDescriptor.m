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
                    missing = cd.missingValueByField{iF};
                    assert(isscalar(missing));
                    nVals = cellfun(@numel, data);
                    if any(nVals > 1)
                        throwError('Data must contain scalar values for each trial');
                    end
                    [data{nVals==0}] = deal(missing);
                    nel = numel(data);
                    newClass = ChannelDescriptor.cellDetermineCommonClass(data);
                    data = ChannelDescriptor.cellCast(data, newClass);
                    data = cell2mat(data);
                    assert(isvector(data) && numel(data) == nel);
                end
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
            % data(i) will be converted to double(1). However, if memClass
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
                    newClass = determineCommonClass(data, memClass);
                    data = cast(data, newClass);
                    
                case cd.VECTOR
                    okay = cellfun(@(x) isempty(x) || isvector(x), data);
                    if ~all(okay)
                        throwError('Data cell contents must be vectors or empty');
                    end
                    newClass = ChannelDescriptor.cellDetermineCommonClass(data, memClass); 
                    data = ChannelDescriptor.cellCast(data, newClass);
                    data = cellfun(@makecol, data, 'UniformOutput', false);
                        
                case cd.NUMERIC
                    okay = cellfun(@(x) isempty(x) || isnumeric(x), data);
                    if ~all(okay)
                        throwError('Data cell contents must be numeric');
                    end
                    newClass = ChannelDescriptor.cellDetermineCommonClass(data, memClass); 
                    data = ChannelDescriptor.cellCast(data, newClass);
                                        
                case cd.STRING
                    okay = cellfun(@(x) isempty(x) || (ischar(x) && isvector(x)), data);
                    if ~all(okay)
                        throwError('Data cell contents must be strings');
                    end
                    data = cellfun(@makerow, data, 'UniformOutput', false);
                    newClass = 'char';
                    
                otherwise
                    throwError('Unknown element type')
            end
            
            cd.originalDataClassByField{iF} = newClass;
            
            function newClass = determineCommonClass(data, memClass)
                convertedData = cast(data, memClass);
                if ~isempty(data) && ~isequal(convertedData, data)
                    newClass = class(data);
                else
                    newClass = memClass;
                end
            end
       
            function throwError(varargin)
                error(['Error in channel %s, field %d: ' varargin{1}], cd.name, fieldIdx, varargin{2:end});
            end
           
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
            else
                % convert to cell anyway
                if iscell(data) &&  ~cd.collectAsCellByField(fieldIdx)
                    data = cell2mat(data);
                end
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
        
        function data = correctMissingValueInData(cd, fieldIdx, data)
            missingValue = cd.missingValueByField{iF};
            
            if iscell(data)
                replace = cellfun(@(x) isempty(x) || (isscalar(x) && isnan(x)), values);
                [data(replace)] = deal(missingValue);
            else
                replace = isnan(values);
                data(replace) = missingValue;
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
            missingVals = {false, NaN, [], [], '', NaN, []};
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
        function newClass = cellDetermineCommonClass(data, origClass)
            if nargin > 1
                newClass = origClass;
            else
                newClass = '';
            end
            for iV = 1:numel(data)
                if isempty(data{iV}), continue; end
                if isempty(newClass)
                    newClass = class(data{iV});
                else
                    if ~isa(data{iV}, newClass)
                        convertedData = cast(data{iV}, newClass);
                        if ~isequal(convertedData, data{iV})
                            newClass = class(data{iV});
                        end
                    end
                end
            end      
        end
        
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
        
        function data = cellCast(data, newClass)
            for i = 1:numel(data)
                if ~isa(data{i}, newClass)
                    data{i} = cast(data{i}, newClass);
                end
            end
        end
    end
end
