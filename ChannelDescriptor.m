classdef ChannelDescriptor < matlab.mixin.Heterogeneous
% Use the factory builder methods in subclasses rather than constructing
% directly, or use inferAttributesFromData
    
    properties
        groupName = ''; % name of group to which this channel belongs 

        description = ''; % extended description 

        meta % anything you'd like

        special = false; % whether this channel is a "special" identifier channel used by TrialData
    end
    
    % set by factory builder methods or inferAttributesFromData
    % EACH OF THE FOLLOWING PROPERTIES MUST HAVE THE SAME SIZE (nFields x 1)
    properties(SetAccess=protected)
        name = ''; % short name, must be valid field name

        % one of the element type Constants defined below
        elementTypeByField = [];

        % name of each field in the data struct
        dataFields = {}
        
        % string describing units of each field
        unitsByField = {}
        
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
    end
    
    properties(Dependent)
        nFields
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
        
         % in memory data class of each field (or empty if unknown / mixed)
        dataClassByField % cell array specifying data class to convert each 
        
        % persistent storage data class of each field (or empty if unknown / mixed)
        storageDataClassByField
        
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
        
        % build a channel descriptor from the data
        cd = inferAttributesFromData(cd, dataCell);
    end
    
    methods % Dependent property definitions
        function n = get.nFields(cd)
            n = numel(cd.elementTypeByField);
        end
        
        function vals = get.missingValueByField(cd)
            missingVals = {false, NaN, [], [], '', {}};
            vals = missingVals(cd.elementTypeByField);
        end
        
        function c = get.dataClassByField(cd)
            c = cellvec(cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.DATENUM}
                        c{iF} = 'double';
                    case cd.STRING
                        c{iF} = 'char';
                    otherwise
                        c{iF} = '';
                end
            end
        end
        
        function c = get.storageDataClassByField(cd)
            c = cellvec(cd.nFields);
            for iF = 1:cd.nFields
                switch cd.elementTypeByField(iF)
                    case cd.BOOLEAN
                        c{iF} = 'logical';
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC}
                        c{iF} = cd.originalDataClassByField{iF};
                    case cd.DATENUM
                        c{iF} = 'double';
                    case cd.STRING
                        c{iF} = 'char';
                    otherwise
                        c{iF} = '';
                end
            end
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

    methods % Constructor
        function cd = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.parse(varargin{:});

            cd.name = p.Results.name;
        end

        function str = getAxisLabelPrimary(cd)
            if isempty(cd.unitsByField{1})
                str = sprintf('%s', cd.name);
            else
                str = sprintf('%s (%s)', cd.name, cd.unitsPrimary);
            end
        end
        
        function cd = renameDataField(cd, iF, newName)
            cd.dataFields{iF} = newName;
        end
        
        function cd = setPrimaryUnits(cd, unitName)
            assert(ischar(unitName));
            cd.unitsByField{1} = unitName;
        end
        
        function tf = getIsShareableByField(cd)
            tf = true(cd.nFields, 1);
            tf(1) = false;
        end

        % look at the data struct and check that the correct data fields
        % are present in the R struct
        function [ok, msg] = checkData(cd, data)
            fields = cd.dataFields;
            if ~all(isfield(data, fields))
                ok = false;
                msg = sprintf('Missing fields %s', strjoin(fields(~isfield(data, fields)), ', '));
            else
                ok = true;
                msg = '';
            end
        end

        % do any replacement of missing values, etc.
        function data = repairData(cd, data)
            % replace empty values with appropriate missing value
            for iF = 1:cd.nFields
                fld = cd.dataFields{iF};
                missingValue = cd.missingValueByField{iF};
                
                values = {data.(fld)};
                replace = cellfun(@isempty, values);
                [data(replace).(fld)] = deal(missingValue);
            end
        end

        function data = convertDataToMemoryClass(cd, data)
            % change the storage class type of the data to the in-memory
            % class
            for iF = 1:cd.nFields
                from = cd.storageDataClassByField{iF};
                to = cd.dataClassByField{iF};
                if ~isempty(from) && ~isempty(to) && ~strcmp(from, to)
                    data = structConvertFieldValues(data, to, cd.dataFields{iF});
                end
            end
        end 
        
        function data = convertDataToStorageClass(cd, data)
            for iF = 1:cd.nFields
                to = cd.storageDataClassByField{iF};
                from = cd.dataClassByField{iF};
                if ~isempty(from) && ~isempty(to) && ~strcmp(from, to)
                    data = structConvertFieldValues(data, to, cd.dataFields{iF});
                end
            end
        end 
    end
    
    methods(Static) % Utility methods
        function cls = getCellElementClass(dataCell)
            if isempty(dataCell)
                cls = 'double';
            else
                nonEmpty = cellfun(@isempty, dataCell);
                if ~any(nonEmpty)
                    cls = 'double';
                else
                    first = find(nonEmpty, 1);
                    cls = class(dataCell{first});
                end
            end
        end
    end
end
