classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    % Use the factory builder methods in subclasses rather than constructing
    % directly
    
    properties
        description = ''; % extended description
        
        meta % anything you'd like
        
        displayGroup string = "" % solely for display purposes or logical grouping, not actual data field grouping
    end
    
    properties(SetAccess={?TrialDataInterface,?TrialData,?ChannelDescriptor})
        name char = ''; % short name, must be valid field name
    end
    
    properties(SetAccess={?TrialDataInterface,?TrialData,?ChannelDescriptor})
        groupName = ''; % name of group to which this channel belongs
        
        special = false; % whether this channel is a "special" identifier channel used by TrialData
        
        required = true; % whether data for this channel is required or not
    end
    
    % set by factory builder methods or inferAttributesFromData
    % EACH OF THE FOLLOWING PROPERTIES MUST HAVE THE SAME SIZE (nFields x 1)
    properties(SetAccess={?ChannelDescriptor,?ChannelImpl})
        % these are fixed strings like "data" or "time" that describes what each field contains
        fieldIds = {};
        
        % one of the element type Constants defined below
        elementTypeByField = [];
        
        % name of each field in the data struct
        dataFields = {}
        
        % string describing units of each field
        unitsByField = {}
        
        % original class() of each field, used to maintain the same type in memory .memoryClassByFIeld
        originalDataClassByField = {};
        
        catAlongFirstDimByField = []; % manual override specifying that values may be concatenated along dim 1, this will set collectAsCell to be false
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
        CATEGORICAL = 8;
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
        isVectorizableByField % elements may be concatenated into vector
        isCategoricalByField
        
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
        
        % return a subclass of ChannelImpl
        cimpl = getImpl(cd);
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
            
            if isempty(cd.catAlongFirstDimByField)
                cd.catAlongFirstDimByField = false(cd.nFields, 1);
            end
        end
        
        function cd = validate(cd)
            cd.warnIfNoArgOut(nargout);
            
            % check that all field values are okay
            nFields = numel(cd.dataFields); %#ok<*PROP>
            
            assert(numel(cd.fieldIds) == nFields, 'fieldIds has wrong size');
            assert(numel(cd.elementTypeByField) == nFields, 'elementTypeByField has wrong size');
            assert(numel(cd.unitsByField) == nFields, 'unitsByField has wrong size');
            assert(numel(cd.originalDataClassByField) == nFields, 'originalDataClassByField has wrong size');
            assert(numel(cd.catAlongFirstDimByField) == nFields, 'catAlongFirstDimByField has wrong size');
        end
        
        function [tf, idx] = hasFieldId(cd, fieldId)
            fieldId = string(fieldId);
            if isnumeric(fieldId)
                tf = fieldId >= 1 & fieldId <= cd.nFields;
                idx = fieldId;
            else
                [tf, idx] = ismember(fieldId, cd.fieldIds);
            end
        end
        
        function idx = lookupFieldId(cd, fieldId)
            if isnumeric(fieldId)
                idx = fieldId;
            else
                [tf, idx] = ismember(fieldId, cd.fieldIds);
                if ~all(tf)
                    error('Field with id %s not found for channel %s', fieldId{find(~tf, 1)}, cd.name);
                end
            end
        end
        
        function info = getFieldInfoById(cd, fieldId)
            ind = cd.getFieldIndexById(fieldId);
            info.id = fieldId;
            info.dataField = cd.dataFields{ind};
            info.elementType = cd.elementTypesByField{ind};
            info.units = cd.unitsByField{ind};
            info.originalDataClass = cd.originalDataClassByField{ind};
            info.catAlongFirstDim = cd.catAlongFirstDimByField(ind);
            info.accessClass = cd.accessClassByField{ind};
            info.memoryClass = cd.memoryClassByField{ind};
            info.storageClass = cd.storageClassByField{ind};
            info.isShareable = cd.isShareableByField(ind);
            
        end
        
        % used by trial data when it needs to change field names
        function name = suggestFieldName(cd, fieldIdx)
            fieldIdx = cd.lookupFieldId(fieldIdx);
            if fieldIdx == 1
                name = cd.name;
            else
                name = sprintf('%s_f%s', cd.fieldIds{fieldIdx});
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
                    case cd.CATEGORICAL
                        c{iF} = 'categorical';
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
                    case {cd.SCALAR, cd.VECTOR, cd.NUMERIC, cd.DATENUM, cd.CELL}
                        if ~isempty(cd.originalDataClassByField{iF})
                            c{iF} = cd.originalDataClassByField{iF};
                        else
                            c{iF} = 'double';
                        end
                    case cd.STRING
                        c{iF} = 'string';
                    case cd.CATEGORICAL
                        c{iF} = 'categorical';
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
                    case cd.CATEGORICAL
                        c{iF} = 'categorical';
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
            iF = cd.lookupFieldId(iF);
            assert(iF ~= 1, 'Should not rename primary field using this method');
            cd.dataFields{iF} = newName;
            cd = cd.initialize();
        end
        
        function cd = setPrimaryUnits(cd, unitName)
            assert(ischar(unitName));
            cd.unitsByField{1} = unitName;
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
        
        function cd = buildSubChannelDescriptor(~, nameOrIdx) %#ok<STOUT,INUSD>
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
        function groupName = getGroupName(cd, groupName) %#ok<INUSL>
            
        end
        
        function n = get.nFields(cd)
            n = numel(cd.dataFields);
        end
        
        function tf = getIsShareableByField(cd)
            tf = true(1, cd.nFields);
            tf(1) = false;
        end
        
        function vals = get.missingValueByField(cd)
            vals = cd.getMissingValueByField();
        end
        
        function vals = getMissingValueByField(cd)
            missingVals = {false, NaN, [], [], string(missing), NaN, [], categorical(missing)};
            vals = missingVals(cd.elementTypeByField);
            accClasses = string(cd.accessClassByField); % was memory class, changed to access class so that int types converted to single get filled as NaN
            for iF = 1:numel(vals)
                if ismember(cd.elementTypeByField, [cd.VECTOR, cd.NUMERIC])
                    if ~cd.collectAsCellByField(iF)
                        % replace [] with NaN since this will be collected as
                        % matrix
                        vals{iF} = nan(1, accClasses{iF});
                    else
                        vals{iF} = nan(0, 1, accClasses{iF});
                    end
                    
                elseif cd.elementTypeByField(iF) == cd.SCALAR
                    if strcmp(accClasses{iF}, 'logical')
                        vals{iF} = false;
                    elseif ismember(accClasses(iF), ["string", "char"])
                        vals{iF} = string(missing);
                    else
                        vals{iF} = ChannelImpl.icast(vals{iF}, accClasses{iF});
                    end
                    
                elseif cd.elementTypeByField(iF) == cd.CELL
                    vals{iF} = {};
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
            tf = ~cd.isVectorizableByField;
            if ~isempty(cd.catAlongFirstDimByField)
                tf = tf & makerow(~cd.catAlongFirstDimByField);
            end
        end
        
        function tf = get.isStringByField(cd)
            tf = cd.elementTypeByField == cd.STRING;
        end
        
        function tf = get.isCategoricalByField(cd)
            tf = cd.elementTypeByField == cd.CATEGORICAL;
        end
        
        function tf = get.isVectorizableByField(cd) % returns true for boolean and scalar including string
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM, cd.STRING, cd.CATEGORICAL]);
        end
        
        function tf = get.isNumericScalarByField(cd)
            tf = ismember(cd.elementTypeByField, [cd.BOOLEAN, cd.SCALAR, cd.DATENUM]) & ...
                ismember(cd.accessClassByField, {'logical', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64', 'double', 'single'});
        end
        
        function tf = get.isVectorByField(cd)
            tf = ismember(cd.elementTypeByField, cd.VECTOR);
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
        function cd = loadobj(s)
            cd = builtin('loadobj', s);
            cd = cd.initialize();
            cd = cd.validate();
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
