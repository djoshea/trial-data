function str = structToString(s, varargin)
    p = inputParser();
    p.addOptional('separator', ' ', @ischar);
    
    p.addParameter('includeFieldNames', true, @islogical);
    p.addParameter('fieldValueSeparator', "=", @isstringlike); % if includeFieldNames is true, e.g. Field=Value
    
    p.addParameter('fieldNameSubstitutions', struct(), @isstruct); % for includeFieldNames = true or useFieldNameForBoolean, if a field is set, use that field's value as the string instead of the field name itself
    p.addParameter('useFieldNameForBoolean', true, @islogical); % if the value is true/false, then use the attribute name or ''
    p.addParameter('suffixByField', struct(), @isstruct);
    p.addParameter('logicalNotPrefix', 'Not ', @ischar);
    p.addParameter('removeIsForLogical', true, @islogical); % "Not Is Something" --> "Not Something"
    p.parse(varargin{:});
    
    separator = p.Results.separator;
    includeFieldNames = p.Results.includeFieldNames;
    suffixByField = p.Results.suffixByField;
    fieldNameSubstitutions = p.Results.fieldNameSubstitutions;
    fieldValueSeparator = string(p.Results.fieldValueSeparator);
    
    fields = fieldnames(s);
    if isempty(fields)
        str = "";
        return;
    end
    vals = cellfun(@(fld) convertToString(s.(fld), fld), fieldnames(s), 'UniformOutput', true);

    str = string(TrialDataUtilities.String.strjoin(vals, separator));
    
    function str = convertToString(v, fld)
        if isfield(suffixByField, fld) && strlength(suffixByField.(fld)) > 0
            suffix = " " + string(suffixByField.(fld));
        else
            suffix = "";
        end
        
        if isfield(fieldNameSubstitutions, fld) && ~isempty(isfield(fieldNameSubstitutions, fld))
            fldSub = string(fieldNameSubstitutions.(fld));
        else
            fldSub = string(fld);
        end
        
        if includeFieldNames
            prefix = fldSub + fieldValueSeparator;
        else
            prefix = "";
        end
        
        if ischar(v) || (isstring(v) && numel(v) == 1)
            str = prefix + string(v) + suffix;
        elseif islogical(v) && p.Results.useFieldNameForBoolean
            if p.Results.removeIsForLogical
                if strncmp(fldSub, 'Is ', 3)
                    fldSub = extractAfter(fldSub, "Is ");
                end
            end
            if v
                str = fldSub + suffix;
            else
                str = string(p.Results.logicalNotPrefix) + fldSub + suffix;
            end
        elseif isempty(v)
            str = prefix + "[]";
        elseif isnumeric(v) || islogical(v)
            if int32(v) == v
                vstr = mat2str(v);
            else
                vstr = mat2str(v, 3);
            end
            str = prefix + string(vstr) + suffix;
        elseif iscategorical(v)
            if isscalar(v)
                str = prefix + string(v);
            else
                str = prefix + "[" + strjoin(string(v), ", ") + "]";
            end
        elseif iscellstr(v) || isstring(v) % scalar string handled above
            str = prefix + "{" + strjoin(string(v), suffix + ",") + "}";
        else
            error('Could not convert struct field value');
        end
    end        

end