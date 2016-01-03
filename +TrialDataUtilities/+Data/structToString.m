function str = structToString(s, varargin)
    p = inputParser();
    p.addOptional('separator', ' ', @ischar);
    p.addParameter('includeFieldNames', true, @islogical);
    p.addParameter('useFieldNameForBoolean', true, @islogical); % if the value is true/false, then use the attribute name or ''
    p.parse(varargin{:});
    
    separator = p.Results.separator;
    includeFieldNames = p.Results.includeFieldNames;

    fields = fieldnames(s);
    if isempty(fields)
        str = '';
        return;
    end
    vals = cellfun(@(fld) convertToString(s.(fld), fld), fieldnames(s), 'UniformOutput', false);

    if includeFieldNames
        str = strjoin(cellfun(@(fld, val) strcat(fld, '=', val), fields, vals, 'UniformOutput', false), separator);
    else
        str = strjoin(vals, separator);
    end
    
    return;
    
    function str = convertToString(v, fld)
        if ischar(v)
            str = v;
        elseif islogical(v) && p.Results.useFieldNameForBoolean && ~includeFieldNames
            if v
                str = fld;
            else
                str = strcat('~', fld);
            end
        elseif isempty(v)
            str = '[]';
        elseif isnumeric(v) || islogical(v)
            str = mat2str(v, 3);
        elseif iscellstr(v)
            str = ['{', strjoin(v, ','), '}'];
        else
            error('Could not convert struct field value');
        end
    end        

end