function str = structToString(s, varargin)
    p = inputParser();
    p.addOptional('separator', ' ', @ischar);
    p.addParameter('includeFieldNames', true, @islogical);
    p.parse(varargin{:});
    
    separator = p.Results.separator;
    includeFieldNames = p.Results.includeFieldNames;

    fields = fieldnames(s);
    if isempty(fields)
        str = '';
        return;
    end
    vals = structfun(@convertToString, s);

    if includeFieldNames
        str = strjoin(cellfun(@(fld, val) [fld '=' val], fields, vals, 'UniformOutput', false), separator);
    else
        str = strjoin(vals, separator);
    end
    
    return;
    
    function str = convertToString(v)
        if ischar(v)
            str = v;
        elseif isempty(v)
            str = '[]';
        elseif isnumeric(v) || islogical(v)
            str = mat2str(v, 3);
        elseif iscellstr(v)
            str = ['{', strjoin(v, ','), '}'];
        else
            error('Could not convert struct field value');
        end
        
        str = {str};
    end        

end