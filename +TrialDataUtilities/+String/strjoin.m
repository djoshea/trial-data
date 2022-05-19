function str = strjoin(strCell, join)
    % str = strjoin(strCell, join)
    % creates a string by concatenating the elements of strCell, separated by the string
    % in join (default = ', ')
    %
    % e.g. strCell = {'a','b'}, join = ', ' [ default ] --> str = 'a, b'

    if isempty(strCell)
        str = "";
        return;
    end
    if ischar(strCell)
        strCell = {strCell};
    end
    if nargin < 2
        join = ', ';
    end

    if iscell(strCell) && ~iscellstr(strCell)
        strCell = string(strCell);
    end

    assert(isstring(strCell) || iscellstr(strCell));
    if isstring(strCell) || isstring(join)
        asString = true;
    else
        asString = false;
    end
    

    if isempty(strCell)
        str = '';
    elseif ischar(strCell)
        str = strCell;
    else
        if isnumeric(strCell) || islogical(strCell)
            % convert numeric vectors to strings
            strCell = arrayfun(@num2str, strCell, 'UniformOutput', false);
        elseif iscell(strCell)
            for i = 1:numel(strCell)
                if isnumeric(strCell{i})
                    strCell{i} = num2str(strCell{i});
                end
                assert(ischar(strCell{i}) || isstring(strCell{i}), 'Contents of strCell must be strings');
            end
        elseif isstring(strCell)
            strCell = cellstr(strCell);
        end
        join = char(join);
        str = cellfun(@(str) [str join], strCell, ...
            'UniformOutput', false);
        str = [str{:}]; 
        str = str(1:end-length(join));
    end
    if asString
        str = string(str);
    end
end
