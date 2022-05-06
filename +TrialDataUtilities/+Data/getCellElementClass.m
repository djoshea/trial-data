function cls = getCellElementClass(dataCell, default)
    if nargin < 2
        default = 'double';
    end
    if isempty(dataCell)
        cls = default;
    elseif ~iscell(dataCell)
        cls = class(dataCell);
    else
        nonEmpty = ~cellfun(@isempty, dataCell);
        if ~any(nonEmpty)
            cls = default; % assume double if no values found
        else
            first = find(nonEmpty, 1);
            cls = class(dataCell{first});
        end
    end
 end
 