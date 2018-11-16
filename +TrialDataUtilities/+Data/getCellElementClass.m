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
 