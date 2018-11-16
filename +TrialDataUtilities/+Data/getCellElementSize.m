function sz = getCellElementSize(dataCell)
    if isempty(dataCell)
        sz = [];
    elseif ~iscell(dataCell)
        sz = size(dataCell);
    else
        nonEmpty = ~cellfun(@isempty, dataCell);
        if ~any(nonEmpty)
            sz = [];
        else
            first = find(nonEmpty, 1);
            sz = size(dataCell{first});
        end
    end
end