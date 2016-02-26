function contents = getFirstNonEmptyCellContents(data)
    contents = [];
    for i = 1:numel(data)
        if isempty(data{i}), continue; end
        contents = data{i};
        break;
    end
end