function newClass = cellDetermineCommonClass(data, origClass)
    if nargin > 1
        newClass = origClass;
    else
        newClass = '';
    end
    nonEmpty = ~cellfun(@isempty, data);
    if ~any(nonEmpty)
        newClass = origClass;
    else
        classes = cellfun(@class, data(nonEmpty), 'UniformOutput', false);
        newClass = TrialDataUtilities.Data.determineCommonClass(classes{:}, newClass);
    end
end