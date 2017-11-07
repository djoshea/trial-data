function fnames = orderNevFilesNumerically(nevNames)
% returns nevNames ordered by numerical suffix,
% specifically sort by monkeyInitial ascending, dateString ascending, number ascending, protocol ascending

[~, sortInds] = sortStructArray(parseNevName(nevNames), {'monkeyInitial', 'dateString', 'number', 'protocol'});

fnames = nevNames(sortInds);

end

function [S, finalSortInds] = sortStructArray(S, fldNames)

    if isempty(S)
        S = [];
        finalSortInds = [];
        return;
    end

    finalSortInds = 1:length(S);

    % sort sequentially in reverse order (so that first field specifed is
    % absolutely in order, then by second field, then by third field)
    for iFld = length(fldNames):-1:1
        if ischar(S(1).(fldNames{iFld}))
            vals = {S(finalSortInds).(fldNames{iFld})};
        else
            vals = [S(finalSortInds).(fldNames{iFld})];
        end
        [~, sortInds] = sort(vals);
        finalSortInds = finalSortInds(sortInds);
    end

    % apply the sort
    S = S(finalSortInds);

end