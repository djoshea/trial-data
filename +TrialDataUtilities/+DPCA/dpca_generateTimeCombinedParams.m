function [combined, names] = dpca_generateTimeCombinedParams(nDims, dimNames)
% builds combinedParams argument assuming dim 1 is time and all pure param
% mixtures should be added to param x time.

import(getPackageImportString);

    if nargin < 2
        dimNames = arrayfun(@(i) sprintf('Param %d', i), ...
            1:nDims, 'UniformOutput', false);
    end

    sets = subsets(2:nDims+1);
    %sets = subsets(2:nDims);

    combined = cell(numel(sets) + 1, 1);
    names = cell(numel(sets) + 1, 1);
    combined{1} = {1};
    names{1} = 'time';

    for i = 1:numel(sets)
        pure = sets{i};
        withTime = [1, pure];
        combined{i+1} = {pure, withTime};
        names{i+1} = strjoin(dimNames(pure-1), ' x ');
    end

end

function S = subsets(X)

    % S = subsets(X) returns a cell array of all subsets of vector X apart
    % from the empty set. Subsets are ordered by the number of elements in
    % ascending order.
    %
    % subset([1 2 3]) = {[1], [2], [3], [1 2], [1 3], [2 3], [1 2 3]}

    d = length(X);
    pc = dec2bin(1:2^d-1) - '0';
    [~, ind] = sort(sum(pc, 2));
    pc = fliplr(pc(ind,:));
    for i=1:length(pc)
        S{i} = X(find(pc(i,:)));
    end
end