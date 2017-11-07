function [eachAsMatchInB eachBsMatchInA] = matchStructArraysbyField(A, B, fields, nameA, nameB, varargin)
% [eachAsMatchInB eachBsMatchInA] = matchStructArraysbyField(A, B, fields, nameA, nameB)

if ~exist('nameA', 'var') || ~exist('nameB', 'var')
    nameA = '';
    nameB = '';
end

def.printStats = true;
def.warnThreshA = 0.95;
def.warnThreshB = 0.95;
assignargs(def, varargin);

if ischar(fields)
    fields = {fields};
end

% stats on how the matching goes
multipleMatchesFound = 0;
eachAsMatchInB = nan(length(A), 1);
eachBsMatchInA = nan(length(B), 1);

% loop over B, find matches in p
for iB = 1:length(B)
    % a match A(iA) must sastisfy A(iA).(fields{i}) == B(iB).(fields{i}) for all i
    matchesInA = true(1,length(A));
    for ifld = 1:length(fields)
        matchesInA = matchesInA & B(iB).(fields{ifld}) == [A.(fields{ifld})];
    end

    % is there a match?
    indInA = find(matchesInA);
    if ~isempty(indInA)
        if numel(indInA) > 1
            multipleMatchesFound = multipleMatchesFound + 1;
        end
        eachBsMatchInA(iB) = indInA(1);
        eachAsMatchInB(indInA(1)) = iB;
    end
end

if printStats & ~isempty(nameA) & ~isempty(nameB)
    % calculate summary stats on matching
    numBTotal = length(B);
    numBMatched = nnz(~isnan(eachBsMatchInA));
    numATotal = length(A);
    numAMatched = nnz(~isnan(eachAsMatchInB));

    % include a warning if > some percentage aren't matched
    if(numBMatched >= warnThreshB*numBTotal)
        incompleteStrB = ' ';
    else
        incompleteStrB = '*';
    end
    if(numAMatched >= warnThreshA*numATotal)
        incompleteStrA = ' ';
    else
        incompleteStrA = '*';
    end

    if multipleMatchesFound
        warningStr = sprintf('[ warning: %d ambiguous matches encountered! ]', multipleMatchesFound);
    else
        warningStr = '';
    end

    fprintf('\t\b\b%s %d / %d %s, %s %d / %d %s matched %s\n', ...
         incompleteStrA, numAMatched, numATotal, nameA, ...
         incompleteStrB, numBMatched, numBTotal, nameB, warningStr);
end
