function [fullList, names] = dpca_generateTimeCombinedParams(dims, varargin)
% function [fullList, names] = dpca_generateTimeCombinedParams(dims, varargin)
% builds combinedParams argument for dpca, which specifies the covariate and
% covariate interactions to combine into f
%
% dims: 
% 
% 'combine': optional cell of cells, each element is a cell vector specifying
%   the covariates or sets of interacting covariates to be defined. Each element of the
%   inner cells is a covariate index. These indices should be 1-indexed, and 0 
%   should be used to represent time
%
% 'combineEachWithTime': true (default) or false. If true, all
%   covariates will be automatically combined with the time+covariate
%   interaction. this is equivalent to adding to 'combine':
%     { {1, {0 1}}, {2, {0 2}, ..., {nCov, {0 nCov} }
%
% 'covariateNames': nCov cellstr. optional, specifies the names for each
%   covariate.

    import(getPackageImportString);
    
    p = inputParser();
    p.addParameter('combine', {}, @(x) isempty(x) || iscell(x));
    p.addParameter('combineEachWithTime', true, @islogical);
    p.addParameter('covariateNames', {}, @iscellstr);
    p.parse(varargin{:});
    
    dimListsToCombineList = p.Results.combine;
    
    if isempty(dimListsToCombineList)
        dimListsToCombineList = {};
    end
        
    if ~iscell(dimListsToCombineList)
        error('Dim Lists to combine must be cell of cells of dim indices');
    end
    
    % ensure all contents are cell arrays using num2cell so [1 2] becomes
    % {1, 2}
    for i = 1:numel(dimListsToCombineList)
        if ~iscell(dimListsToCombineList{i})
            error('Dim Lists to combine must be cell of cells of dim indices');
        end
    end
    
    dimListsToCombineList = makecol(dimListsToCombineList);

    % get unique dims inside
    flattened = cellvec(numel(dimListsToCombineList));
    for i = 1:numel(dimListsToCombineList)
        inner = cellfun(@makecol, dimListsToCombineList{i}, 'UniformOutput', false) ;
        flattened{i} = cat(1, inner{:});
    end
    uniqueDims = unique(cat(1, flattened{:}));
    
    assert(all(ismember(uniqueDims, dims)), 'Some dims in lists not found in dims');
    
    % generate names for each dim if not provided
    nDims = max(dims);
    if isempty(p.Results.covariateNames)
        covariateNames = arrayfun(@(i) sprintf('Axis %d', i), ...
            1:nDims, 'UniformOutput', false);
    else
        covariateNames = p.Results.covariateNames;
        assert(numel(covariateNames) >= nDims, 'Provided covariateNames needs at least %d entries', nDims);
    end
    covariateNames = makecol(covariateNames);
    
    % add dim, time+dim combinations to the list 
    if p.Results.combineEachWithTime
        % add all the {0, dim} combinations
        timeList = cellfun(@(dim) {dim; [0; dim]}, subsets(dims), 'UniformOutput', false); 
        dimListsToCombineList = [timeList; dimListsToCombineList];
    end
    
    % remove singular lists (nothing being combined)
    nToCombine = cellfun(@numel, dimListsToCombineList);
    dimListsToCombineList = dimListsToCombineList(nToCombine > 1);
    
    % generate full list
    fullList = num2cell(subsets(union(0, dims)));
    
    for iC = 1:numel(dimListsToCombineList)
        % search for any superset of the elements of the combine list
        combine = dimListsToCombineList{iC};
        searchFor = supersets(combine);
        
        % search for rows of full list that match
        matches = falsevec(numel(fullList));
        for iV = 1:numel(searchFor)
            value = searchFor{iV};
            rowMatchesFn = @(list) any(cellfun(@(v) isequal(v, value), list));
            matches = matches | cellfun(rowMatchesFn, fullList);
        end
       
        % combine all the matches
        combinedMatches = uniqueCell(cat(1, fullList{matches}));
        
        % remove the matching rows and add the combined row
        fullList = [fullList(~matches); {combinedMatches}];
    end
    
    if nargout > 1
        covariateNames = [{'time'}; covariateNames];
        names = cellvec(numel(fullList));
        for iF = 1:numel(fullList)
            pieceStr = cellvec(numel(fullList{iF}));
            for iP = 1:numel(fullList{iF})
                idx = fullList{iF}{iP} + 1;
                if p.Results.combineEachWithTime && ~isequal(idx, 1) && ismember(1, idx)
                    pieceStr{iP} = '';
                else
                    pieceStr{iP} = strjoin(covariateNames(idx), ' x ');
                end
            end
            pieceStr = pieceStr(~cellfun(@isempty, pieceStr));
            names{iF} = strjoin(pieceStr, ', ');
        end
    end
    
    % change to being time = 1, dim 1 = 2 indexed
    for iF = 1:numel(fullList)
        for iP = 1:numel(fullList{iF})
            fullList{iF}{iP} = fullList{iF}{iP} + 1;
        end
    end
end

function S = subsets(X)

    % S = subsets(X) returns a cell array of all subsets of vector X apart
    % from the empty set. Subsets are ordered by the number of elements in
    % ascending order.
    %
    % subset([1 2 3]) = {[1], [2], [3], [1 2], [1 3], [2 3], [1 2 3]}

    X = makecol(X);    
    d = length(X);
    pc = dec2bin(1:2^d-1) - '0';
    [~, ind] = sort(sum(pc, 2));
    pc = fliplr(pc(ind,:));
    S = cellvec(length(pc));
    for iP=1:length(pc)
        S{iP} = makecol(X(logical(pc(iP,:))));
    end
end

function S = supersets(X)
    % S = supersets(X) returns a cell array of all supersets of the cell contents of X apart
    % from the empty set.
    %
    % superset({1, [1 2], [1 3]}) = {[1], [1 2], [1 3], [1 2 3]}

    X = cellfun(@makecol, X, 'UniformOutput', false);
    idxSubsets = subsets(1:numel(X));
    
    S = cellvec(numel(idxSubsets));
    for iS = 1:numel(idxSubsets)
        S{iS} = unique(cat(1, X{idxSubsets{iS}}));
    end
    
    S = uniqueCell(S);
end

function [B, I, J] = uniqueCell(A, varargin)
    removeEmpty = false;
    assignargs(varargin);

    B = cell(0, 1);
    I = [];
    J = zeros(numel(A),1);
    for iA = 1:numel(A)
        if removeEmpty && isempty(A{iA})
            I(iA) = NaN;
            J(iA) = NaN;
            continue;
        end
        
        idx = find(cellfun(@(b) isequal(A{iA}, b), B), 1, 'first');
        if isempty(idx)
            B{end+1} = A{iA};
            I(end+1) = iA;
            J(iA) = numel(B); 
        else
            J(iA) = idx;
        end
    end
    B = makecol(B);
end