function [joinIdxInPre, nextIdxInPost, dataCat] = computeBestJoinPoint(dataPre, dataPost, nOverlap, varargin)
% time along dim 2, bases along dim 1, joinIdx will have size of dims 3 and
% beyond.
%
% If either of joinAfterIndexPre or joinBeforeIndexPost is NaN, you can specify joinWithin to indicate, 
% splice within joinWithin samples of the specified bound. Both cannot be NaN though.


    p = inputParser();
    p.addParameter('joinAfterIndexPre', -Inf, @isscalar);
    p.addParameter('joinBeforeIndexPost', Inf, @isscalar);
    p.addParameter('joinAfterIndexPost', -Inf, @isscalar);
    p.addParameter('joinBeforeIndexPre', Inf, @isscalar);
    p.addParameter('joinAfterWithin', NaN, @isscalar);
    p.addParameter('joinBeforeWithin', NaN, @isscalar);
    p.addParameter('commonJoinAcrossTrajectories', false, @islogical);
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    nPre =  size(dataPre, 2);
    nPost = size(dataPost, 2);
    showPlot = p.Results.showPlot;

    sz = size(dataPre);
    nTraj = prod(sz(3:end));
    
    function v = rep_ntraj(v)
        if isscalar(v)
            v = repmat(v, nTraj, 1);
        else
            assert(numel(v) == nTraj);
        end
    end
    joinAfterIndexPre = rep_ntraj(p.Results.joinAfterIndexPre);
    joinBeforeIndexPost = rep_ntraj(p.Results.joinBeforeIndexPost);
    joinAfterIndexPost = rep_ntraj(p.Results.joinAfterIndexPost);
    joinBeforeIndexPre = rep_ntraj(p.Results.joinBeforeIndexPre);
    joinAfterWithin = rep_ntraj(p.Results.joinAfterWithin);
    joinBeforeWithin = rep_ntraj(p.Results.joinBeforeWithin);
    
    [joinIdxInPre, nextIdxInPost] = deal(nan(nTraj, 1));
        
    % concatenate the data for each trajectory
    szCat = sz;
    szCat(2) = nPre + nPost - min(nOverlap(:));
    dataCat = nan(szCat);
    
    % quick debug plotting
%     clf;
%     plot(sum(dataPre, [1 3]))
%     hold on
%     plot((1:size(dataPost, 2)) + nPre - nOverlap, sum(dataPost, [1 3]));
    
    if p.Results.commonJoinAcrossTrajectories
        nOverlap = nOverlap(1);
        X = dataPre(:, end-nOverlap+1:end, :);
        Y = dataPost(:, 1:nOverlap, :);
        
%         nCat = nPre + nPost - nOverlap;
        preIndexX = nPre-nOverlap + 1: nPre;
        postIndexY = 1:nOverlap;
        joinAfterIndexPre = max(joinAfterIndexPre);
        joinBeforeIndexPost = min(joinBeforeIndexPost);
        joinBeforeWithin = min(joinBeforeWithin);
        joinAfterWithin = min(joinAfterWithin);

        % take a squared distance between the overlapping points
        cost = sum(sum((X - Y).^2, 1, 'omitnan'), 3, 'omitnan');

        valid = computeValid(preIndexX, postIndexY, joinAfterIndexPre, joinBeforeIndexPost, ...
            joinBeforeIndexPre, joinAfterIndexPost, joinBeforeWithin, joinAfterWithin);     
        cost(~valid) = Inf;

        [~, joinIdxInOverlap] = min(cost);

        joinIdxInPre(:) = preIndexX(joinIdxInOverlap);
        nextIdxInPost(:) = postIndexY(joinIdxInOverlap) + 1;

        insert = cat(2, dataPre(:, 1:joinIdxInPre(1), :), dataPost(:, nextIdxInPost(1):end, :));
        dataCat(:, 1:size(insert, 2), :) = insert;
        
    else
        % per-trajectory joining
        if isscalar(nOverlap)
            nOverlap = repmat(nOverlap, nTraj, 1);
        end
        
        % invalidate too early or too late on a per-trajectory basis
        for c = 1:nTraj
            % extract the overlapping pieces as N x nOverlap
            X = dataPre(:, end-nOverlap(c)+1:end, c);
            Y = dataPost(:, 1:nOverlap(c), c);

%             nCat = nPre + nPost - nOverlap(c);
            preIndexX = nPre-nOverlap(c) + 1: nPre;
            postIndexY = 1:nOverlap(c);

            % take a squared distance between the overlapping points
            cost = sum((X - Y).^2, 1, 'omitnan');

            valid = computeValid(preIndexX, postIndexY, joinAfterIndexPre(c), joinBeforeIndexPost(c), ...
                joinBeforeIndexPre(c), joinAfterIndexPost(c), joinBeforeWithin(c), joinAfterWithin(c));
            cost(~valid) = Inf;

            [~, joinIdxInOverlap] = min(cost);

            joinIdxInPre(c) = preIndexX(joinIdxInOverlap);
            nextIdxInPost(c) = postIndexY(joinIdxInOverlap) + 1;

            insert = cat(2, dataPre(:, 1:joinIdxInPre(c), c), dataPost(:, nextIdxInPost(c):end, c));
            dataCat(:, 1:size(insert, 2), c) = insert;
        end
    end
    
    if showPlot
        figure(); clf;
        if isscalar(nOverlap)
            nOverlap = repmat(nOverlap, nTraj, 1);
        end
        for c = 1:nTraj
            tPre = 1:nPre;
            tPost = (nPre+1:nPre+nPost) - nOverlap(c);
            plot(tPre(1:joinIdxInPre(c)), dataPre(1, 1:joinIdxInPre(c), c), 'k-');
            hold on;
            plot(tPost(nextIdxInPost(c):end), dataPost(1, nextIdxInPost(c):end, c), 'r-');
        end
    end
end

function valid = computeValid(preIndex, postIndex, joinAfterIndexPre, joinBeforeIndexPost, joinBeforeIndexPre, joinAfterIndexPost, joinBeforeWithin, joinAfterWithin)
    valid = true(numel(preIndex), 1);
    if ~isnan(joinAfterIndexPre)
        valid(preIndex < min(joinAfterIndexPre)) = false;
        if ~any(valid)
            error('joinAfterIndexPre is too large');
        end
    end
    if ~isnan(joinBeforeIndexPost)
        valid(postIndex > max(joinBeforeIndexPost)) = false;
        if ~any(valid)
            error('joinBeforeIndexPost is too small');
        end
    end
    if ~isnan(joinBeforeIndexPre)
        valid(preIndex > max(joinBeforeIndexPre)) = false;
        if ~any(valid)
            error('joinBeforeIndexPre is too small');
        end
    end
    if ~isnan(joinAfterIndexPost)
        valid(postIndex < min(joinAfterIndexPost)) = false;
        if ~any(valid)
            error('joinAfterIndexPost is too large');
        end
    end
    
    if ~isnan(joinBeforeWithin)
        % assume joinBefore is specified, we need to be in one of the last joinBeforeWithin samples that precede
        assert(~isnan(joinAfterWithin));
        indFirst = find(valid, 1, 'first');
        indWithin = min(numel(valid), indFirst + joinBeforeWithin);
        valid(indWithin+1:end) = false;
    elseif ~isnan(joinAfterWithin)
        % assume joinAfter is specified, we need to splice in one of the first joinAfterWithin samples that follow
        indLast = find(valid, 1, 'first');
        indWithin = min(numel(valid), indLast + joinAfterWithin);
        valid(indWithin+1:end) = false;
    end
end
