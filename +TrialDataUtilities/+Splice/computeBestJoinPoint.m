function [joinIdxInPre, nextIdxInPost, dataCat] = computeBestJoinPoint(dataPre, dataPost, nOverlap, varargin)
% time along dim 2, bases along dim 1, joinIdx will have size of dims 3 and
% beyond. joinIndex is guaranteed to correspond to 

    p = inputParser();
    p.addParameter('joinAfterIndexPre', 0, @isscalar);
    p.addParameter('joinBeforeIndexPost', Inf, @isscalar);
    p.addParameter('commonJoinAcrossTrajectories', false, @islogical);
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    nPre =  size(dataPre, 2);
    nPost = size(dataPost, 2);
    joinAfterIndexPre = p.Results.joinAfterIndexPre;
    joinBeforeIndexPost = p.Results.joinBeforeIndexPost;
    showPlot = p.Results.showPlot;

    sz = size(dataPre);
    nTraj = prod(sz(3:end));
    
    if isscalar(joinAfterIndexPre)
        joinAfterIndexPre = repmat(joinAfterIndexPre, nTraj, 1);
    end
    if isscalar(joinBeforeIndexPost)
        joinBeforeIndexPost = repmat(joinBeforeIndexPost, nTraj, 1);
    end
    
    [joinIdxInPre, nextIdxInPost] = deal(nan(nTraj, 1));
        
    % concatenate the data for each trajectory
    szCat = sz;
    szCat(2) = nPre + nPost - min(nOverlap(:));
    dataCat = nan(szCat);
    
    if p.Results.commonJoinAcrossTrajectories
        nOverlap = nOverlap(1);
        X = dataPre(:, end-nOverlap+1:end, :);
        Y = dataPost(:, 1:nOverlap, :);
        
        nCat = nPre + nPost - nOverlap;
        preIndexX = nPre-nOverlap + 1: nPre;
        postIndexY = 1:nOverlap;

        % take a squared distance between the overlapping points
        cost = sum(sum((X - Y).^2, 1, 'omitnan'), 3, 'omitnan');

        valid = true(nOverlap, 1);
        valid(preIndexX < min(joinAfterIndexPre)) = false;
        if ~any(valid)
            error('joinAfterIndexPre is too large');
        end
        valid(postIndexY > max(joinBeforeIndexPost)) = false;
        if ~any(valid)
            error('joinBeforeIndexPost is too small');
        end
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

            nCat = nPre + nPost - nOverlap(c);
            preIndexX = nPre-nOverlap(c) + 1: nPre;
            postIndexY = 1:nOverlap(c);

            % take a squared distance between the overlapping points
            cost = sum((X - Y).^2, 1, 'omitnan');

            valid = true(nOverlap(c), 1);
            valid(preIndexX < joinAfterIndexPre(c)) = false;
            if ~any(valid)
                error('joinAfterIndexPre is too large');
            end
            valid(postIndexY > joinBeforeIndexPost(c)) = false;
            if ~any(valid)
                error('joinBeforeIndexPost is too small');
            end
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
