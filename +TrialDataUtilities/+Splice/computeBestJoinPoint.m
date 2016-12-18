function [joinIdxInPre, nextIdxInPost, dataCat] = computeBestJoinPoint(dataPre, dataPost, nOverlap, joinAfterIndexPre, joinBeforeIndexPost)
% time along dim 2, bases along dim 1, joinIdx will have size of dims 3 and
% beyond. joinIndex is guaranteed to correspond to 

    nPre =  size(dataPre, 2);
    nPost = size(dataPost, 2);
        
    if nargin < 4
        joinAfterIndexPre = 0;
    end
    if nargin < 5
        joinBeforeIndexPost = Inf;
    end
    
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
    
%     showPlot = false;
%     if showPlot
%         clf;
%         tPre = 1:nPre;
%         tPost = (nPre+1:nPre+nPost) - nTimepointsOverlap;
%         for c = 1:nTraj
%             plot(tPre(1:joinIdxInPre(c)), dataPre(1, 1:joinIdxInPre(c), c), 'k-');
%             hold on;
%             plot(tPost(nextIdxInPost(c):end), dataPost(1, nextIdxInPost(c):end, c), 'r-');
%         end
%     end
end
