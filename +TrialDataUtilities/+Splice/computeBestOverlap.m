function nTimepointsOverlap = computeBestOverlap(dataPre, dataPost, minOverlap, maxOverlap, varargin)
% nTimepointsOverlap = computeBestOverlap(dataPre, dataPost, minOverlap, maxOverlap)
%
% figure out the amount of overlap along dimension 2 (time) between the
% right edge of dataPre and the left edge of dataPost. This is done by
% minimizing mean pointwise squared difference over the set of overlapping
% points. The range of overlap can be set with minOverlap and maxOverlap.
%
% dataPre and dataPost should not have lagging or leading NaNs
%    
    
    p = inputParser();
    p.addParameter('commonAcrossTrajectories', false, @islogical); % enforce all trajectories (along dim 3) have the same nTimepointsOverlap
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    sz = size(dataPre);
    nTraj = prod(sz(3:end));
    nTrajOrig = nTraj;
    
    if p.Results.commonAcrossTrajectories
        % treat all other dimensions as just more along dimension 1
        if ~ismatrix(dataPre)
            otherDims = TensorUtils.otherDims(size(dataPre), 2);
            dataPre = TensorUtils.reshapeByConcatenatingDims(dataPre, {otherDims, 2});
            dataPost = TensorUtils.reshapeByConcatenatingDims(dataPost, {otherDims, 2});
            nTraj = 1;
        end
    end
    
    % we need to explore distance between overlapping trajectories
    nPre = size(dataPre, 2);
    nPost = size(dataPost, 2);
    maxOverlap = min([nPre nPost maxOverlap]);
    minOverlap = max(1, minOverlap);
    if minOverlap > maxOverlap
        error('Min overlap exceeds max overlap or length of data');
    end

    nTimepointsOverlap = nan(nTraj, 1);
    
    for c = 1:nTraj
        % grab the touching edges of dataPre and dataPost
        X = dataPre(:, end-maxOverlap+1:end, c);
        Y = dataPost(:, 1:maxOverlap, c);

        % dist doesn't play with NaNs so we strip any rows from X and Y
        % that have NaNs
        keep = ~any(isnan(X), 2) & ~any(isnan(Y), 2);
        X = X(keep, :);
        Y = Y(keep, :);
        
        % X and Y are (bases x maxOverlap time) ,  or (bases * conditions x maxOverlap time) if commonAcrossTrajectories
        % dist(i, j) is the squared distance between X(:, i) and Y(:, j) (note transposed X and Y in call to pdist2 below)
        % dist is maxOverlap x maxOverlap
        dist = pdist2(X', Y', 'squaredeuclidean');

        % the sum of the subdiagonals of dist corresponds to the total distance
        % between the overlapping fragments of X and Y
        totalDistByOverlap = sum(spdiags(dist, -maxOverlap+minOverlap:0), 1);
        nOverlap = minOverlap:maxOverlap;

        meanDistOverlap = totalDistByOverlap ./ nOverlap;
        [~, idx] = min(meanDistOverlap);
        nTimepointsOverlap(c) = nOverlap(idx);
        
        if p.Results.showPlot
            stem(nOverlap, meanDistOverlap)
            hold on;
        end
    end

    if p.Results.commonAcrossTrajectories
        nTimepointsOverlap = repmat(nTimepointsOverlap, nTrajOrig, 1);
    end
end