function dataSpliced = spliceSingleTrajectoryPair(dataPre, dataPost, varargin)
    p = inputParser();
    % p.addParameter('tPre', 1:size(dataPre, 2), @isvector);
    % p.addParameter('tPost', size(dataPre, 2) + (1:size(dataPost, 2)), @isvector);
%     p.addParameter('dim', 2, @isscalar);
    p.addParameter('searchPre', 100, @isscalar);
    p.addParameter('searchPost', 100, @isscalar);
    p.addParameter('splineIgnore', 5, @isscalar);
    p.addParameter('splineFit', 15, @isscalar);
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});

    windowPre = p.Results.searchPre;
    windowPost = p.Results.searchPost;
    splineReplaceWindow = p.Results.splineIgnore;
    splineIncludeWindow = p.Results.splineFit;
    
    nPre = size(dataPre, 2);
    nPost = size(dataPost, 2);

    X = fliplr(dataPre(:, end-windowPre+1:end)); % goes from end --> end-windowPre+1
    Y = dataPost(:, 1:windowPost);

    % X = X ./ nanstd(dataPre(:));
    % Y = Y ./ nanstd(dataPost(:));

    dist = pdist2(X', Y', 'squaredeuclidean');

    % scale to [0 1]
    dist = (dist - nanmin(dist(:))) ./ (nanmax(dist(:)) - nanmin(dist(:)));
    % [Xadd,Yadd] = ndgrid(1:windowPre, 1:windowPost);
    % dist = dist + (Xadd + Yadd) ./ (windowPre + windowPost);

    % dist is windowPre x windowPost
    [~, idxMin] = nanmin(dist(:));
    [idxPre, idxPost] = ind2sub(size(dist), idxMin);

    % ensure at least 1 point is dropped
    idxPre = max(idxPre, 2);

    dataPreKeep = dataPre(:, 1:end-idxPre+1);
    dataPostKeep = dataPost(:, idxPost:end);
    nPreKeep = size(dataPreKeep, 2);
    nPostKeep = size(dataPostKeep, 2);
    
    % take splineIncludeWindow points to either side for fitting
    % excluding the last splineReplaceWindow points to either side
    dataFitSpline = cat(2, dataPreKeep(:, end-splineReplaceWindow-splineIncludeWindow+1:end-splineReplaceWindow), ...
        dataPostKeep(:, splineReplaceWindow+1:splineReplaceWindow+splineIncludeWindow));

    % use technique from cscvn but use smoothing spline rather than
    % interpolation
    if size(dataFitSpline, 1)==1
        dt = 0;
    else
        dt = sum((diff(dataFitSpline').^2).'); 
    end
    tFitSpline = cumsum([0,dt.^(1/4)]);
    % cs = csaps(tFitSpline, dataFitSpline, 'variational');
    cs = csaps(tFitSpline, dataFitSpline);

    % evaluate on even time grid that matches original data

    nOverwritePre = nPre - nPreKeep + splineIncludeWindow + splineReplaceWindow;
    nOverwritePost = nPost - nPostKeep + splineIncludeWindow + splineReplaceWindow;
    tEvalSpline = linspace(min(tFitSpline), max(tFitSpline), nOverwritePre + nOverwritePost);
    dataEvalSpline = fnval(cs, tEvalSpline);
 
    dataPreSmoothed = dataPre;
    dataPreSmoothed(:, end-nOverwritePre+1:end) = dataEvalSpline(:, 1:nOverwritePre);

    dataPostSmoothed = dataPost;
    dataPostSmoothed(:, 1:nOverwritePost) = dataEvalSpline(:, nOverwritePre+1:end);

    dataSpliced = cat(2, dataPreSmoothed, dataPostSmoothed);

    if p.Results.showPlot
        % plot traces and traces ends in red
        plot3(dataPre(1, :), dataPre(2, :), dataPre(3, :), 'k.');
        hold on;
        plot3(dataPre(1, end), dataPre(2, end), dataPre(3, end), 'ro', 'MarkerFaceColor', 'r');

        plot3(dataPost(1, :), dataPost(2, :), dataPost(3, :), 'k.');
        plot3(dataPost(1, 1), dataPost(2, 1), dataPost(3, 1), 'ro', 'MarkerFaceColor', 'r');

        % plot last timepoints retained for splicing in green
        plot3(dataPreKeep(1, end), dataPreKeep(2, end), dataPreKeep(3, end), 'go', 'MarkerFaceColor', 'g');
        plot3(dataPostKeep(1, 1), dataPostKeep(2, 1), dataPostKeep(3, 1), 'go', 'MarkerFaceColor', 'g');

        % plot splice results
        plot3(dataSpliced(1, :), dataSpliced(2, :), dataSpliced(3, :), 'r');

        set(findall(gca, 'Type', 'line'), 'Clipping', 'off');
    end
    
end