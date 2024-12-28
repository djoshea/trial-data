function [dataSpliced, info, opts] = splicePair(dataPre, dataPost, varargin)
% Splicing together high dimensional state trajectories along the time
% dimension dim 2. dim 1 is bases. dim 3 and beyond are different trajectories,
% presumably from different conditions.
%
% First, we determine the amount of temporal overlap
% between the two trajectories, which works better if there is sufficient
% real overlap between the two. This can be constrained using minOverlap
% and maxOverlap. It will be a common value for all trajectories. 
%
% Next, we find the best "splice point" at which to stop taking data from
% dataPre and start taking data from dataPost. The range of possible splice 
% points here can be constrained using searchPre and searchPost, which both
% default to Inf (any splice point is valid).
%
% Lastly, optionally, we chop off some data to either side and replace it
% with a spline interpolation.
%
% Returns:
%   dataSpliced
%   info - information about where the splicing took place
%   opts - a structure that can be passed in to ensure splicing with new
%      data occurs at the exact same locations. 

    p = inputParser();
    
    %%%
    % Dimension parameteris
    %%%

    p.addParameter('basisDim', 1, @isscalar);
    p.addParameter('spliceDim', 2, @isscalar);

    %%%
    % OVERLAP Parameters
    %%%
    
    % constrain the number of timepoints that could be overlapping between
    % the two trajectories
    p.addParameter('commonOverlapAcrossTrajectories', true, @islogical);
        
    % EITHER specify exactly with this:
    p.addParameter('nTimepointsOverlap', [], @(x) isempty(x) || isnumeric(x));
    % OR specify a range to search over
    p.addParameter('minOverlap', 0, @isscalar); 
    p.addParameter('maxOverlap', Inf, @isscalar);
    
    %%%
    % Join Parameters
    %%%
    
    % common timepoint where the splicing occurs, true requires
    % commonOverlapAcrossTrajectories to be true also
    p.addParameter('commonJoinAcrossTrajectories', false, @islogical);
    
    % EITHER specify these join points directly (as scalar or matrix)
    p.addParameter('joinIdxInPre', [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('nextIdxInPost', [], @(x) isempty(x) || isnumeric(x));
    
    % OR specify the edges of the range we should search over to
    % constrain the search window for the splice break point
    p.addParameter('joinAfterIndexPre', NaN, @(x) isempty(x) || isnumeric(x));
    p.addParameter('joinBeforeIndexPost', NaN, @(x) isempty(x) || isnumeric(x)); % indices to search into the post trajectory to find the splice point
    p.addParameter('joinAfterIndexPost', NaN, @(x) isempty(x) || isnumeric(x));
    p.addParameter('joinBeforeIndexPre', NaN, @(x) isempty(x) || isnumeric(x));
    p.addParameter('joinAfterWithin', NaN, @isnumeric); % indices to search into the post trajectory to find the splice point
    p.addParameter('joinBeforeWithin', NaN, @isnumeric); % indices to search into the post trajectory to find the splice point
    
    % for spline interpolation
    p.addParameter('interpolateMethod', 'spline', @isstringlike);
    p.addParameter('interpIgnoreWindow', 15, @isscalar); % ignore these last points around the trjaectories when fitting the spline
    p.addParameter('interpFitWindow', 30, @isscalar); % include these last points as waypoints for the splines, should be bigger than splineIgnore since splineIgnore will be cut out of the middle of this time window pre and post
    
    %%%
    % PCA preprocessing parameters
    %%%
    
    p.addParameter('usePCA', false, @islogical); % if true, use PCA and splice based on capturing 80% of the variance 
    p.addParameter('nPCs', 6, @(x) isempty(x) || isscalar(x));
    
    % specify both of these to use the projections
    p.addParameter('pcaCoeff', [], @(x) isempty(x) || ismatrix(x)); % specify if you want to manually specify the projection
    p.addParameter('pcaMeans', [],  @(x) isempty(x) || isvector(x));% specify if you want to manaully speicify the means subtracted off
    
    p.addParameter('showPlot', false, @islogical);
    
    p.parse(varargin{:});

    % handle permute to put bases along dim 1, splice along dim 2
    basisDim = p.Results.basisDim;
    spliceDim = p.Results.spliceDim;
    dimConditions = TensorUtils.otherDims(size(dataPre), [basisDim, spliceDim]);
    dataPre = permute(dataPre, [basisDim, spliceDim, dimConditions]);
    dataPost = permute(dataPost, [basisDim, spliceDim, dimConditions]);

    opts.basisDim = basisDim;
    opts.spliceDim = spliceDim;

    nPre = size(dataPre, 2); %#ok<NASGU>
    nPost = size(dataPost, 2);
    %nBases = size(dataPre, 1);
    nTraj = prod(TensorUtils.sizeOtherDims(dataPre, [1 2])); %#ok<NASGU>

    joinIdxInPre = p.Results.joinIdxInPre;
    nextIdxInPost= p.Results.nextIdxInPost;
    joinBeforeIndexPost = p.Results.joinBeforeIndexPost;
    joinAfterIndexPre = p.Results.joinAfterIndexPre;
    joinBeforeIndexPre = p.Results.joinBeforeIndexPre;
    joinAfterIndexPost = p.Results.joinAfterIndexPost;
    joinAfterWithin = p.Results.joinAfterWithin;
    joinBeforeWithin = p.Results.joinBeforeWithin;
    showPlot = p.Results.showPlot;

    % optional PCA preprocessing of the data to make splicing easier
    if p.Results.usePCA
        dataPCA = cat(2, dataPre, dataPost);
        if isempty(p.Results.pcaCoeff) || isempty(p.Results.pcaMeans)
            nPCs = min(p.Results.nPCs, size(dataPCA, 1));
            sWarn = warning('off', 'stats:pca:ColRankDefX'); % we're not using tsquared
            [coeff, ~, ~, ~, ~, mu] = TensorUtils.pcaAlongDim(dataPCA, 1, 'NumComponents', nPCs);
            warning(sWarn);
        else
            coeff = p.Results.pcaCoeff;
            mu = p.Results.pcaMeans;
        end
            
        dataPreProj = TensorUtils.linearCombinationAlongDimension(dataPre - mu, 1, coeff', 'replaceNaNWithZero', true);
        dataPostProj= TensorUtils.linearCombinationAlongDimension(dataPost - mu, 1, coeff', 'replaceNaNWithZero', true);
        
        opts.usePCA = true;
        opts.nPCs = p.Results.nPCs;
        opts.pcaCoeff = coeff;
        opts.pcaMeans = mu;
    else
        coeff = eye(size(dataPre, 1));
        mu = zerosvec(size(dataPre, 1));
        dataPreProj = dataPre;
        dataPostProj = dataPost;
        opts.usePCA = false;
    end
    
    % first, we compute the best amount of temporal overlap of the edges of
    % the trajectories, unless specified as input
    if isempty(p.Results.nTimepointsOverlap)
        nTimepointsOverlap_auto = true;
        nTimepointsOverlap = TrialDataUtilities.Splice.computeBestOverlap(dataPreProj, dataPostProj, ...
            p.Results.minOverlap, p.Results.maxOverlap, 'commonAcrossTrajectories', p.Results.commonOverlapAcrossTrajectories, 'showPlot', showPlot);
    else
        nTimepointsOverlap_auto = false;
        nTimepointsOverlap = p.Results.nTimepointsOverlap;
    end
    opts.commonOverlapAcrossTrajectories = p.Results.commonOverlapAcrossTrajectories;
        
    if isempty(joinIdxInPre) || isempty(nextIdxInPost)
        % then we assume that overlap and use the best splice point on a per
        % trajectory basis (i.e. for each traj along dims 3)
        [opts.joinIdxInPre, opts.nextIdxInPost] = TrialDataUtilities.Splice.computeBestJoinPoint(dataPreProj, dataPostProj, nTimepointsOverlap, ...
            'joinAfterIndexPre', joinAfterIndexPre, 'joinBeforeIndexPost', joinBeforeIndexPost, ...
            'joinBeforeIndexPre', joinBeforeIndexPre, 'joinAfterIndexPost', joinAfterIndexPost, ...
            'joinAfterWithin', joinAfterWithin, 'joinBeforeWithin', joinBeforeWithin, ...
            'commonJoinAcrossTrajectories', p.Results.commonJoinAcrossTrajectories && p.Results.commonOverlapAcrossTrajectories, ...
            'showPlot', showPlot);
    else
        % use the specified splice point (optionally per trajectory)
        szPre = TensorUtils.sizeNDims(dataPre, 3);
        opts.joinIdxInPre = TensorUtils.scalarExpandToSize(joinIdxInPre, szPre(3:end));
        opts.nextIdxInPost = TensorUtils.scalarExpandToSize(nextIdxInPost, szPre(3:end));

        nOverPre = nPre - opts.joinIdxInPre;
        nOverPost = opts.nextIdxInPost - 1;
            
        if ~nTimepointsOverlap_auto
            % check nTimepointsOverlap and joinIdxInPre / nextIdxInPost are compatible
            assert(all(nOverPre + nOverPost == nTimepointsOverlap), "Specified nTimepointsOverlap inconsistent with specified joinIdxInPre, nextIdxInPost");
        else
            % replace nTimepointsOverlap with computed value
            nTimepointsOverlap = nOverPre(1) + nOverPost(1); % must be same for all entries
            assert(all(nTimepointsOverlap == nTimepointsOverlap(1)), "Specified joinIdxInPre and nextIdxInPost induce differing nTimepointsOverlap across conditions");
        end
    end
    opts.nTimepointsOverlap = nTimepointsOverlap; % C
    info.nTimepointsOverlap = nTimepointsOverlap;
    
    info.joinIdxInPre = opts.joinIdxInPre;
    info.nextIdxInPost = opts.nextIdxInPost;
    opts.commonJoinAcrossTrajectories = p.Results.commonJoinAcrossTrajectories;
        
    % do the concatenation here (since the best join point may be done on
    % the PCs)
    dataCat = TrialDataUtilities.Splice.simpleCatJoin(dataPre, dataPost, info.joinIdxInPre, info.nextIdxInPost);    
    T = size(dataCat, 2);
    C = numel(info.joinIdxInPre);
    szCat = size(dataCat);
    
    % build info matrices describing where each point of the concatenated
    % timeseries pulls its data from for pre and post
    [info.idxFromPre, info.idxFromPost] = deal(nan([T, szCat(3:end)]));
    for c = 1:C
        info.idxFromPre(1:info.joinIdxInPre(c), c) = 1:info.joinIdxInPre(c);
        info.idxFromPost(info.joinIdxInPre(c)+1:end, c) = info.nextIdxInPost(c) : nPost;
    end
    
    opts.interpolateMethod = string(p.Results.interpolateMethod);
    opts.interpIgnoreWindow = p.Results.interpIgnoreWindow;
    opts.interpFitWindow = p.Results.interpFitWindow;
    
    % do smooth interpolation at the join
    if strcmp(p.Results.interpolateMethod, "spline")
        dataSpliced = TrialDataUtilities.Splice.interpolateSpline(dataCat, info.joinIdxInPre+1, ...
            p.Results.interpFitWindow, p.Results.interpIgnoreWindow, 'showPlot', showPlot);

    elseif ismember(char(p.Results.interpolateMethod), {'linear', 'nearest', 'next', 'previous', 'cspline', 'pchip', 'pchip', 'cubic', 'v5cubic'})
        % do linear interpolation
        dataSpliced = TrialDataUtilities.Splice.interpolateLinear(dataCat, info.joinIdxInPre+1, ...
            p.Results.interpolateMethod, p.Results.interpFitWindow, p.Results.interpIgnoreWindow, ...
            'showPlot', showPlot);
    elseif ismember(char(p.Results.interpolateMethod), {'', 'none'})
        % do simple concatenation
        dataSpliced = dataCat;
    else
        error('Unknown interpolateMethod %s', p.Results.interpolateMethod);
    end
    
    info.spliceStart = info.joinIdxInPre+1 - p.Results.interpFitWindow - p.Results.interpIgnoreWindow;
    info.spliceStop = info.joinIdxInPre+1 + p.Results.interpFitWindow + p.Results.interpIgnoreWindow;
    
    if showPlot 
        
        c = 1;
        preColor = [0 0 0];
        endColor = 'r';
        startColor = 'g';
        postColor = [0.5 0.5 0.5];
        
        clf
        % plot pre-traces and pre trace ends in red
        plot3(dataPreProj(1, :, c), dataPreProj(2, :, c), dataPreProj(3, :, c), '.', 'Color', preColor);
        hold on;
        plot3(dataPreProj(1, info.joinIdxInPre(c), c), dataPreProj(2, info.joinIdxInPre(c), c), dataPreProj(3, info.joinIdxInPre(c), c), 'o', 'MarkerSize', 15, 'MarkerFaceColor', endColor, 'MarkerEdgeColor', endColor);
        
        % plot pre-traces and pre trace ends in red
        plot3(dataPostProj(1, :, c), dataPostProj(2, :, c), dataPostProj(3, :, c), '.', 'Color', postColor);
        plot3(dataPostProj(1, info.nextIdxInPost(c), c), dataPostProj(2, info.nextIdxInPost(c), c), dataPostProj(3, info.nextIdxInPost(c), c), 'o', 'MarkerSize', 15,  'MarkerFaceColor', startColor, 'MarkerEdgeColor', startColor);

        % plot last timepoints retained for splicing in green
        %plot3(dataPreProj(1, end, iC), dataPreProj(2, end, iC), dataPreProj(3, end, iC), 'go', 'MarkerFaceColor', 'g');
        %plot3(dataPostProj(1, 1, iC), dataPostProj(2, 1, iC), dataPostProj(3, 1, iC), 'go', 'MarkerFaceColor', 'g');
        
        % smooth lines through the portion before and after splicing in light red (start) and light green (post)
        dataPreProjSmooth = smoothTraces(dataPreProj);
        dataPostProjSmooth = smoothTraces(dataPostProj);
        plot3(dataPreProjSmooth(1, 1:info.joinIdxInPre(c), c), dataPreProjSmooth(2, 1:info.joinIdxInPre(c), c), dataPreProjSmooth(3, 1:info.joinIdxInPre(c), c), '-', 'LineWidth', 2, 'Color', [1 0.5 0.5]);
        plot3(dataPostProjSmooth(1,info.nextIdxInPost(c):end, c), dataPostProjSmooth(2, info.nextIdxInPost(c):end, c), dataPostProjSmooth(3, info.nextIdxInPost(c):end, c), '-', 'LineWidth', 2, 'Color', [0.5 1 0.5]);

        % plot splice results
        dataSplicedProj = TensorUtils.linearCombinationAlongDimension(dataSpliced - mu, 1, coeff', 'replaceNaNWithZero', true);
        plot3(dataSplicedProj(1, :, c), dataSplicedProj(2, :, c), dataSplicedProj(3, :, c), 'b');

        xlabel('pc1'); ylabel('pc2'); zlabel('pc3');
        tv = ThreeVector(gca);
        set(findall(gca, 'Type', 'line'), 'Clipping', 'off');
        
        figure();
        ptstack(2, 1, cat(2, dataPreProj, dataPostProj))
        xline(nPre);
        
        figure();
        ptstack(2, 1, dataSplicedProj)
        xline(info.joinIdxInPre(c));
    end

    dataSpliced = ipermute(dataSpliced, [basisDim, spliceDim, dimConditions]);

end

function in = smoothTraces(in)
    % in is Bases x Time x conditions
    for b = 1:size(in, 1)
        for c = 1:size(in, 3)
            in(b, :, c) = smooth(in(b, :, c), 9, 'sgolay');
        end
    end
end
