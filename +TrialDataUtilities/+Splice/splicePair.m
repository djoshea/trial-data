function [dataSpliced, info] = splicePair(dataPre, dataPost, varargin)
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

    nPre = size(dataPre, 2);
    nPost = size(dataPost, 2);
    nBases = size(dataPre, 1);
    nTraj = prod(TensorUtils.sizeOtherDims(dataPre, [1 2]));
    
    p = inputParser();
    
    % constrain the number of timepoints that could be overlapping between
    % the two trajectories
    p.addParameter('commonOverlapAcrossTrajectories', true, @islogical);
    
    % common timepoint where the splicing occurs, true requires
    % commonOverlapAcrossTrajectories to be true also
    p.addParameter('commonJoinAcrossTrajectories', false, @islogical);
    p.addParameter('minOverlap', 0, @isscalar); 
    p.addParameter('maxOverlap', Inf, @isscalar);
    
    % constrain the search window for the splice break point
    p.addParameter('joinAfterIndexPre', 1, @isnumeric);
    p.addParameter('joinBeforeIndexPost', nPost, @isnumeric); % indices to search into the post trajectory to find the splice point
    
    % for spline interpolation
    p.addParameter('interpolateMethod', 'spline', @ischar);
    p.addParameter('interpIgnoreWindow', 15, @isscalar); % ignore these last points around the trjaectories when fitting the spline
    p.addParameter('interpFitWindow', 30, @isscalar); % include these last points as waypoints for the splines, should be bigger than splineIgnore
    
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    % first, we compute the best amount of temporal overlap of the edges of
    % the trajectories
    nTimepointsOverlap = TrialDataUtilities.Splice.computeBestOverlap(dataPre, dataPost, ...
        p.Results.minOverlap, p.Results.maxOverlap, 'commonAcrossTrajectories', p.Results.commonOverlapAcrossTrajectories, 'showPlot', false);
    info.nTimepointsOverlap = nTimepointsOverlap;
    
    % then we assume that overlap and use the best splice point on a per
    % trajectory basis (i.e. for each traj along dims 3)
    [info.joinIdxInPre, info.nextIdxInPost, dataCat] = TrialDataUtilities.Splice.computeBestJoinPoint(dataPre, dataPost, nTimepointsOverlap, ...
        'joinAfterIndexPre', p.Results.joinAfterIndexPre, 'joinBeforeIndexPost', p.Results.joinBeforeIndexPost, ...
        'commonJoinAcrossTrajectories', p.Results.commonJoinAcrossTrajectories && p.Results.commonOverlapAcrossTrajectories);
    
    % build info matrices describing where each point of the concatenated
    % timeseries pulls its data from for pre and post
    T = size(dataCat, 2);
    C = numel(info.joinIdxInPre);
    szCat = size(dataCat);
    [info.idxFromPre, info.idxFromPost] = deal(nan([T, szCat(3:end)]));
    for c = 1:C
        info.idxFromPre(1:info.joinIdxInPre(c), c) = 1:info.joinIdxInPre(c);
        info.idxFromPost(info.joinIdxInPre(c)+1:end, c) = info.nextIdxInPost(c) : nPost;
    end
    
    % do smooth interpolation at the join
    if strcmp(p.Results.interpolateMethod, 'spline')
        dataSpliced = TrialDataUtilities.Splice.interpolateSpline(dataCat, info.joinIdxInPre+1, ...
            p.Results.interpFitWindow, p.Results.interpIgnoreWindow);

    elseif ismember(p.Results.interpolateMethod, {'linear', 'nearest', 'next', 'previous', 'cspline', 'pchip', 'pchip', 'cubic', 'v5cubic'})
        % do linear interpolation
        dataSpliced = TrialDataUtilities.Splice.interpolateLinear(dataCat, info.joinIdxInPre+1, ...
            p.Results.interpolateMethod, p.Results.interpFitWindow, p.Results.interpIgnoreWindow);
    else
        
        dataSpliced = dataCat;
    end
end
