function [dataSpliced, spliceBreaks] = spliceTrajectories(dataCell, varargin)
% Splices together high dimensional state trajectories along the time
% dimension dim 2. 
%
% First, we determine the amount of temporal overlap
% between the two trajectories, which works better if there is sufficient
% real overlap between the two. This can be constrained using minOverlap
% and maxOverlap. This is done once for all trajectories (each position along dims 3 and
% beyond).
%
% Next, we find the best "splice point" at which to stop taking data from
% dataPre and start taking data from dataPost. The range of possible splice 
% points here can be constrained using searchPre and searchPost, which both
% default to Inf (any splice point is valid).
%
% Lastly, optionally, we chop off some data to either side and replace it
% with a spline interpolation.
%
% dataCell is nPeriods x 1 cell , which each cell containing a time window
% all of which are to be spliced together
% each cell inside is size
% D (typically 2 or 3 dimensions) x T x ...
% splicing will occur for each pair of trajectories
%
% alternatively, dataCell is size D x T x ... and splicing is done at the
% indices along dimension 2 given by newSpliceStart. For example, to splice
% times 1:3 to 4:6 to 7:10, specify newSpliceStart = [4 7]
%
% numTimepointsDrop is a vector indicating how many timepoints should be
% dropped at each splice point. 

p = inputParser();
p.addParameter('newSpliceStart', [], @isvector); % specify when dataCell is a matrix
p.addParameter('basisMask', [], @isvector);
p.KeepUnmatched = true;
p.parse(varargin{:});


if iscell(dataCell)
    assert(isvector(dataCell));
    dataCell = makecol(dataCell);
    nSplice = numel(dataCell);

    tPerSplice = cellfun(@(x) size(x, 2), dataCell);
    
    spliceStart = [1; cumsum(tPerSplice(1:end-1))+1];
    spliceStop = cumsum(tPerSplice);
    dataCat = cat(2, dataCell{:});
else
    dataCat = dataCell;
    assert(isnumeric(dataPieces));
    spliceStart = [1; makecol(p.Results.newSpliceStart)];
    spliceStop = [spliceStart(2:end)-1; size(dataCat, 2)];
end

spliceSizes = spliceStop - spliceStart + 1;

basisMask = p.Results.basisMask;
if isempty(basisMask)
    basisMask = truevec(size(dataCat, 1));
else
    basisMask = TensorUtils.vectorIndicesToMask(basisMask, size(dataCat, 1));
end

% first, we compute the best amount of temporal overlap of the edges of
% the trajectories
nTimepointsOverlap = computeBestOverlap(dataPre, dataPost, p.Results.minOverlap, p.Results.maxOverlap);

% figure out how big the resulting pieces will be
insertSizes = spliceSizes(1:end-1) + spliceSizes(2:end) - nTimepointsOverlap;
insertStart = [1; cumsum(insertSizes(1:end-1))+1];
insertStop = cumsum(insertSizes);

% start by raw splicing data
sz = size(dataCat);
nTraj = prod(sz(3:end));

dataSpliced = nan(sz(1), sum(insertSizes), sz(3:end), 'like', dataCat);

spliceBreaks = nan(nSplice - 1, nTraj);
prog = ProgressBar(nTraj * (nSplice-1));
for i = 1:nTraj
    for s = 2:nSplice
        prog.increment();
        [dataSpliced(basisMask, insertStart(s-1):insertStop(s-1), i), spliceBreaks(s-1, i)] = TrialDataUtilities.Data.spliceSingleTrajectoryPair(...
            dataCat(basisMask, spliceStart(s-1):spliceStop(s-1), i), ...
            dataCat(basisMask, spliceStart(s):spliceStop(s), i), ...
            'numTimepointsDrop', numTimepointsDrop(s-1), ...
            p.Unmatched);
    end
end
prog.finish();


