function dataSpliced = spliceTrajectories(dataCell, varargin)
% find a way to splice N-d timeseries by deleteing points to the left and
% to the right. 
% dataCell is nPeriods x 1 cell , which each cell containing a time window
% all of which are to be spliced together
% each cell inside is size
% D (typically 2 or 3 dimensions) x T x ...
% splicing will occur for each pair of trajectories
%
% alternatively, dataCell is size D x T x ... and splicing is done at the
% indices along dimension 2 given by newSpliceStart. For example, to splice
% times 1:3 to 4:6 to 7:10, specify newSpliceStart = [4 7]

p = inputParser();
p.addParameter('newSpliceStart', [], @isvector); % specify when dataCell is a matrix
p.addParameter('basisMask', [], @isvector);
p.KeepUnmatched = true;
p.parse(varargin{:});

if iscell(dataCell)
    assert(isvector(dataCell));
    dataCell = makecol(dataCell);
    nSplice = numel(dataCell);

    sz = size(dataCell{1});
    tPerSplice = cellfun(@(x) size(x, 2), dataCell);
    sz(2) = sum(tPerSplice);

    spliceStart = [1; cumsum(tPerSplice(1:end-1))+1];
    spliceStop = cumsum(tPerSplice);
    dataSpliced = cat(2, dataCell{:});
else
    dataSpliced = dataCell;
    assert(isnumeric(dataSpliced));
    spliceStart = [1 makecol(p.Results.newSpliceStart)];
    spliceStop = [spliceStart(2:end)-1; size(dataSpliced, 2)];
end

basisMask = p.Results.basisMask;
if isempty(basisMask)
    basisMask = truevec(size(dataSpliced, 1));
else
    basisMask = TensorUtils.vectorIndicesToMask(basisMask, size(dataSpliced, 1));
end

% start by raw splicing data
nTraj = prod(sz(3:end));
    
prog = ProgressBar(nTraj * (nSplice-1));
for i = 1:nTraj
    for s = 2:nSplice
        prog.increment();
        dataSpliced(basisMask, spliceStart(s-1):spliceStop(s), i) = TrialDataUtilities.Data.spliceSingleTrajectoryPair(...
            dataSpliced(basisMask, spliceStart(s-1):spliceStop(s-1), i), ...
            dataSpliced(basisMask, spliceStart(s):spliceStop(s), i), p.Unmatched);
    end
end
prog.finish();


