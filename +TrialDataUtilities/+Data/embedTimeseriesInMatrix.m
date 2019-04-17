function [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeCell, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
% N is trial count, G and C are channel counts (depending on the
% configuration of dataCell, you can have multiple channels inside each element or
% have dataCell be a cell matrix where each column (G) is a channel (C==1)
% or group of channels (C > 1)). In the latter case, the timestamps are
% shared across the C channels. The G channels in the former case need not
% share a time vector.
%
%   dataCell is N x G cell of T_n x C matrices (where C can be vector if
%      contents are tensors
%   timeCell is N x G cell of T_n x 1 vectors
%
% mat will be N x T x C(:) x G, where T is the length of the time vector tvec to
% which all data have been resampled to the same sampling frequency.
% 
% Combines a cell of N timeseries into a matrix with size N x T with a
% common time vector tvec. dataCell is a N x 1 cell of numeric vectors.
% timeData is either a N x 1 cell of corresponding time vectors.
%
% tvec will be chosen to span the minimum and maximum timestamps across
% all elements of dataCell with the specified timeDelta, i.e.
% min(tMin):timeDelta:max(tMax).
%
% % Optional args: 
%    tvec: if specified, will serve as the time vector to interpolate 
%       data onto. if not specified, will be computed automatically.
%
%    timeDelta: scalar indicating the time delta between successive
%       timestamps. defaults to the minimum delta automatically determined
%       across all G channels in the outer cell
%
%    timeReference: scalar indicating the reference point for generating
%       time vectors. Time vectors will be generated using the
%       tMin:timeDelta:tMax formula, but tMin and tMax will be adjusted so
%       that timeReference would be an exact integer multiple of timeDelta
%       away. default = 0.
% 
%    fixDuplicateTimes: remove timestamps which are duplicates of previous
%       timestamps, which can result from rounding and breaks the
%       monotonicity required by interpolation
%
%    assumeUniformSampling: if true, we will only shift adjacent trials to
%       make the time vectors line up in the matrix. if false, interpolation
%       will be used to ensure the signals from each trial are sampled at the
%       correct, evenly spaced times.
%    
%    interpolateMethod: string. See interp1 help for description of 
%       interpolation methods. default = 'pchip'.
%
%    minTrials : if > 0, will shrink the edges of mat and tvec to timepoints where
%    at least this many trials have data
%
%    minTrialFraction : if > 0, will shrink the edges of mat and tvec to timepoints where
%    at least this fraction of trials have data. the fraction will be taken
%    of valid trials only.
%
%    trialValid : mask of which trials to consider when computing minTrials
%    and minTrialFraction

    p = inputParser();
    p.addRequired('dataCell', @(x) iscell(x));
    p.addRequired('timeCell', @(x) iscell(x));
    p.addParameter('assumeUniformSampling', false, @islogical); % allows us to skip the interp/resample step if the data are known to be uniform
%     p.addParameter('tvec', [], @(x) isempty(x) || isvector(x));
    p.addParameter('interpolateMethod', 'linear', @ischar);
    
    p.addParameter('fixNonmonotonicTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('resampleMethod', 'filter', @ischar); % valid modes are filter, average, repeat , interp   
    
    p.addParameter('showProgress', true, @islogical);
    p.addParameter('minTrials', 0, @isscalar);
    p.addParameter('minTrialFraction', 0, @isscalar);
    p.addParameter('trialValid', [], @islogical);
    
    % these are used as a secondary guard to truncate data within tMin :
    % tMax, when the input data includes padded edges to facilitate
    % resampling
    p.addParameter('tMinExcludingPadding', [], @ismatrix);
    p.addParameter('tMaxExcludingPadding', [], @ismatrix);
    
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x)); % specify manually to save time if known, otherwise will be inferred
    p.addParameter('ignoreNaNSamples', false, @islogical); % ignore NaN data samples when inferring origDelta (time skips among successive non-nan samples)
    
    p.PartialMatching = false;
    p.parse(dataCell, timeCell, varargin{:});
    
    % check sizes match
    % okay to have one empty and the other not, simply ignore
    szData = cellfun(@(x) size(x, 1), dataCell);
    szTime = cellfun(@numel, timeCell);
    
    if isempty(szData) || isempty(szTime)
        mat = nan(size(dataCell, 1), 0, 0);
        tvec = zeros(0, 1);
        return;
    end
    
    empty = all(szData == 0, 2) | all(szTime == 0, 2);
    if all(empty)
        mat = nan(size(dataCell, 1), 0, size(dataCell, 2));
        tvec = zeros(0, 1);
        return;
    end
    
    szData(empty, :) = 0;
    szTime(empty, :) = 0;
    dataCell(empty, :) = {[]};
    timeCell(empty, :) = {[]};
    assert(all(TensorUtils.flatten(bsxfun(@eq, szData, szTime))), 'Sizes of dataCell and timeCell contents must match');
    
    % check column counts match (after converting to 2d matrix)
    cData = cellfun(@(x) size(x(:, :), 2), dataCell);
    uniqueColCounts = unique(cData(~empty));
    assert(numel(uniqueColCounts) == 1, 'All entries must have ');
    C = uniqueColCounts(1);
    
    if isempty(dataCell)
        mat = [];
        tvec = zeros(0, 1);
        return;
    end
    
    % figure out the shape of dataCell contents so we can reshape later
    first = find(~empty, 1, 'first');
    Cvec = size(dataCell{first});
    Cvec = Cvec(2:end);
    
    N = size(dataCell, 1);
    trialValid = p.Results.trialValid;
    if isempty(trialValid)
        trialValid = truevec(N);
    end
    
    % fix duplicate timestamps
    if p.Results.fixNonmonotonicTimes
        [timeCell, dataCell] = TrialDataUtilities.Data.fixNonmonotonicTimeseries(timeCell, dataCell);
    end

    if size(timeCell, 2) == 1 && size(dataCell, 2) > 1
        timeCell = repmat(timeCell, 1, size(dataCell, 2));
    end
    assert(size(timeCell, 2) == size(dataCell, 2), 'Column counts of data cell and time cell must match');
    
    timeDelta = p.Results.timeDelta;
    
    
    if isempty(p.Results.origDelta)
        origDelta = TrialDataUtilities.Data.inferTimeDeltaFromSampleTimes(timeCell, dataCell, 'ignoreNaNSamples', p.Results.ignoreNaNSamples);
    else
        origDelta = p.Results.origDelta;
    end
    [tvec, tMin, tMax] = TrialDataUtilities.Data.inferCommonTimeVectorForTimeseriesData(timeCell, dataCell, ...
        'timeDelta', timeDelta, 'timeReference', p.Results.timeReference, ...
        'binAlignmentMode', p.Results.binAlignmentMode, ...
        'origDelta', origDelta, ...
        'fixNonmonotonicTimes', false, ... % consider already fixed
        'tMinExcludingPadding', p.Results.tMinExcludingPadding, 'tMaxExcludingPadding', p.Results.tMaxExcludingPadding);
    
    timeCell = TrialDataUtilities.Data.removeSmallTimeErrors(timeCell, origDelta, p.Results.timeReference);
    
    % tMin and tMax are now vectors of the start and stop times of each
    % trial
    tvec = makecol(tvec);
    
    tMinGlobal = nanmin(tvec);
    timeDelta = p.Results.timeDelta;
    if isempty(timeDelta)
        if numel(tvec) == 1
            % special case with one sample
            timeDelta = 1;
        else
            timeDelta = nanmin(origDelta);
        end
    end
    
    % clean up small inconsistencies due to floating point
    tMin = TrialDataUtilities.Data.removeSmallTimeErrors(tMin, timeDelta, p.Results.timeReference);
    tMax = TrialDataUtilities.Data.removeSmallTimeErrors(tMax, timeDelta, p.Results.timeReference);
    
    % build the data matrix by inserting the interpolated segment of each timeseries
    % in the appropriate location in each row, keeping the non-spanned timepoints as NaN
    T = numel(tvec);
    
    G = size(dataCell, 2);
    
    dclass = TrialDataUtilities.Data.getCellElementClass(dataCell);
    if ~ismember(dclass, {'single', 'double'})
        dclass = 'single';
    end
    mat = nan([N, T, C, G], dclass); % we'll reshape this later, C is channels per matrix of dataCell, G is over columns of dataCell    
    
    indPutStart = TrialDataUtilities.Stats.floortol((tMin - tMinGlobal) / timeDelta, timeDelta/1000) + 1;
    indPutStop  = TrialDataUtilities.Stats.floortol((tMax - tMinGlobal) / timeDelta, timeDelta/1000) + 1;
    
    % if specified, we can skip the resampling step which makes this quick
    uniformSampling = p.Results.assumeUniformSampling && all(timeDelta == origDelta);
    
    err = timeDelta / 1000;
    for i = 1:N % loop over trials
        if ~trialValid(i), continue; end
%         if mod(i, 100) && p.Results.showProgress, prog.update(i); end

        for g = 1:G % loop over second dim of cell, which is over channels that need not share common time 

            if ~isnan(tMin(i,g)) && ~isnan(tMax(i,g))    
                if uniformSampling
                    % just figure out where to insert this into the matrix,
                    % no need for resampling
                    indTake = timeCell{i, g} + err >= tvec(indPutStart(i, g)) & timeCell{i,g} <= tvec(indPutStop(i, g)) + err;
                    mat(i, indPutStart(i,g):indPutStop(i,g), :, g) = dataCell{i,g}(indTake, :);
                else
                    mask = ~all(isnan(dataCell{i, g}), 2);

                    % use resampleTensorInTime to support multiple resampling
                    % methods
                    vals = TrialDataUtilities.Data.resampleTensorInTime(dataCell{i,g}(mask, :), 1, timeCell{i, g}(mask), ...
                        'origDelta', origDelta(g), 'timeDelta', timeDelta, 'interpolateMethod', p.Results.interpolateMethod, ...
                        'binAlignmentMode', p.Results.binAlignmentMode, 'resampleMethod', p.Results.resampleMethod, ...
                        'origDelta', origDelta(g), 'tMin', tMin(i,g), 'tMax', tMax(i,g));
                    
                    mat(i, indPutStart(i,g):indPutStop(i,g), :, g) = vals;
                end
            end
        end
    end
    
    % pare down time points from edges with insufficient trial counts
    if p.Results.minTrials > 0 || p.Results.minTrialFraction > 0
        nTrialsOverTime = sum(all(~isnan(mat(trialValid, :, :)), 3), 1);
        minTrials = max(p.Results.minTrials, nnz(trialValid) * p.Results.minTrialFraction);
        
        tMask= falsevec(T);
        tMask(find(nTrialsOverTime >= minTrials, 1, 'first') : find(nTrialsOverTime >= minTrials, 1, 'last')) = true;
    else
        tMask = truevec(T);
    end
    
    if C == 1
        mat = TensorUtils.squeezeDims(mat(:, tMask, :, :), 3);
    else
        mat = reshape(mat(:, tMask, :, :), [N nnz(tMask) Cvec G]);
    end
    tvec = tvec(tMask);
end

