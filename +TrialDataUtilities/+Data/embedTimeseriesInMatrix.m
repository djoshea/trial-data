function [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeCell, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
% N is trial count, G and C are channel counts (depending on the
% configuration of dataCell, you can have multiple channels inside each element or
% have dataCell be a cell matrix where each column (G) is a channel (C==1) or group of channels (C > 1))
%   dataCell is N x G cell of T_n x C matrices
%   timeCell is N x G cell of T_n vectors
%
% mat will be N x T x (C*G), where T is the length of the time vector tvec to
% which all data have been interpolated
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
%       timestamps. default = 1.
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

    p = inputParser();
    p.addRequired('dataCell', @(x) iscell(x));
    p.addRequired('timeCell', @(x) iscell(x));
    p.addParameter('assumeUniformSampling', false, @islogical);
    p.addParameter('tvec', [], @(x) isempty(x) || isvector(x));
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('fixDuplicateTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('showProgress', true, @islogical);
%     p.addParameter('sparse', false, @islogical);
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
    
    % check column counts match
    cData = cellfun(@(x) size(x, 2), dataCell);
    uniqueColCounts = unique(cData(~empty));
    assert(numel(uniqueColCounts) == 1, 'All entries must have ');
    C = uniqueColCounts(1);
    
    if isempty(dataCell)
        mat = [];
        tvec = zeros(0, 1);
        return;
    end
    
    % fix duplicate timestamps
    if p.Results.fixDuplicateTimes
        [timeCell, dataCell] = TrialDataUtilities.Data.fixNonmonotonicTimeseries(timeCell, dataCell);
    end

    if size(timeCell, 2) == 1 && size(dataCell, 2) > 1
        timeCell = repmat(timeCell, 1, size(dataCell, 2));
    end
    assert(size(timeCell, 2) == size(dataCell, 2), 'Column counts of data cell and time cell must match');
    
    if isempty(p.Results.tvec)
        % auto-compute appropriate time vector
        [tvec, tMin, tMax] = TrialDataUtilities.Data.inferCommonTimeVectorForTimeseriesData(timeCell, dataCell, ...
            'timeDelta', p.Results.timeDelta, 'timeReference', p.Results.timeReference, ...
            'interpolate', ~p.Results.assumeUniformSampling);
    else
        tvec = p.Results.tvec;
        % compute the global min / max timestamps or each trial
        [tMinRaw, tMaxRaw] = TrialDataUtilities.Data.getValidTimeExtents(timeCell, dataCell);
%         [tMinRaw, tMaxRaw] = cellfun(@minmax, timeCell);
        
        % then shift these to lie within tvec, without overwriting nans
        tMin = max(min(tvec), tMinRaw);
        tMin(isnan(tMinRaw)) = NaN;
        
        tMax = min(max(tvec), tMaxRaw);
        tMax(isnan(tMaxRaw)) = NaN;
    end
    
    % tMin and tMax are now vectors of the start and stop times of each
    % trial
    tvec = makecol(tvec);
    
    tMinGlobal = nanmin(tvec);
    if numel(tvec) == 1
        % special case with one sample
        timeDelta = 1;
    else
        timeDelta = inferTimeDelta(tvec);
    end
    
    % build the data matrix by inserting the interpolated segment of each timeseries
    % in the appropriate location in each row, keeping the non-spanned timepoints as NaN
    T = numel(tvec);
    N = size(dataCell, 1);
    G = size(dataCell, 2);
    
    mat = nan([N, T, C, G]); % we'll reshape this later

    indStart = floor(((tMin - tMinGlobal) / timeDelta) + 1);
    indStop  = floor(((tMax - tMinGlobal) / timeDelta) + 1);
    
    
    
    if p.Results.showProgress
        prog = ProgressBar(N, 'Embedding data over trials into common time vector');
    end
    for i = 1:N
        if p.Results.showProgress, prog.update(i); end
        for g = 1:G
            if ~isnan(indStart(i,g)) && ~isnan(indStop(i,g))
                if numel(indStart(i,g):indStop(i,g)) > 1
                    mask = ~all(isnan(dataCell{i, g}), 2);
                    if p.Results.assumeUniformSampling
                        % in this case, we just make sure the timepoint
                        % closest to zero ends up in gthe right place
                        [thisTimeRef, thisIndRef] = min(abs(timeCell{i, g}));
                        [~, tvecIndRef] = min(abs(thisTimeRef - tvec));

                        locationInMat = (1:numel(timeCell{i, g})) + tvecIndRef-thisIndRef;
                        locationInMat = locationInMat(mask);

                        mask = mask(locationInMat >= 1 & locationInMat <= T);
                        locationInMat = locationInMat(mask);

                        % data cell is T x C
                        mat(i, locationInMat, :, g) = dataCell{i, g}(mask, :);
                    else
                        % don't assume uniform sampling, just interpolate
                        % to the right time vector
                        mat(i, indStart(i,g):indStop(i,g), :, g) = interp1(double(timeCell{i, g}(mask)), dataCell{i, g}(mask, :), ...
                            tvec(indStart(i,g):indStop(i,g)), p.Results.interpolateMethod, 'extrap');
                    end
                else
                    mat(i, indStart(i,g), :, g) = dataCell{i, g};
                end
            end
        end
    end
    if p.Results.showProgress
        prog.finish();
    end
    
    mat = reshape(mat, [N T C*G]);
end

function timeDelta = inferTimeDelta(tvec)
     timeDelta = nanmedian(diff(tvec));
end

function [mn, mx] = minmax(x)
    if isempty(x)
        mn = NaN;
        mx = NaN;
    else
        mn = nanmin(x);
        mx = nanmax(x);
    end
end