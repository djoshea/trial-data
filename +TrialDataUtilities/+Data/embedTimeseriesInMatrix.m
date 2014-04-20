function [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeCell, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
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
%    interpolate: boolean. if false, it will be assumed that the timeStamps
%       within each dataCell have consistent spacing and each timestamp lies
%       an integer multiple away from timeReference. If true, data samples
%       will be interpolated to timestamps which satisfy this requirement.
%       default = true.
%
%    interpolateMethod: string. See interp1 help for description of 
%       interpolation methods. default = 'linear'.

    p = inputParser();
    p.addRequired('dataCell', @(x) iscell(x));
    p.addRequired('timeCell', @(x) iscell(x));
    p.addParamValue('tvec', [], @(x) isempty(x) || isvector(x));
    p.addParamValue('interpolateMethod', 'linear', @ischar);
    p.addParamValue('fixDuplicateTimes', true, @(x) islogical(x) && isscalar(x));
    p.KeepUnmatched = true;
    p.parse(dataCell, timeCell, varargin{:});
    
    % check sizes match
    % okay to have one empty and the other not, simply ignore
    szData = cellfun(@numel, dataCell);
    szTime = cellfun(@numel, timeCell);
    empty = szData == 0 | szTime == 0;
    szData(empty) = 0;
    szTime(empty) = 0;
    dataCell(empty) = {[]};
    timeCell(empty) = {[]};
    assert(all(szData(:) == szTime(:)), 'Sizes of dataCell and timeCell contents must match');
    
    % fix duplicate timestamps
    if p.Results.fixDuplicateTimes
        [dataCell, timeCell] = cellfun(@fix, dataCell, timeCell, 'UniformOutput', false);
    end

    if isempty(p.Results.tvec)
        % auto-compute appropriate time vector
        [tvec, tMin, tMax] = TrialDataUtilities.Data.inferCommonTimeVectorForTimeseriesData(timeCell, p.Unmatched);
    else
        tvec = p.Results.tvec;
        % compute the global min / max timestamps or each trial
        [tMinRaw, tMaxRaw] = cellfun(@minmax, timeCell);
        
        % then shift these to lie within tvec, without overwriting nans
        tMin = max(min(tvec), tMinRaw);
        tMin(isnan(tMinRaw)) = NaN;
        
        tMax = min(max(tvec), tMaxRaw);
        tMax(isnan(tMaxRaw)) = NaN;
    end
    
    tMinGlobal = nanmin(tvec);
    timeDelta = inferTimeDelta(tvec);
    
    % build the data matrix by inserting the interpolated segment of each timeseries
    % in the appropriate location in each row, keeping the non-spanned timepoints as NaN
    T = numel(tvec);
    N = size(dataCell, 1);
    C = size(dataCell, 2);
    mat = nan([N, T, C]);
	indStart = floor(((tMin - tMinGlobal) / timeDelta) + 1);
    indStop  = floor(((tMax - tMinGlobal) / timeDelta) + 1);
    for c = 1:C
        for i = 1:N
            if ~isnan(indStart(i,c)) && ~isnan(indStop(i,c))
                mat(i, indStart(i,c):indStop(i,c), c) = interp1(timeCell{i, c}, dataCell{i, c}, ...
                    tvec(indStart(i,c):indStop(i,c)), p.Results.interpolateMethod, 'extrap');
            end
        end
    end
end

function timeDelta = inferTimeDelta(tvec)
     timeDelta = nanmedian(diff(tvec));
end

function [d, t] = fix(d, t)
    if isempty(t) || isempty(d)
        d = [];
        t = [];
        return;
    end
    
    diffT = diff(t);
    stuck = find(diffT(1:end-1) == 0 & diffT(2:end) == 2);
    t(stuck+1) = t(stuck+1) + 1;
    skip = find(diffT(1:end-1) == 2 & diffT(2:end) == 0);
    t(skip+1) = t(skip+1) - 1;

    tMask = diff(t)>0;
    t = t(tMask);
    d = d(tMask);
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