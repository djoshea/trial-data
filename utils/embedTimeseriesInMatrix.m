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
    p.addRequired('dataCell', @(x) iscell(x) && isvector(x));
    p.addRequired('timeCell', @(x) iscell(x) && isvector(x));
    p.addParamValue('timeDelta', 1, @isscalar);
    p.addParamValue('timeReference', 0, @isscalar);
    p.addParamValue('fixDuplicateTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParamValue('interpolate', true, @(x) islogical(x) && isscalar(x));
    p.addParamValue('interpolateMethod', 'linear', @ischar);
    p.parse(dataCell, timeCell, varargin{:});

    timeDelta = p.Results.timeDelta;
    timeReference = p.Results.timeReference;
    fixDuplicateTimes = p.Results.fixDuplicateTimes;
    interpolate = p.Results.interpolate;
    interpolateMethod = p.Results.interpolateMethod;
    
    szData = cellfun(@numel, dataCell);
    szTime = cellfun(@numel, timeCell);
    
    assert(all(szData == szTime), 'Sizes of dataCell and timeCell contents must match');
    
    if fixDuplicateTimes
        [dataCell, timeCell] = cellfun(@fix, dataCell, timeCell, 'UniformOutput', false);
    end
    
    % compute the global min / max timestamps
    [tMinRaw, tMaxRaw] = cellfun(@minmax, timeCell);

    % expand the global min / max timestamps to align with timeReference
    if interpolate
        %ceilfix = @(x)ceil(abs(x)).*sign(x);
        % these used to be round, changed to match start / stop in Trial
        % Data
        tMin = timeReference + floor((tMinRaw - timeReference) / timeDelta) * timeDelta;
        tMax = timeReference + ceil((tMaxRaw - timeReference) / timeDelta) * timeDelta;
    else
        tMin = tMinRaw;
        tMax = tMaxRaw;
    end
    
    if ~interpolate
        tol = 1e-9;
        if  (any(abs(tMin - tMinRaw) > tol) || any(abs(tMax - tMaxRaw) > tol))
            error('Timestamps do not align with timeReference. Set ''interpolate'' to true');
        end
    end
    
    % build the global time vector
    tMinGlobal = nanmin(tMin);
    tMaxGlobal = nanmax(tMax);
    tvec = tMinGlobal:timeDelta:tMaxGlobal;
    T = numel(tvec);
    N = numel(dataCell);
    
    % build the data matrix by inserting each timeseries in the appropriate
    % location in each row
    mat = nan(N, T);
    
	indStart = floor(((tMin - tMinGlobal) / timeDelta) + 1);
    indStop  = floor(((tMax - tMinGlobal) / timeDelta) + 1);
    if interpolate
        for i = 1:N
            if ~isnan(indStart(i)) && ~isnan(indStop(i))
                mat(i, indStart(i):indStop(i)) = interp1(timeCell{i}, dataCell{i}, ...
                    tvec(indStart(i):indStop(i)), interpolateMethod, 'extrap');
            end
        end
    else
        for i = 1:N
            if ~isnan(indStart(i)) && ~isnan(indStop(i))
                mat(i, indStart(i):indStop(i)) = dataCell{i};
            end
        end
    end
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