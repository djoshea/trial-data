function [tvec, tMinCell, tMaxCell] = inferCommonTimeVectorForTimeseriesData(timeCell, dataCell, varargin)
% Determines a time vector that can be used to interpolate data from
% one or more cells of timeseries data. This vector will have evenly spaced
% timepoints that span the minimum and maximum of the widest (if 'tight' is
% false) or the narrowest (if tight is 'true', default) time window across
% all of the data sets, e.g. timeCellCell{i} and dataCellCell{i}. 
%
% timeCell is a nTrials x nChannels cell matrix of time vectors
%
%    timeDelta: scalar indicating the time delta between successive
%       timestamps. if not provided, will be inferred from the minimum 
%       spacing of each timeCell, though you should provide this for
%       consistent results.
%
%    timeReference: scalar indicating the reference point for generating
%       time vectors. Time vectors will be generated using the
%       tMin:timeDelta:tMax formula, but tMin and tMax will be adjusted so
%       that timeReference would be an exact integer multiple of timeDelta
%       away from tMin or tMax, that is, so that timeReference appears
%       somewhere in the tvec. default = 0.
%

    p = inputParser();
    p.addRequired('timeCell', @(x) iscell(x));
    p.addParamValue('timeDelta', [], @isscalar);
    p.addParamValue('timeReference', 0, @isscalar);
    p.addParamValue('fixDuplicateTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParamValue('interpolate', true, @(x) islogical(x) && isscalar(x));
   % p.addParamValue('interpolateMethod', 'linear', @ischar);
    p.parse(timeCell, varargin{:});

    timeDelta = p.Results.timeDelta;
    timeReference = p.Results.timeReference;
    fixDuplicateTimes = p.Results.fixDuplicateTimes;
    interpolate = p.Results.interpolate;
    %interpolateMethod = p.Results.interpolateMethod;
   
     % compute the global min / max timestamps or each trial
     [tMinRaw, tMaxRaw] = TrialDataUtilities.Data.getValidTimeExtents(timeCell, dataCell);
%      [tMinRaw, tMaxRaw] = cellfun(@minmax, timeCell);

     if isempty(timeDelta)
         % auto-compute time delta
         warning('Auto-computing time-delta from timeseries. Specify time delta for consistent results');
         if fixDuplicateTimes
            timeCell = cellfun(@fixDup, timeCell, 'UniformOutput', false);
        end
         deltaCell = cellfun(@(x) nanmedian(diff(x)), timeCell);
         % compute the median delta for each channel
         timeDelta = nanmedian(deltaCell(:), 1);
         % and use the minimum spacing
         timeDelta = nanmin(timeDelta);
     end
     
    % auto-compute appropriate time vector

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

        tol = 1e-9;
        if  (any(abs(tMin(:) - tMinRaw(:)) > tol) || any(abs(tMax(:) - tMaxRaw(:)) > tol))
            error('Timestamps do not align with timeReference. Set ''interpolate'' to true');
        end
    end

    % build the global time vector
    tMinGlobal = nanmin(tMin);
    tMaxGlobal = nanmax(tMax);
    tvec = makecol(tMinGlobal:timeDelta:tMaxGlobal);
    tMinCell = tMin;
    tMaxCell = tMax;
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

function t = fixDup(t)
    if isempty(t)
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
end
