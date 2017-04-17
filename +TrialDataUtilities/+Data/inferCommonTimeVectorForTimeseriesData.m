function [tvec, tMinCell, tMaxCell, origDelta, indMin, indMax] = inferCommonTimeVectorForTimeseriesData(timeCell, dataCell, varargin)
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
% origDelta is nChannels x 1 indicator of the time between samples for each
% channel
%

    p = inputParser();
    p.addRequired('timeCell', @(x) iscell(x));
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('fixDuplicateTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('interpolate', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    
    % these are used as a secondary guard to truncate data within tMin :
    % tMax, when the input data includes padded edges to facilitate
    % resampling
    p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
    p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
    
   % p.addParamValue('interpolateMethod', 'linear', @ischar);
    p.parse(timeCell, varargin{:});

    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    fixDuplicateTimes = p.Results.fixDuplicateTimes;
    binAlignmentMode = p.Results.binAlignmentMode;
    interpolate = p.Results.interpolate;
    %interpolateMethod = p.Results.interpolateMethod;
   
     % compute the global min / max timestamps or each trial
    if fixDuplicateTimes
        [timeCell, dataCell] = cellfun(@fixDup, timeCell, dataCell, 'UniformOutput', false);
    end
    [tMinRaw, tMaxRaw, origDelta, indMin, indMax] = TrialDataUtilities.Data.getValidTimeExtents(timeCell, dataCell);

    nTimes = cellfun(@numel, timeCell);

    if isempty(timeDelta)
        if all(nTimes == 1)
            timeDelta = 1;
        else
            timeDelta = min(origDelta);
        end
     end

    % expand the global min / max timestamps to align with timeReference
    if interpolate
        %ceilfix = @(x)ceil(abs(x)).*sign(x);
        % these used to be round, changed to match start / stop in Trial
        % Data
        
        % switching from floor here to ceil
        % the logic is this: a sample at time t represents an average of
        % data in interval t-origDelta:t, whereas origDelta is the original sampling rate.
        % Therefore, we can only include a new time k if we have data from k-timeDelta:k, 
        % where timeDelta is the new sampling rate.
        % so we have data from t(1)-origDelta : t(end)
        %
        % e.g. t(1) == 0, origDelta = 1, timeDelta = 10 --> tMin = 10
        % e.g. t(1) == 0, origDelta = 1, timeDelta = 1 --> tMin = 0

        % should work due to vectorization and auto expanding (R2017 only)
%         [tMin, tMax] = binAlignmentMode.getTimeLimitsForRebinning(tMinRaw, tMaxRaw, origDelta, timeDelta, timeReference);
        
        [tMin, tMax] = deal(nan(size(tMinRaw)));
        for g = 1:size(dataCell, 2)
            [tMin(:, g), tMax(:, g)] = binAlignmentMode.getTimeLimitsForRebinning(tMinRaw(:, g), tMaxRaw(:, g), origDelta(g), timeDelta, timeReference);
        end
    else
        tMin = tMinRaw;
        tMax = tMaxRaw;
    end
    
    if ~isempty(p.Results.tMinExcludingPadding)
        tMin = max(tMin, p.Results.tMinExcludingPadding, 'includenan');
    end
    if ~isempty(p.Results.tMaxExcludingPadding)
        tMax = min(tMax, p.Results.tMaxExcludingPadding, 'includenan');
    end

    % build the global time vector
    tMinGlobal = nanmin(tMin(:));
    tMaxGlobal = nanmax(tMax(:));
    tvec = makecol(tMinGlobal:timeDelta:tMaxGlobal);
    tMinCell = tMin;
    tMaxCell = tMax;
end

function [t, d] = fixDup(t, d)
    if isempty(t)
        t = [];
        return;
    end
    
    diffT = diff(t);
    stuck = find(diffT(1:end-1) == 0 & diffT(2:end) == 2);
    t(stuck+1) = t(stuck+1) + 1;
    skip = find(diffT(1:end-1) == 2 & diffT(2:end) == 0);
    t(skip+1) = t(skip+1) - 1;

    tMask = [true; diff(t)>0];
    t = t(tMask);
    d = d(tMask, :);
end
