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
    p.addParameter('fixNonmonotonicTimes', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('interpolate', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    
    % these are used as a secondary guard to truncate data within tMin :
    % tMax, when the input data includes padded edges to facilitate
    % resampling
    p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
    p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
    
    p.addParameter('origDelta', [], @(x) isempty(x) || isvector(x)); % specify manually to save time if known, otherwise will be inferred
    p.addParameter('ignoreNaNSamples', false, @islogical); % ignore NaN data samples when inferring origDelta (time skips among successive non-nan samples)
    
   % p.addParamValue('interpolateMethod', 'linear', @ischar);
    p.parse(timeCell, varargin{:});

    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    binAlignmentMode = p.Results.binAlignmentMode;
    interpolate = p.Results.interpolate;
    %interpolateMethod = p.Results.interpolateMethod;
   
    % compute the global min / max timestamps or each trial
    if p.Results.fixNonmonotonicTimes
        [timeCell, dataCell] = TrialDataUtilities.Data.fixNonmonotonicTimeseries(timeCell, dataCell);
    end
    
    if isempty(p.Results.origDelta)
        origDelta = TrialDataUtilities.Data.inferTimeDeltaFromSampleTimes(timeCell, dataCell, 'ignoreNaNSamples', p.Results.ignoreNaNSamples);
    else
        origDelta = p.Results.origDelta;
    end
    [tMinRaw, tMaxRaw, indMin, indMax] = TrialDataUtilities.Data.getValidTimeExtents(timeCell, dataCell);
    
    nTimes = cellfun(@numel, timeCell);
    if isempty(timeDelta)
        if all(nTimes == 1)
            timeDelta = 1;
        else
            timeDelta = min(origDelta);
        end
    end
    
    % clean up small discrepancies
    tMinRaw = TrialDataUtilities.Data.removeSmallTimeErrors(tMinRaw, timeDelta, p.Results.timeReference);
    tMaxRaw = TrialDataUtilities.Data.removeSmallTimeErrors(tMaxRaw, timeDelta, p.Results.timeReference);
    
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
    if timeDelta == 0 && tMinGlobal == tMaxGlobal
        tvec = tMinGlobal;
    else
        tvec = makecol(tMinGlobal:timeDelta:tMaxGlobal);
    end
    tMinCell = tMin;
    tMaxCell = tMax;
end
