function [data, time] = resampleTensorInTime(data, timeDim, time, timeNew, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.parse(varargin{:});
    
    assert(isvector(time));
    nTime = numel(time);
    assert(size(data, timeDim) == nTime);
    
    % figure out new time vector
    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    fixDuplicateTimes = p.Results.fixDuplicateTimes;
    interpolate = p.Results.interpolate;
    %interpolateMethod = p.Results.interpolateMethod;
   
     % compute the global min / max timestamps or each trial
     [tMinRaw, tMaxRaw] = TrialDataUtilities.Data.getValidTimeExtents(timeCell, dataCell);
%      [tMinRaw, tMaxRaw] = cellfun(@minmax, timeCell);

     nTimes = cellfun(@numel, timeCell);

     if isempty(timeDelta)
         if all(nTimes == 1)
            timeDelta = 1;
         else
             % auto-compute time delta
             %warning('Auto-computing time-delta from timeseries. Specify time delta for consistent results');

             if fixDuplicateTimes
                 timeCell = cellfun(@fixDup, timeCell, 'UniformOutput', false);
             end
             deltaCell = cellfun(@(x) nanmedian(diff(x)), timeCell);
             % compute the median delta for each channel
             timeDelta = nanmedian(deltaCell(:), 1);
             % and use the minimum spacing
             timeDelta = double(nanmin(timeDelta));
         end
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

        % not sure how this was supposed to work
%         tol = 1e-9;
%         if  (any(abs(tMin(:) - tMinRaw(:)) > tol) || any(abs(tMax(:) - tMaxRaw(:)) > tol))
%             error('Timestamps do not align with timeReference. Set ''interpolate'' to true');
%         end
    end

    % build the global time vector
    tMinGlobal = nanmin(tMin);
    tMaxGlobal = nanmax(tMax);
    tvec = makecol(tMinGlobal:timeDelta:tMaxGlobal);
    tMinCell = tMin;
    tMaxCell = tMax;
    

    tensor = TensorUtils.shiftdimToFirstDim(tensor, timeDim);
    tensor = interp1(tensor, time, timeNew
    
    N = size(dataCell, 1);
    trialValid = p.Results.trialValid;
    if isempty(trialValid)
        trialValid = truevec(N);
    end
   
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
    
    G = size(dataCell, 2);
    
    mat = nan([N, T, C, G]); % we'll reshape this later

    indStart = floor(((tMin - tMinGlobal) / timeDelta) + 1);
    indStop  = floor(((tMax - tMinGlobal) / timeDelta) + 1);
    
%     if p.Results.showProgress
%         prog = ProgressBar(N, 'Embedding data over trials into common time vector');
%     end
    for i = 1:N
        if ~trialValid(i), continue; end
%         if mod(i, 100) && p.Results.showProgress, prog.update(i); end
        for g = 1:G
            if ~isnan(indStart(i,g)) && ~isnan(indStop(i,g))
                if numel(indStart(i,g):indStop(i,g)) > 1
                    mask = ~all(isnan(dataCell{i, g}), 2);
                    if assumeUniformSampling
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
%     if p.Results.showProgress
%         prog.finish();
%     end
    
    % pare down time points from edges with insufficient trial counts
    if p.Results.minTrials > 0 || p.Results.minTrialFraction > 0
        nTrialsOverTime = sum(all(~isnan(mat(trialValid, :, :)), 3), 1);
        minTrials = max(p.Results.minTrials, nnz(trialValid) * p.Results.minTrialFraction);
        
        tMask= falsevec(T);
        tMask(find(nTrialsOverTime >= minTrials, 1, 'first') : find(nTrialsOverTime >= minTrials, 1, 'last')) = true;
    else
        tMask = truevec(T);
    end
    
    mat = reshape(mat(:, tMask, :, :), [N nnz(tMask) C*G]);
    tvec = tvec(tMask);
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