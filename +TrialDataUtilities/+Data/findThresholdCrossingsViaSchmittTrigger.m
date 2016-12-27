function [timeCrossHighCell, timeCrossLowCell, dataThresh] = findThresholdCrossingsViaSchmittTrigger(data, time, threshHighLow, varargin)
% [timeCrossHighCell, timeCrossLowCell, dataThresh] = findThresholdCrossingsViaSchmittTrigger(data, time, threshHighLow, varargin)
%
% either: 
%   data is time x trials matrix, time is vector
%   data is trials cell of vectors, time is cell of time vectors
% threshHighLow is [threshHigh threshLow]
%
% data is nTrials x 1 cell or nTrials x 
    if iscell(data)
        nTrials = numel(data);
    elseif ismatrix(data)
        nTrials = size(data, 2);
        assert(isvector(time));
    end
    
    lo = threshHighLow(2);
    hi = threshHighLow(1);
    
    if iscell(data)
        [dataThresh, timeCrossHighCell, timeCrossLowCell] = cellfun(@schmidtTrigger, data, time, 'UniformOutput', false);
    else
        [dataThresh, timeCrossHighCell, timeCrossLowCell] = arrayfun(@(trial) schmidtTrigger(data(:, trial), time), 1:nTrials, 'UniformOutput', false);
    end

    function [thresh, timeCrossHigh, timeCrossLow] = schmidtTrigger(sig, time)
        state = 0;
        thresh = false(numel(sig), 1);
        for i = 1:numel(sig)
            if sig(i) >= hi
                state = 1;
            elseif sig(i) <= lo
                state = 0;
            end

            thresh(i) = state;
        end

        timeCrossHigh = time(diff(thresh) == 1);
        timeCrossLow = time(diff(thresh) == -1);
    end

end