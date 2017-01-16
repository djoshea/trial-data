function [timeCrossLow, idxCrossLow] = findThresholdCrossingsLowThenHigh(data, time, threshLow, threshHigh)
    % [timeCrossHighCell, timeCrossLowCell] = findThresholdCrossingsLowThenHigh(data, time, threshLow, threshHigh)
    % Find threshold crossings above high threshold, then return the time at
    % which the data most proximally crosses the low threshold. This enables a
    % very low "low threshold" which does not capture transients that do not
    % proceed to cross the high threshold.
    %
    % either: 
    %   data is time x trials matrix, time is vector
    %   data is trials cell of vectors, time is cell of time vectors
    % threshHighLow is [threshHigh threshLow]
    %
    % data is nTrials x 1 cell or nTrials x 1
    
    if iscell(data)
        nTrials = numel(data);
    elseif ismatrix(data)
        nTrials = size(data, 2);
        assert(isvector(time));
    end
    
    lo = threshLow;
    hi = threshHigh;
    
    if iscell(data)
        [timeCrossLow, idxCrossLow] = cellfun(@doThresh, data, time, 'UniformOutput', true);
    else
        [timeCrossLow, idxCrossLow] = arrayfun(@(trial) doThresh(data(:, trial), time), 1:nTrials, 'UniformOutput', true);
    end

    function [timeCrossLow, idxCrossLow] = doThresh(sig, time)
        idxCrossHigh = find(sig > hi, 1, 'first');
        if isempty(idxCrossHigh)
            idxCrossLow = NaN;
            timeCrossLow = NaN;
            return;
        end
        sig(idxCrossHigh : end) = Inf;
        idxCrossLow = find(sig < lo, 1, 'last');
        if isempty(idxCrossLow)
            idxCrossLow = NaN;
            timeCrossLow = NaN;
            return;
        end
        timeCrossLow = time(idxCrossLow);
    end

end