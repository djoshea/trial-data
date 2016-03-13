function [timeCellByTrial] = findThresholdCrossings(data, timeCell, thresh, lockoutPeriod)
% timeCellByTrial = findThresholdCrossings(data, timeCell, thresh, lockoutPeriod)
% either: 
%   data is trials x time, time is vector
%   data is trials cell of vectors, time is cell of time vectors
%
% thresh is threshold
% lockoutPeriod is minimum spacing between successive crossings

    if iscell(data)
        nTrials = numel(data);
    elseif ismatrix(data)
        nTrials = size(data, 1);
        assert(isvector(timeCell));
    end
     
    if nargin < 4
        lockoutPeriod = 0;
    end
    
    % utility for finding all crossing times of thresh
    function time = findCrossings(dataVec, timeVec)
        if thresh < 0 
            % falling threshold < 0
            idx = find(diff(dataVec <= thresh) == 1);
        else
            % rising threshold > 0
            idx = find(diff(dataVec >= thresh) == 1);
        end
        time = timeVec(idx);
    end
    
    if iscell(data)
        crossings = cellfun(@findCrossings, data, timeCell, 'UniformOutput', false);
    else
        crossings = arrayfun(@(trial) findCrossings(data(trial, :), timeCell), 1:nTrials, 'UniformOutput', false);
    end
    
    % ensure they are spaced by at least 1 waveform
    % loop through and extract snippets
    if lockoutPeriod <= 0
        timeCellByTrial = makecol(crossings);
        
    else
        timeCellByTrial = cell(nTrials, 1);
        for iR = 1:nTrials
            lastCrossing = -Inf;
            nCrossingsKept = 0;
            mask = true(size(crossings{iR}));
            
            for iC = 1:numel(crossings{iR})
                thisCross = crossings{iR}(iC);
                if thisCross - lastCrossing < lockoutPeriod
                    % too close to lastCrossing or to ends of data trace --> throw away
                    mask(iC) = false;
                else
                    % keep this crossing, sample snippet
                    lastCrossing = thisCross;
                    nCrossingsKept = nCrossingsKept + 1;
                end
            end

            timeCellByTrial{iR} = makecol(crossings{iR}(mask));
        end
    end
%     prog.finish();
end
        