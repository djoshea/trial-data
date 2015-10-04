function [timeCellByTrial, waveformsByTrial] = thresholdExtractSnippets(data, tvec, thresh, snippetWindow)
% [timeCellByTrial, waveformsByTrial] = thresholdExtractSnippets(data, thresh, snippetWindow)
% either: 
%   data is time x trials, time is vector
%   data is trials cell of vectors, time is cell of time vectors
%
% thresh is threshold
% snippetWindow is [samplesPreCrossing, samplesPostCrossing]
% total snippet
% length will be sum(snippetWindow)

    if iscell(data)
        nTrials = numel(data);
    elseif ismatrix(data)
        nTrials = size(data, 2);
    end
    
    snippetPre = snippetWindow(1);
    snippetPost = snippetWindow(2);
    snippetLength = snippetPre + snippetPost;
    
    % utility for finding all crossing times of thresh
    function idx = findCrossings(dataVec)
        if thresh < 0 
            % falling threshold < 0
            idx = find(diff(dataVec <= thresh) == 1);
        else
            % rising threshold > 0
            idx = find(diff(dataVec >= thresh) == 1);
        end
    end
    
    if iscell(data)
        crossings = cellfun(@findCrossings, data, 'UniformOutput', false);
    else
        crossings = arrayfun(@(trial) findCrossings(data(:, trial)), 1:nTrials, 'UniformOutput', false);
    end
    
    % ensure they are spaced by at least 1 waveform
    % loop through and extract snippets
    timeCellByTrial = cell(nTrials, 1);
    waveformsByTrial = cell(nTrials, 1);
    prog = ProgressBar(nTrials, 'Extracting threhold crossings from continuous data');
    for iR = 1:nTrials
        prog.update(iR);
        lastCrossing = -Inf;
        nCrossingsKept = 0;
        mask = true(size(crossings{iR}));
        waveMat = nan(snippetLength, numel(crossings{iR}));
        
        if iscell(data)
            nTimeThisTrial = numel(data{iR});
        else
            nTimeThisTrial = size(data, 1);
        end
        
        for iC = 1:numel(crossings{iR})
            thisCross = crossings{iR}(iC);
            if thisCross - lastCrossing < snippetLength || thisCross <= snippetPre || thisCross + snippetPost - 1 > nTimeThisTrial
                % too close to lastCrossing or to ends of data trace --> throw away
                mask(iC) = false;
            else
                % keep this crossing, sample snippet
                lastCrossing = thisCross;
                nCrossingsKept = nCrossingsKept + 1;
                if iscell(data)
                    waveMat(:, nCrossingsKept) = data{iR}(thisCross + (-snippetPre : snippetPost-1));
                else
                    waveMat(:, nCrossingsKept) = data(thisCross + (-snippetPre : snippetPost-1), iR);
                end
            end
        end
        
        timeCellByTrial{iR} = makecol(tvec(crossings{iR}(mask)));
        waveformsByTrial{iR} = waveMat(:, 1:nCrossingsKept);
    end
    prog.finish();
end
        