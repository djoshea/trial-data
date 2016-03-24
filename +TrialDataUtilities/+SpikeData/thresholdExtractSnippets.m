function [timeCellByTrial, waveformsByTrial, idxCellByTrial] = thresholdExtractSnippets(data, time, thresh, snippetWindow, varargin)
% [timeCellByTrial, waveformsByTrial] = thresholdExtractSnippets(data, thresh, snippetWindow)
% either: 
%   data is time x trials, time is vector
%   data is trials cell of vectors, time is cell of time vectors
%
% thresh is threshold
% snippetWindow is [samplesPreCrossing, samplesPostCrossing]
% total snippet
% length will be sum(snippetWindow)
% waveformsByTrial will be trials x 1 cell, each is nSpikes x nSamples

    p = inputParser();
    p.addParameter('mode', 'leftToRight', @ischar); % or 'largestFirst'
    p.addParameter('lockout', snippetWindow, @isscalar);
    p.addParameter('extremumWithin', 7, @islogical);
    p.parse(varargin{:});
    lockout = p.Results.lockout;
    extremumWithin = p.Results.extremumWithin;
    
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
    idxCellByTrial = cell(nTrials, 1);
    timeCellByTrial = cell(nTrials, 1);
    waveformsByTrial = cell(nTrials, 1);
    
    if strcmp(p.Results.mode, 'leftToRight')
        % sweep left to right, thresholding as you go
        prog = ProgressBar(nTrials, 'Extracting threhold crossings from continuous data');
        for iR = 1:nTrials
            prog.update(iR);
            lastCrossing = -Inf;
            nCrossingsKept = 0;
            mask = true(size(crossings{iR}));
            waveMat = nan(numel(crossings{iR}), snippetLength);

            if iscell(data)
                nTimeThisTrial = numel(data{iR});
            else
                nTimeThisTrial = size(data, 1);
            end

            for iC = 1:numel(crossings{iR})
                thisCross = crossings{iR}(iC);
                if thisCross - lastCrossing < lockout || thisCross <= snippetPre || thisCross + snippetPost - 1 > nTimeThisTrial
                    % too close to lastCrossing or to ends of data trace --> throw away
                    mask(iC) = false;
                else
                    % keep this crossing, sample snippet
                    lastCrossing = thisCross;
                    nCrossingsKept = nCrossingsKept + 1;
                    if iscell(data)
                        waveMat(nCrossingsKept, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1))';
                    else
                        waveMat(nCrossingsKept, :) = data(thisCross + (-snippetPre : snippetPost-1), iR)';
                    end
                end
            end

            if iscell(time)
                tvec = time{iR};
            else
                tvec = time;
            end
            idxCellByTrial{iR} = makecol(crossings{iR}(mask));
            timeCellByTrial{iR} = makecol(tvec(crossings{iR}(mask)));
            waveformsByTrial{iR} = waveMat(1:nCrossingsKept, :);
        end
        prog.finish();
        
    elseif strcmp(p.Results.mode, 'largestFirst')
        % take biggest waveforms first, and ignore those within lockout
        % window
        prog = ProgressBar(nTrials, 'Extracting threhold crossings from continuous data');
        for iR = 1:nTrials
            prog.update(iR);
            cross = crossings{iR};
            C = numel(cross);
            if iscell(data)
                nTimeThisTrial = numel(data{iR});
            else
                nTimeThisTrial = size(data, 1);
            end
            waveMat = nan(C, snippetLength);
            
            % figure out the extreme value of each snipet
            waveExt = nanvec(C);
            for iC = 1:C
                thisCross = cross(iC);    
                if thisCross - snippetPre < 1, continue, end
                if thisCross + snippetPost - 1 > nTimeThisTrial, continue; end
                if iscell(data)
                    snip = data{iR}(thisCross+(0 : extremumWithin-1));
                    waveMat(iC, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1))';
                else
                    snip = data(thisCross+(0 : extremumWithin-1), iR)';
                    waveMat(iC, :) = data(thisCross + (-snippetPre : snippetPost-1), iR)';
                end
                if thresh > 0
                    waveExt(iC) = max(snip);
                else
                    waveExt(iC) = min(snip);
                end
            end
            
            maskPicked = falsevec(C);
            maskEligible = truevec(C);
            
            while(any(maskEligible))
                % waveExt(ineligible) will be NaN
                [~, pick] = max(abs(waveExt));
                    
                if isnan(waveExt(pick))
                    break;
                end
                maskPicked(pick) = true;
                maskEligible(pick) = false;
                
                % mark ineligible other crossings within lockout window
                tooClose = abs(cross - cross(pick)) < lockout;
                maskEligible(tooClose) = false;
                waveExt(tooClose) = NaN;
            end
            if iscell(time)
                tvec = time{iR};
            else
                tvec = time;
            end
            idxCellByTrial{iR} = makecol(crossings{iR}(maskPicked));
            timeCellByTrial{iR} = makecol(tvec(crossings{iR}(maskPicked)));
            waveformsByTrial{iR} = waveMat(maskPicked, :);
        end
        prog.finish();
    else
        error('Unknown mode');
    end
        
end
        