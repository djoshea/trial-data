function [spikes, waveforms, idxCell] = thresholdExtractSnippetsGroup(data, time, thresh, snippetWindow, varargin)
% [timeCellByTrial, waveformsByTrial] = thresholdExtractSnippets(data, thresh, snippetWindow)
%   data is trials cell of times x nChannels
%   time is trials cell of times x 1
%   thresh is trials x nChannels
%
%   spike, waveforms, and idxCell are trials x nChannels cell
%
% thresh is threshold
% snippetWindow is [samplesPreCrossing, samplesPostCrossing]
% total snippet
% length will be sum(snippetWindow)
% timeCellByTrial, waveformsByTrial will be trials x nUnits cell, each is nSpikes x nSamples

    p = inputParser();
    p.addParameter('mode', 'leftToRight', @ischar); % or 'largestFirst'
    p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
    p.addParameter('extremumWithin', 7, @islogical);
    p.parse(varargin{:});
    if isempty(p.Results.lockoutPrePost)
      lockoutPrePost = snippetWindow;
    else
      lockoutPrePost = p.Results.lockoutPrePost;
    end
    lockoutPre = lockoutPrePost(1);
    lockoutPost = lockoutPrePost(2);
    extremumWithin = p.Results.extremumWithin;

    nTrials = size(data, 1);
    emptyMask = cellfun(@isempty, data);
    nonEmpty = find(~emptyMask, 1);
    if ~any(nonEmpty)
        error('No non-empty data elements');
    end
    nChannels = size(data{nonEmpty}, 2);
    
    assert(size(thresh, 1) == nTrials);
    assert(size(thresh, 2) == nChannels);
    
    snippetPre = snippetWindow(1);
    snippetPost = snippetWindow(2);
    snippetLength = snippetPre + snippetPost;

    % utility for finding all crossing times of thresh
    function idx = findCrossings(dataVec, thresh)
        if thresh < 0
            % falling threshold < 0
            idx = find(diff(dataVec <= thresh) == 1);
        else
            % rising threshold > 0
            idx = find(diff(dataVec >= thresh) == 1);
        end
    end

    crossings = cell(nTrials, nChannels);
    prog = ProgressBar(nTrials, 'Finding threhold crossings in continuous data');
    for iR = 1:nTrials
        prog.update(iR);
        for iC = 1:nChannels
            crossings{iR, iC} = findCrossings(data{iR}(:, iC), thresh(iR, iC));
        end
    end
    prog.finish();
    
    % ensure they are spaced by at least 1 waveform
    % loop through and extract snippets
    [spikes, waveforms, idxCell] = deal(cell(nTrials, nChannels));
    
    if strcmp(p.Results.mode, 'leftToRight')
        % sweep left to right, thresholding as you go
        prog = ProgressBar(nTrials, 'Extracting threhold crossings from continuous data');
        for iR = 1:nTrials
            prog.update(iR);
            
            for iC = 1:nChannels
                lastCrossing = -Inf;
                nCrossingsKept = 0;
                mask = true(size(crossings{iR, iC}));
                waveMat = nan(numel(crossings{iR, iC}), snippetLength);

                nTimeThisTrial = numel(data{iR}(:, iC));
                
                for iT = 1:numel(crossings{iR})
                    thisCross = crossings{iR}(iT);
                    if thisCross - lastCrossing < lockoutPost || thisCross <= snippetPre || thisCross + snippetPost - 1 > nTimeThisTrial
                        % too close to lastCrossing or to ends of data trace --> throw away
                        mask(iT) = false;
                    else
                        % keep this crossing, sample snippet
                        lastCrossing = thisCross;
                        nCrossingsKept = nCrossingsKept + 1;
                        waveMat(nCrossingsKept, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1), iC)';
                    end
                end
                tvec = time{iR};
                idxCell{iR, iC} = makecol(crossings{iR, iC}(mask));
                spikes{iR, iC} = makecol(tvec(crossings{iR, iC}(mask)));
                waveforms{iR, iC} = waveMat(1:nCrossingsKept, :);
            end
        end
        prog.finish();

    elseif strcmp(p.Results.mode, 'largestFirst')
        % take biggest waveforms first, and ignore those within lockout
        % window
        prog = ProgressBar(nTrials, 'Extracting threhold crossings from continuous data');
        for iR = 1:nTrials
            prog.update(iR);
            for iC = 1:nChannels
                cross = crossings{iR, iC};
                C = numel(cross);
                nTimeThisTrial = numel(data{iR}(:, iC));
                waveMat = nan(C, snippetLength);

                % figure out the extreme value of each snipet
                waveExt = nanvec(C);
                for iT = 1:C
                    thisCross = cross(iT);
                    if thisCross - snippetPre < 1, continue, end
                    if thisCross + snippetPost - 1 > nTimeThisTrial, continue; end
                    snip = data{iR}(thisCross+(0 : extremumWithin-1));
                    waveMat(iT, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1), iC)';
                    if thresh > 0
                        waveExt(iT) = max(snip);
                    else
                        waveExt(iT) = min(snip);
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
                    tooClose = cross - cross(pick) >= -lockoutPre & cross-cross(pick) < lockoutPost;
                    maskEligible(tooClose) = false;
                    waveExt(tooClose) = NaN;
                end
                tvec = time{iR};
                idxCell{iR, iC} = makecol(crossings{iR, iC}(maskPicked));
                spikes{iR, iC} = makecol(tvec(crossings{iR, iC}(maskPicked)));
                waveforms{iR, iC} = waveMat(maskPicked, :);
            end
        end
        prog.finish();
    else
        error('Unknown mode');
    end

end
