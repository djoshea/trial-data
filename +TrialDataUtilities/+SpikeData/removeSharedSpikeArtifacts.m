function [td, artifactCounts] = removeSharedSpikeArtifacts(td, varargin)

p = inputParser();
p.addOptional('unitNames', td.listSpikeChannels, @iscellstr);
p.addOptional('minChannels', 5, @isscalar);
p.addOptional('minFracChannels', 0, @isscalar);
p.addOptional('timeWindow', 0.5, @isscalar);
p.parse(varargin{:});

unitNames = p.Results.unitNames;
if isempty(unitNames)
    unitNames = td.listSpikeChannels();
end

C = numel(unitNames);

threshC = max(p.Results.minChannels, ceil(p.Results.minFracChannels * C));

prog = ProgressBar(C, 'Scanning channels for shared artifact');

td = td.reset();

% nTrials x C cell of spike times
data = td.getRawSpikeTimes(td.listSpikeChannels());
artifactCounts = nan(size(data));
timeWindowHalf = p.Results.timeWindow / 2;
matchCounts = cellfun(@(times) zeros(size(times)), data, 'UniformOutput', false);
    
for c = 1:C
    unitName = unitNames{c};
    prog.update(c, 'Scanning %s for shared artifacts', unitName);
    
    maskKeep = cellvec(td.nTrials);
    
%     progI = ProgressBar(td.nTrials, 'Scanning for matches over trials');
    for t = 1:td.nTrials
%         progI.update(t);
        
        thisTrial = data{t, c};
        
        if ~isempty(thisTrial)
            for oc = c+1:C
                thisTrialOC = data{t, oc};

                % detect my matches on this other channel, i.e. whether i have
                % a matching spike within timeWindow/2 of each of my spikes
                if ~isempty(thisTrialOC)
                    [minDist, idxOC] = pdist2(thisTrialOC, thisTrial, 'cityblock', 'Smallest', 1);
                    maskHasMatch = minDist' <= timeWindowHalf;
                    % increment count for this channel
                    matchCounts{t, c} = matchCounts{t, c} + maskHasMatch;
                    % increment count for other channel
                    matchCounts{t, oc}(idxOC(maskHasMatch)) = matchCounts{t, oc}(idxOC(maskHasMatch)) + 1;
                end
            end

            maskKeep{t} = matchCounts{t, c} < threshC; 
        else
            maskKeep{t} = [];
        end
        artifactCounts(t, c) = nnz(~maskKeep{t});
    end
%     progI.finish();
    
    td = td.maskSpikeChannelSpikesRaw(unitName, maskKeep);    
end
prog.finish();