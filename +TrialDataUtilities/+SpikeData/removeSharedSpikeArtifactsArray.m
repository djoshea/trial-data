function [td, artifactCounts] = removeSharedSpikeArtifactsArray(td, array, varargin)
% detects spikes coincident on at least a certain number or fraction of
% trials and removes them from ALL channels. This is for removing
% corrupting artifacts, e.g. from a solenoid.

p = inputParser();
p.addOptional('minChannels', 5, @isscalar);
p.addOptional('minFracChannels', 0, @isscalar);
p.addOptional('timeWindow', 0.1, @isscalar);
p.addParameter('keepRemovedSpikes', false, @islogical);    

p.parse(varargin{:});

cd = td.getChannelDescriptor(array);
nTrials = td.nTrials;
% check all channels against all other channels
nC = cd.nChannels;

threshC = max(p.Results.minChannels, ceil(p.Results.minFracChannels * nC));

prog = ProgressBar(nC*td.nTrials, 'Scanning channel pairs for shared artifact');
prog.enableParallel();

td = td.reset();

% nTrials x C cell of spike times
dataCheck = td.getSpikeTimes(array);
data = dataCheck;

artifactCounts = nan(size(data));
timeWindowHalf = p.Results.timeWindow / 2;

    
maskKeep = cell(td.nTrials, nC);

parfor c = 1:nC
    prog.update(c);
    
    matchCounts = cellfun(@(times) zeros(size(times)), data(:, c), 'UniformOutput', false);
    
%     progI = ProgressBar(td.nTrials, 'Scanning for matches over trials');
    for t = 1:nTrials
        thisTrial = data{t, c};
        
        if ~isempty(thisTrial)
            for oc = c+1:nC
                thisTrialOC = dataCheck{t, oc};

                % detect my matches on this other channel, i.e. whether i have
                % a matching spike within timeWindow/2 of each of my spikes
                if ~isempty(thisTrialOC)
                    [minDist, ~] = pdist2(thisTrialOC, thisTrial, 'cityblock', 'Smallest', 1);
                    maskHasMatch = minDist' <= timeWindowHalf;
                    % increment count for this channel
                    matchCounts{t} = matchCounts{t} + maskHasMatch;
                end
            end

            maskKeep{t, c} = matchCounts{t} < threshC; 
        else
            maskKeep{t, c} = [];
        end
        artifactCounts(t, c) = nnz(~maskKeep{t, c});
    end
end
prog.finish();

for c = 1:nC
    fprint('\t%d artifacts in %d_%d\n', sum(artifactCounts(:, c)), cd.electrodes(c), cd.units(c));
end

if p.Results.keepRemovedSpikes
    keepAs = append(array, 'Artifact');
else
    keepAs = '';
end

td = td.maskSpikeChannelSpikes(array, maskKeep, 'keepRemovedSpikes', keepAs);    