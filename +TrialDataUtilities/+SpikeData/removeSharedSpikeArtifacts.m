function [td, artifactCounts] = removeSharedSpikeArtifacts(td, varargin)
% detects spikes coincident on at least a certain number or fraction of
% trials and removes them from ALL channels. This is for removing
% corrupting artifacts, e.g. from a solenoid.

p = inputParser();
p.addOptional('unitNames', td.listSpikeChannels, @iscellstr);
p.addOptional('minChannels', 5, @isscalar);
p.addOptional('minFracChannels', 0, @isscalar);
p.addOptional('timeWindow', 0.1, @isscalar);
p.addParameter('keepRemovedSpikes', false, @islogical);
p.addParameter('checkUnitNames', td.listSpikeChannels, @iscellstr);            

p.parse(varargin{:});


unitNames = p.Results.unitNames;
if isempty(unitNames)
    unitNames = td.listSpikeChannels();
end

checkUnitNames = p.Results.checkUnitNames;
if isempty(checkUnitNames)
    checkUnitNames = td.listSpikeChannels();
end

nC = numel(unitNames);
nOC = numel(checkUnitNames);

processed = false(nC, nOC);

threshC = max(p.Results.minChannels, ceil(p.Results.minFracChannels * nC));

prog = ProgressBar(nC, 'Scanning channels for shared artifact');

td = td.reset();

% nTrials x C cell of spike times
dataCheck = td.getSpikeTimes(checkUnitNames);

if ~isequal(checkUnitNames, unitNames)
    data = td.getSpikeTimes(unitNames);
else
    data = dataCheck;
end
artifactCounts = nan(size(data));
timeWindowHalf = p.Results.timeWindow / 2;
matchCounts = cellfun(@(times) zeros(size(times)), data, 'UniformOutput', false);
    
for c = 1:nC
    unitName = unitNames{c};
    prog.update(c, 'Scanning %s for shared artifacts', unitName);
    
    maskKeep = cellvec(td.nTrials);
    
%     progI = ProgressBar(td.nTrials, 'Scanning for matches over trials');
    for t = 1:td.nTrials
%         progI.update(t);
        
        thisTrial = data{t, c};
        
        if ~isempty(thisTrial)
            for oc = 1:nOC
                if processed(c, oc), continue; end
                
                checkName = checkUnitNames{oc};
                if isequal(checkName, unitName), continue; end
                
                thisTrialOC = dataCheck{t, oc};

                % detect my matches on this other channel, i.e. whether i have
                % a matching spike within timeWindow/2 of each of my spikes
                if ~isempty(thisTrialOC)
                    [minDist, idxOC] = pdist2(thisTrialOC, thisTrial, 'cityblock', 'Smallest', 1);
                    maskHasMatch = minDist' <= timeWindowHalf;
                    % increment count for this channel
                    matchCounts{t, c} = matchCounts{t, c} + maskHasMatch;
                    
                    % check if we can count this for the reverse direcition too, where we were removing
                    % spikes from unitName <--> checkName 
                    % increment count for other channel
                    [tf1, idx_check_in_unit] = ismember(checkName, unitNames);
                    [tf2, idx_unit_in_check] = ismember(unitName, checkUnitNames);
                    if tf1 && tf2
                        matchCounts{t, idx_check_in_unit}(idxOC(maskHasMatch)) = matchCounts{t, idx_check_in_unit}(idxOC(maskHasMatch)) + 1;
                        processed(idx_check_in_unit, idx_unit_in_check) = true;
                    end
                end
            end

            maskKeep{t} = matchCounts{t, c} < threshC; 
        else
            maskKeep{t} = [];
        end
        artifactCounts(t, c) = nnz(~maskKeep{t});
    end
    processed(c, :) = true;
%     progI.finish();

    debug('%d artifacts in %s\n', sum(artifactCounts(:, c)), unitName);
    
    if p.Results.keepRemovedSpikes
        [a,e] = SpikeChannelDescriptor.parseArrayElectrodeUnit(unitName);
        keepAs = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(a, e, 255);
    else
        keepAs = '';
    end

    td = td.maskSpikeChannelSpikes(unitName, maskKeep, 'keepRemovedSpikes', keepAs);    
end
prog.finish();