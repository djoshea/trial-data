function Q = addSegmentedSpikes(Q, spikeData, varargin)
% Q = segmentDataIntoTrials(trialInfo, spikeData)
%	Adds segmented spike data and waveforms into a Q struct with .CerebusInfo
%
%   spikeData : as returned by nevExtractSpikesEvents or with optional field .ratin

def.fieldSuffix = 'Raw'; % adds spikes to spikesRaw, waveforms to waveformsRaw, etc.
def.excludeZeroUnit = false; 
def.spikeWindowPre = 0; % time before trial start to include
def.spikeWindowPost = 0; % time after trial end to include
def.prefix = 'unit';
assignargs(def, varargin);

nTrials = length(Q);

if excludeZeroUnit
	filterInds = spikeData.unit > 0;
	spikeData.timestamp = spikeData.timestamp(filterInds);
	spikeData.electrode = spikeData.electrode(filterInds);
	spikeData.unit = spikeData.unit(filterInds);
	spikeData.waveform = spikeData.waveform(:,filterInds);
	if isfield(spikeData, 'rating')
		spikeData.rating = spikeData.rating(filterInds);
	end
end

prog = ProgressBar(nTrials, 'Segmenting spike times');

nSpikes = numel(spikeData.electrode);
trialInd = zeros(nSpikes, 1, 'uint32');
times = nan(nSpikes, 1);
for iq = 1:nTrials
    prog.update(iq);
    
    % get trial times
    zeroTime = Q(iq).CerebusInfo.startTime;
    startTime = Q(iq).CerebusInfo.startTime - spikeWindowPre;
    endTime = Q(iq).CerebusInfo.endTime + spikeWindowPost;

    % add window extensions to the CerebusInfo struct
    Q(iq).CerebusInfo.spikeWindowPre = spikeWindowPre;
    Q(iq).CerebusInfo.spikeWindowPost = spikeWindowPost;
    
    Q(iq).CerebusInfo.waveformScaleBy = spikeData.waveformScaleBy;

    % grab spike data and waveforms, subtract time offset
    mask = spikeData.timestamp >= startTime & spikeData.timestamp <= endTime;
    trialInd(mask) = iq;
    times(mask) = spikeData.timestamp(mask) - zeroTime;
end
prog.finish();
    
debug('Splitting spike times and waveforms\n');
mask = trialInd > 0;
trialInd = trialInd(mask);
times = times(mask);
electrode = spikeData.electrode(mask);
unit = spikeData.unit(mask);
waveforms = spikeData.waveform(:, mask);

[pairs, ~, idxPair] = unique([makecol(electrode) makecol(uint16(unit))], 'rows');
nPairs = size(pairs, 1);

% here we get clever in the name of efficiency, using accum array to
% group the times
times = accumarray([trialInd, idxPair], times, [nTrials, nPairs], @(x) {x});

% and even more complicated, the waves keeping columns together
nWave = size(waveforms, 1);
trialIndMat = repmat(trialInd', [nWave 1]);
idxPairMat = repmat(idxPair', [nWave 1]);
waveforms = accumarray([trialIndMat(:), idxPairMat(:)], waveforms(:), [nTrials, nPairs], @(x) {reshape(x, nWave, [])});

% get unit ratings by spike (these can change over time due to sorting epochs)
if isfield(spikeData, 'rating')
    rating = spikeData.rating(spikeInds);
    rating = accumarray([trialInd idxPair], rating, [nTrials nPairs], @(x) {x});
end
	
prog = ProgressBar(nTrials, 'Inserting segmented spike times');
if ~isempty(pairs)
    for iq = 1:nTrials
        prog.update(iq);
        for ipair = 1:nPairs
            % convert to 'elec#.unit#'
            key = LoadNev.getElecUnitStr(prefix, pairs(ipair, 1), pairs(ipair, 2));

            if ~isempty(times{iq, ipair})
                Q(iq).spikes.(key) = times{iq, ipair};
                Q(iq).waves.(key) = waveforms{iq, ipair};

                % add the rating of the first spike for this unit
                if isfield(spikeData, 'rating')
                    Q(iq).sortRatings.(key) = rating{iq, ipair}(1);
                end
            end
        end
    end
end
prog.finish();



end
