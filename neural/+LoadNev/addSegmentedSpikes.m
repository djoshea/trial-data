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
[spikeMaps waveforms] = deal(cell(nTrials,1));

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
for iq = 1:nTrials
    prog.update(iq);
    
    % get trial times
    zeroTime = Q(iq).CerebusInfo.startTime;
    startTime = Q(iq).CerebusInfo.startTime - spikeWindowPre;
    endTime = Q(iq).CerebusInfo.endTime + spikeWindowPost;

    % add window extensions to the CerebusInfo struct
    Q(iq).CerebusInfo.spikeWindowPre = spikeWindowPre;
    Q(iq).CerebusInfo.spikeWindowPost = spikeWindowPost;

    % grab spike data and waveforms, subtract time offset
    spikeInds = spikeData.timestamp >= startTime & spikeData.timestamp <= endTime;
	times = spikeData.timestamp(spikeInds) - zeroTime;
	electrode = spikeData.electrode(spikeInds);
	unit = spikeData.unit(spikeInds);
	waveforms = spikeData.waveform(:, spikeInds);

    % get unit ratings by spike (these can change over time due to sorting epochs)
    if isfield(spikeData, 'rating')
        rating = spikeData.rating(spikeInds);
    end

	% build map from 'elec#.unit#' to spike times and to waveforms
	spikeMaps{iq} = containers.Map('ValueType', 'any', 'KeyType', 'char');
	waveformMaps{iq} = containers.Map('ValueType', 'any', 'KeyType', 'char');
    if isfield(spikeData, 'rating')
        ratingMaps{iq} = containers.Map('ValueType', 'any', 'KeyType', 'char');
    end

	pairs = unique([makecol(electrode) makecol(uint16(unit))], 'rows');
    if ~isempty(pairs)
        for ipair = 1:size(pairs,1)
            % convert to 'elec#.unit#'
            key = LoadNev.getElecUnitStr(prefix, pairs(ipair, 1), pairs(ipair, 2));

            % add spikeMaps for this electrode / unit pair 
            matches = electrode==pairs(ipair, 1) & unit==pairs(ipair,2);
           
            % add this unit to the maps if it has any spikes
            if nnz(matches) > 0
                spikeMaps{iq}(key) = times(matches);
                waveformMaps{iq}(key) = waveforms(:, matches);

                % add the rating of the first spike for this unit
                if isfield(spikeData, 'rating')
                    ratingMaps{iq}(key) = rating(find(matches,1));
                end
            end
        end
    end
end
prog.finish();

% add them to the struct array
if ~isempty(Q)
    [Q.(['spikes' fieldSuffix])] = deal(spikeMaps{:});
    [Q.(['waveforms' fieldSuffix])] = deal(waveformMaps{:});
    
    if isfield(spikeData, 'rating')
        [Q.(['sortRatings' fieldSuffix])] = deal(ratingMaps{:});
    end
end


end
