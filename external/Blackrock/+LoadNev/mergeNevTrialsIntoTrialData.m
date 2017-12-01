function [td, numMerged] = mergeNevTrialsIntoTrialData(td, Q, varargin)
% given a Q struct returned by loadNevNoSerial and a ChannelReach trial data instance,
% merge the spikes from the nev into the trial data with the same unit
% names. spikes from the td will not be deleted, but they will be
% overwritten on merged trials

p = inputParser;
p.addParameter('includePhotobox', false, @islogical); % copy Q photobox to td.photoboxNev
p.addParameter('overwriteSpikes', true, @islogical);
p.addParameter('includeWaveforms', true, @islogical);
p.addParameter('includeLFP', true, @islogical);
p.addParameter('includeBroadband', true, @islogical);
p.parse(varargin{:});
includeSpikes = p.Results.overwriteSpikes;
includeWaveforms = p.Results.includeWaveforms;
includeLFP = p.Results.includeLFP;
includeBroadband = p.Results.includeBroadband;

% align td to Q based on computed delay periods
qMatchInTD = LoadNev.findNevTrialsMatchInTrialData(td, Q);

prog = ProgressBar(numel(Q), 'Merging nev trials into TrialData');
nTD = td.nTrials;
spikeData = struct();
waveData = struct();
lfpData = cellvec(nTD);
lfpTime = cellvec(nTD);
broadbandData = cellvec(nTD);
broadbandTime = cellvec(nTD);

td = td.reset();
nevShort = cellvec(nTD);

convertName = @(dotName) strrep(dotName, '.', '_');

mergedMask = falsevec(td.nTrials);

for iQ = 1:numel(Q)
    prog.update(iQ);
    if isnan(qMatchInTD(iQ))
        continue;
    end
    
    q = Q(iQ);
    iR = qMatchInTD(iQ);
    mergedMask(iR) = true;
    units = q.spikesRaw.keys;
    
    tOffsetQ = 0;

    if ~isempty(nevShort{iR})
        warning('Trial %d in TrialData has already been merged with file %s', iR, nevShort{iR});
    end 
    nevShort{iR} = q.CerebusInfo.nevNameShort;

%     % get photobox
%     photoboxNev{iR} = q.analog.photobox;
%     photoboxNevTime{iR} = q.analog.time;
%     photoboxScaleLims = q.analog.scaleLims;
    
    % grab spikes into spikeData struct
    for iK = 1:numel(units)
        chName = convertName(units{iK});
        if ~isfield(spikeData, chName) % add field if needed
            spikeData.(chName) = cellvec(nTD);
            waveData.(chName) = cellvec(nTD);
        end
            
        spikeData.(chName){iR} = q.spikesRaw(units{iK}) - tOffsetQ;
        if includeWaveforms
            waveData.(chName){iR} = q.waveformsRaw(units{iK});
        end
    end
    
    % overwrite lfp
    if isfield(q, 'lfp') && includeLFP
        lfpData{iR} = q.lfp.data'; % need each channel as successive column
        lfpTime{iR} = makecol(q.lfp.time - tOffsetQ);
        lfpScaleFromLims = q.lfp.scaleLims(1:2);
        lfpScaleToLims = q.lfp.scaleLims(3:4);
    end
    
    if isfield(q, 'broadband') && includeBroadban
        broadbandData{iR} = q.broadband.data'; % need each channel as successive column
        broadbandTime{iR} = makecol(q.broadband.time - tOffsetQ);
        broadbandScaleFromLims = q.broadband.scaleLims(1:2);
        broadbandScaleToLims = q.broadband.scaleLims(3:4);
    end
end
prog.finish();

% if p.Results.includePhotobox
%     td = td.addAnalog('photoboxNev', photoboxNev, photoboxNevTime, 'units', 'V', ...
%         'isAligned', false, 'clearForInvalid', false, ...
%         'scaleFromlims', photoboxScaleLims(1:2), 'scaleToLims', photoboxScaleLims(3:4));
% end

numMerged = nnz(mergedMask);
if numMerged == 0
    return
end

tdMask = removenan(qMatchInTD);

if includeSpikes
    units = fieldnames(spikeData);
    prog = ProgressBar(numel(units), 'Adding spike data to TrialData');
    for iU = 1:numel(units)
        chName = strrep(units{iU}, 'unit', 'ch');
        prog.update(iU, 'Adding %s to TrialData', chName);
        if includeWaveforms
            wd = waveData.(units{iU});
            nonEmpty = find(~cellfun(@isempty, wd), 1);
            if isempty(nonEmpty)
                wd = [];
            else
                % swap dims of each waveform array (samples along dim 2)
                wd = cellfun(@(x) x', wd, 'UniformOutput', false);
                nWave = size(wd{nonEmpty}, 2);
                wavetvec = (-10 : nWave-11) / 30;
            end
            td = td.addOrUpdateSpikeChannel(chName, spikeData.(units{iU}), 'mask', tdMask, 'waveforms', wd, 'waveformsTime', wavetvec);
        else
            td = td.addOrUpdateSpikeChannel(chName, spikeData.(units{iU}), 'mask', tdMask);
        end
    end
    prog.finish();
end

if includeLFP && isfield(q, 'lfp')
    chLookup = Q(1).lfp.lookup;
    prefix = 'lfp';
    lfpCh = arrayfun(@(ch) sprintf('%s_ch%03d', prefix, ch), chLookup, 'UniformOutput', false);

    debug('Adding %s data to TrialData as continuous neural group\n', prefix);
    td = td.addOrUpdateContinuousNeuralChannelGroup(prefix, lfpCh, ...
     lfpData, lfpTime, 'mask', tdMask, 'units', 'uV', 'scaleFromLims', lfpScaleFromLims, 'scaleToLims', lfpScaleToLims, ...
     'dataInMemoryScale', true);
    
    hasLfp = cellfun(@(x) size(x, 2) > 0, lfpData);
    td = td.addOrUpdateBooleanParam(sprintf('has%s', upperFirst(prefix)), hasLfp, 'mask', tdMask);
end
if includeBroadband && isfield(q, 'broadband')
    chLookup = Q(1).broadband.lookup;
    prefix = 'broadband';
    broadbandCh = arrayfun(@(ch) sprintf('%s_ch%03d', prefix, ch), chLookup, 'UniformOutput', false);
        
    debug('Adding %s data to TrialData as continuous neural group\n', prefix);
    td = td.addOrUpdateContinuousNeuralChannelGroup(prefix, broadbandCh, ...
        broadbandData, broadbandTime, 'mask', tdMask, 'units', 'uV', 'scaleFromLims', broadbandScaleFromLims, 'scaleToLims', broadbandScaleToLims, ...
        'dataInMemoryScale', true);
    hasBroadband = cellfun(@(x) size(x, 2) > 0, broadbandData);
    td = td.addOrUpdateBooleanParam(sprintf('has%s', upperFirst(prefix)), hasBroadband, 'mask', tdMask);
end

prog.finish();

td = td.addOrUpdateBooleanParam('nevMerged', mergedMask, 'mask', tdMask);
td = td.addOrUpdateStringParam('nevMergedFile', nevShort,  'mask', tdMask);
