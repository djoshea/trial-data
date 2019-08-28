function [td, numMerged, empty] = mergeNevTrialsIntoTrialData(td, Q, varargin)
% given a Q struct returned by loadNevNoSerial and a ChannelReach trial data instance,
% merge the spikes from the nev into the trial data with the same unit
% names. spikes from the td will not be deleted, but they will be
% overwritten on merged trials

p = inputParser;
p.addParameter('includePhotobox', false, @islogical); % copy Q photobox to td.photoboxNev
p.addParameter('overwriteSpikes', true, @islogical);
p.addParameter('includeWaveforms', true, @islogical);
p.addParameter('dropChannelsWithSameName', false, @islogical);
p.addParameter('array', 'ch', @ischar); % prefix spike channels with this name
p.addParameter('neuralChannelGroupList', {'lfp', 'broadband'}, @iscellstr);

p.parse(varargin{:});
includeSpikes = p.Results.overwriteSpikes;
includeWaveforms = p.Results.includeWaveforms;

empty = [];
array = p.Results.array;

% align td to Q based on computed delay periods
qMatchInTD = LoadNev.findNevTrialsMatchInTrialData(td, Q);

prog = ProgressBar(numel(Q), 'Merging nev trials into TrialData');
nTD = td.nTrials;

td = td.reset();
nevShort = cellvec(nTD);

waveformScaleBy = Q(1).CerebusInfo.waveformScaleBy;

mergedMask = falsevec(td.nTrials);
nsxGroupData = struct();
units =[];
electrodes = [];

for iQ = 1:numel(Q)
    prog.update(iQ);
    if isnan(qMatchInTD(iQ))
        continue;
    end
    
    q = Q(iQ);
    iR = qMatchInTD(iQ);
    mergedMask(iR) = true;
    
    tOffsetQ = 0;
    if isempty(units) && ~isempty(Q(iQ).units)
        units = Q(iQ).units;
        electrodes = Q(iQ).electrodes;
    end

    if ~isempty(nevShort{iR})
        warning('Trial %d in TrialData has already been merged with file %s', iR, nevShort{iR});
    end 
    nevShort{iR} = q.CerebusInfo.nevNameShort;
    
    if tOffsetQ ~= 0
        Q(iQ).spikes = cellfun(@(v) v - tOffsetQ, Q(iQ).spikes, 'UniformOutput', false);
    end
   
    % handle nsx analog signals
    if isfield(q, 'nsxData')
        nsxGroupsThisTrial = fieldnames(q.nsxData);
    else
        nsxGroupsThisTrial = {};
    end
    for iG = 1:numel(nsxGroupsThisTrial)
        grp = nsxGroupsThisTrial{iG};
        if ~isfield(nsxGroupData, grp)
            nsxGroupData.(grp).isGroup = q.nsxData.(grp).isGroup;
            nsxGroupData.(grp).scaleFromLims = q.nsxData.(grp).scaleLims(1:2);
            nsxGroupData.(grp).scaleToLims  = q.nsxData.(grp).scaleLims(3:4);
            nsxGroupData.(grp).data = cell(nTD, 1);
            nsxGroupData.(grp).time = cell(nTD, 1);
            nsxGroupData.(grp).names = q.nsxData.(grp).names;
            nsxGroupData.(grp).units = q.nsxData.(grp).units;
            nsxGroupData.(grp).isNeural = ismember(grp, p.Results.neuralChannelGroupList);
            nsxGroupData.(grp).channelIds = q.nsxData.(grp).channelIds;
        end
        
        % grab this trial's data
        nsxGroupData.(grp).data{iR} = q.nsxData.(grp).data';
        nsxGroupData.(grp).time{iR} = makecol(q.nsxData.(grp).time - tOffsetQ);
    end
end
prog.finish();

numMerged = nnz(mergedMask);
if numMerged == 0
    return
end

tdMask = removenan(qMatchInTD);

if includeSpikes
    units = Q(1).units;
    
    
    nUnits = numel(units);
    qSpikes = cat(1, Q.spikes);
    spikes = cell(nTD, nUnits);
    maskMatched = ~isnan(qMatchInTD);
    spikes(qMatchInTD(maskMatched)) = qSpikes(maskMatched);
    
    waveforms = [];
    waveLimsFrom = [];
    waveLimsTo = [];
    wavetvec = [];
    
    if includeWaveforms
        qWaveforms = cat(1, Q.waves);
        waveforms = cell(nTD, nUnits);
        waveforms(qMatchInTD(maskMatched)) = qWaveforms(maskMatched);
        
        nonEmpty = find(~cellfun(@isempty, waveforms), 1);
        if ~isempty(nonEmpty)
            % swap dims of each waveform array (spikes along dim 1, samples along dim 2)
            waveforms = cellfun(@(x) x', waveforms, 'UniformOutput', false);
            exampleWaves = waveforms{nonEmpty};
            nWave = size(exampleWaves, 2);
            wavetvec = (-10 : nWave-11) / 30; % hardcoded blackrock offset
            waveClass = class(exampleWaves);
            
            waveLimsFrom = double([intmin(waveClass) intmax(waveClass)]);
            waveLimsTo = waveLimsFrom * waveformScaleBy;
        end
    end
    
    if isfield(Q, 'sortRatings')
        qSortQuality = cat(1, Q.sortQuality);
        sortQuality = cell(nTD, nUnits);
        sortQuality(qMatchInTD(maskMatched)) = qSortQuality(maskMatched);
    else
        sortQuality = [];
    end
    
    td = td.addSpikeArrayChannel(array, electrodes, units, spikes, 'isAligned', false, 'waveforms', waveforms, ...
        'waveformsTime', wavetvec, 'waveformsScaleFromLims', waveLimsFrom, 'waveformsScaleToLims', waveLimsTo, ...
        'waveformsInMemoryScale', true, 'sortQualityEachTrial', sortQuality);
end

% TODO add handling for multiple arrays by renaming the signal
nsxGroupNames = fieldnames(nsxGroupData);
for iG = 1:numel(nsxGroupNames)
    grpName = nsxGroupNames{iG};
    grp = nsxGroupData.(grpName);
    
    if p.Results.dropChannelsWithSameName
        td = td.dropChannels(grpName).dropChannels(grp.names);
    end
    
    if grp.isGroup
        debug('Adding %s data to TrialData as continuous neural group with %d channels\n', grpName, numel(grp.names));
        if grp.isNeural
            electrodes = grp.channelIds;
            td = td.addOrUpdateContinuousNeuralChannelGroup(append(grpName, '_', array), electrodes, ...
                 grp.data, grp.time, 'mask', tdMask, 'units', grp.units, 'scaleFromLims', grp.scaleFromLims, 'scaleToLims', grp.scaleToLims, ...
                 'subChannelNames', grp.names, ...
                 'dataInMemoryScale', true);
        else
            if ischar(grp.units)
                subChannelUnits = [];
                units = grp.units;
            else
                subChannelUnits = grp.units;
                units = '';
            end
            td = td.addOrUpdateAnalogChannelGroup(grpName, ...
                 grp.data, grp.time, 'mask', tdMask, 'units', units, 'scaleFromLims', grp.scaleFromLims, 'scaleToLims', grp.scaleToLims, ...
                 'subChannelNames', grp.names, 'subChannelUnits', subChannelUnits, ...
                 'dataInMemoryScale', true);
        end
    else
        td = td.addOrUpdateAnalog(grpName, grp.data, grp.time, 'mask', tdMask, 'units', grp.units, ...
            'scaleFromLims', grp.scaleFromLims, 'scaleToLims', grp.scaleToLims, ...
            'subChannelNames', grp.names, ...
            'dataInMemoryScale', true);
    end    
    nsxGroupData.(grpName) = [];
end

prog.finish();

if p.Results.dropChannelsWithSameName
    td = td.dropChannels({'nevMerged', 'nevMergedFile'});
end
td = td.addOrUpdateBooleanParam('nevMerged', mergedMask, 'mask', tdMask);
td = td.addOrUpdateStringParam('nevMergedFile', nevShort,  'mask', tdMask);
