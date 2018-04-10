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

p.parse(varargin{:});
includeSpikes = p.Results.overwriteSpikes;
includeWaveforms = p.Results.includeWaveforms;

empty = [];
array = p.Results.array;

% align td to Q based on computed delay periods
qMatchInTD = LoadNev.findNevTrialsMatchInTrialData(td, Q);

prog = ProgressBar(numel(Q), 'Merging nev trials into TrialData');
nTD = td.nTrials;
spikeData = struct();
waveData = struct();

td = td.reset();
nevShort = cellvec(nTD);

waveformScaleBy = Q(1).CerebusInfo.waveformScaleBy;

convertName = @(dotName) strrep(dotName, '.', '_');

mergedMask = falsevec(td.nTrials);

nsxGroupData = struct();

for iQ = 1:numel(Q)
    prog.update(iQ);
    if isnan(qMatchInTD(iQ))
        continue;
    end
    
    q = Q(iQ);
    iR = qMatchInTD(iQ);
    mergedMask(iR) = true;
    if isempty(q.spikes)
        units = {};
    else
        units = fieldnames(q.spikes);
    end
    
    tOffsetQ = 0;

    if ~isempty(nevShort{iR})
        warning('Trial %d in TrialData has already been merged with file %s', iR, nevShort{iR});
    end 
    nevShort{iR} = q.CerebusInfo.nevNameShort;
    
    % grab spikes into spikeData struct
    for iK = 1:numel(units)
        chName = convertName(units{iK});
        if ~isfield(spikeData, chName) % add field if needed
            spikeData.(chName) = cellvec(nTD);
            waveData.(chName) = cellvec(nTD);
        end
            
        spikeData.(chName){iR} = q.spikes.(units{iK}) - tOffsetQ;
        if includeWaveforms
            waveData.(chName){iR} = q.waves.(units{iK});
        end
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
    units = fieldnames(spikeData);
    
    prog = ProgressBar(numel(units), 'Adding spike data to TrialData');
    for iU = 1:numel(units)
        chName = strrep(units{iU}, 'unit', array);
        prog.update(iU, 'Adding %s to TrialData', chName);
        if includeWaveforms
            wd = waveData.(units{iU});
            nonEmpty = find(~cellfun(@isempty, wd), 1);
            if isempty(nonEmpty)
                wd = [];
                waveClass = 'int16';
            else
                % swap dims of each waveform array (samples along dim 2)
                wd = cellfun(@(x) x', wd, 'UniformOutput', false);
                nWave = size(wd{nonEmpty}, 2);
                wavetvec = (-10 : nWave-11) / 30;
                waveClass = class(wd{nonEmpty(1)});
            end
            
            waveLimsFrom = double([intmin(waveClass) intmax(waveClass)]);
            waveLimsTo = waveLimsFrom * waveformScaleBy;
            
            td = td.addOrUpdateSpikeChannel(chName, spikeData.(units{iU}), 'mask', tdMask, ...
                'waveforms', wd, 'waveformsTime', wavetvec, ...
                'waveformsScaleFromLims', waveLimsFrom, ...
                'waveformsScaleToLims', waveLimsTo);
            
            waveData.(units{iU}) = [];
        else
            td = td.addOrUpdateSpikeChannel(chName, spikeData.(units{iU}), 'mask', tdMask);
        end
        spikeData.(units{iU}) = [];
    end
    prog.finish();
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
        td = td.addOrUpdateContinuousNeuralChannelGroup(grpName, grp.names, ...
             grp.data, grp.time, 'mask', tdMask, 'units', grp.units, 'scaleFromLims', grp.scaleFromLims, 'scaleToLims', grp.scaleToLims, ...
             'dataInMemoryScale', true);
    else
        td = td.addOrUpdateAnalog(grpName, grp.data, grp.time, 'mask', tdMask, 'units', grp.units, ...
            'scaleFromLims', grp.scaleFromLims, 'scaleToLims', grp.scaleToLims, ...
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
