% Load the MatUdp data from disk and build trial data 
td = MatUdp.DataLoadEnv.loadSaveTagTrialData('protocol', protocol, 'subject', subj, ...
    'dateStr', dateStr, 'saveTag', saveTag, 'includeSpikeData', true, ...
    'includeWaveforms', true, 'includeContinuousNeuralData', true);
td = TrialDataConditionAlign(td);

% customize datasetName if desired
td.datasetName = sprintf('%s %s %s', subj, dateStr, protocol);

% drop online captured spike channels
list = td.listSpikeChannels();
td = td.dropChannels(list);

% load all nevs in directory into cell array
[Qcell, ~, fileList] = LoadNev.loadNevsInDirectory(nevPath, @LFADS_PierreExport.Nev.loadNev_ArmDynamics, ...
    'excludeFilenameContains', 'sort');

% merge nevs into trial data one by one
numMerged = nanvec(numel(Qcell));
for iQ = 1:numel(Qcell)
    [td, numMerged(iQ)] = LoadNev.mergeNevTrialsIntoTrialData(td, Qcell{iQ});
    debug('Merged %d trials from nev file %s into trialData\n', numMerged(iQ), fileList{iQ});
end

% mergeNevTrialsIntoTrialData will add this flag channel if broadband data
% are found
if ~td.hasParamChannel('hasBroadband')
    error('No broadband data found');
end

% if you want to drop channels with no neural data, include this
hasBroadband = td.getParam('hasBroadband');
td = td.markTrialsPermanentlyInvalid(~hasBroadband, 'Missing broadband data from nev files');

% remap channels according to Cerebus Map File
mapFile = 'v24.cmp';
mapDir = getenvCheckPath('CEREBUS_MAP_PATH');
mapFileFull = fullfile(mapDir, mapFile);
debug('Remapping channels using map file %s\n', mapFile);
td = remapNeuralChannelsUsingCerebusMapFile(td, mapFileFull, 'contiguousContinuousNumbering', false);

% sort broadband channel order according to their new channel numbers
td = td.sortChannelsInAnalogChannelGroup('broadband');

% high pass filter at 250 Hz
td = Broadband.highPassFilterAllChannelsInGroupInPlace(td, 'broadband', 'hpCornerHz', 250);

% rethreshold channels at -3.5x RMS per trial, smoothed over 10 trials
td = td.dropChannels(td.listSpikeChannels);
td = Broadband.rethresholdAllChannelsInGroup(td, 'broadband', ...
    'method', 'rmsThreshold', 'rmsMultiplier', -3.5, ...
    'rmsPerTrial', true, 'smoothOverTrials', 10);

td = td.dropAnalogChannelGroup('broadband');

% remove shared spike artifacts across multiple channels
[td, artifactCounts] = TrialDataUtilities.SpikeData.removeSharedSpikeArtifacts(...
    td, 'minChannels', 6, 'timeWindow', 0.1);
fprintf('Max artifacts removed %d\n', max(artifactCounts(:)));

