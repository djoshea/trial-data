nTrials = 10;
nUnits = 5;
T = 100;
arrayName = 'utahA';

td = TrialData.buildEmptyWithTrialDurations(repmat(T, nTrials, 1));

electrodes = 1:2:(nUnits*2)-1;
units = 2:2:nUnits*2;
spikes = cell(nTrials, nUnits);
for i = 1:numel(spikes)
    spikes{i} = randi(100, 20);
end

td = td.addSpikeChannelArray(arrayName, electrodes, units, spikes);