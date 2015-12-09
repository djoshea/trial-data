function [spikeTimes, waveforms] = continuousHPFilterAndThreshold(...
    dataTimeByTrials, spikeThreshold)

% Design the appropriate high pass filter 
Fs = 30000;
hpCornerHz = 250; 
hpCornerNormalized = hpCornerHz / (Fs/2);
[B, A] = butter(4, hpCornerNormalized, 'high');
dataHP = TrialDataUtilities.Data.filterIgnoreLeadingTrailingNaNs(B, A, dataTimeByTrials);



end



