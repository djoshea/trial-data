function [tdca, spikeChList] = resortChannelUsingSpecifiedThreshold(tdca, chNameAR, varargin)
% extract spiking data from continuous signal by using specific threshold
% or rms value. does not sort the resulting waveforms

p = inputParser();
p.addOptional('threshold', [], @(x) isempty(x) || isscalar(x) || isvector(x));
p.addParameter('thresholdPerTrial', false, @islogical);
p.addParameter('spikeChPrefix', '', @ischar);
p.addParameter('unitNumber', 0, @isscalar);
p.addParameter('samplesPrePost', [10 38], @isvector); % should add to total number of waveform samples
p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
p.addParameter('mode', 'largestFirst', @ischar);
p.parse(varargin{:});

threshold = p.Results.threshold;

if p.Results.thresholdPerTrial
    assert(numel(threshold) == tdca.nTrials);
end

timeDelta = tdca.getAnalogTimeDelta(chNameAR);
waveTvec = makecol((-p.Results.samplesPrePost(1) : p.Results.samplesPrePost(2)-1) * timeDelta);
nWaveSamples = sum(p.Results.samplesPrePost(1:2));
assert(nWaveSamples == numel(waveTvec));

%% extract artifact subtracted version for all trials, threshold, and gather waveforms

tdca = tdca.reset();

% fetch the data to be thresholded
[contData, contTime] = tdca.getAnalog(chNameAR);

[timeCellByTrial, waveformsByTrial] = ...
    TrialDataUtilities.SpikeData.thresholdExtractSnippets(contData, contTime, threshold, ...
    p.Results.samplesPrePost, 'lockoutPrePost', p.Results.lockoutPrePost, 'mode', p.Results.mode, ...
    'thresholdPerTrial', p.Results.thresholdPerTrial);

%% Add in newly sorted waveforms

cd = tdca.getChannelDescriptor(chNameAR);

oldSpikeCh = tdca.listSpikeChannelsOnArrayElectrode(sprintf('%s%s', p.Results.spikeChPrefix, cd.array), cd.electrode);
tdca = tdca.dropChannel(oldSpikeCh);
spikeCh = sprintf('%s%s%02d_%d', p.Results.spikeChPrefix, cd.array, cd.electrode, p.Results.unitNumber);

nSpikes = sum(cellfun(@numel, timeCellByTrial));
debug('Adding %.1f spikes / trial for channel %s\n', nSpikes / tdca.nTrialsValid, spikeCh)
tdca = tdca.addSpikeChannel(spikeCh, timeCellByTrial, ...
    'waveforms', waveformsByTrial, ...
    'waveformsTime', waveTvec);

spikeChList = {spikeCh};

