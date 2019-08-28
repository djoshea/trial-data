function td = resortChannelGroupUsingSpecifiedThreshold(td, channelGroup, varargin)
% extract spiking data from continuous signal by using specific threshold
% or rms value. does not sort the resulting waveforms

p = inputParser();
p.addOptional('threshold', [], @(x) isempty(x) || isscalar(x) || ismatrix(x));
p.addParameter('spikeChPrefix', '', @ischar);
p.addParameter('unitNumber', 0, @isscalar);
p.addParameter('samplesPrePost', [10 38], @isvector); % should add to total number of waveform samples
p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
p.addParameter('mode', 'largestFirst', @ischar);
p.parse(varargin{:});

threshold = p.Results.threshold;

cd = td.getChannelDescriptor(channelGroup);
assert(isa(cd, 'ContinuousNeuralChannelGroupDescriptor'));
array = cd.array;
impl = cd.getImpl();

assert(size(threshold, 2) == cd.nChannels);
if size(threshold, 1) == 1
    threshold = repmat(threshold, td.nTrials, 1);
end
assert(size(threshold, 1) == td.nTrials);

timeDelta = td.getAnalogTimeDelta(channelGroup);
assert(timeDelta > 0);
waveTvec = makecol((-p.Results.samplesPrePost(1) : p.Results.samplesPrePost(2)-1) * timeDelta);
nWaveSamples = sum(p.Results.samplesPrePost(1:2));
assert(nWaveSamples == numel(waveTvec));

%% extract artifact subtracted version for all trials, threshold, and gather waveforms

td = td.reset();

threshold_unscaled = impl.convertAccessDataSingleToMemory(1, threshold);

% fetch the data to be thresholded
[contData, contTime] = td.getAnalogChannelGroup(channelGroup, 'applyScaling', false);

[spikes, waveforms] = ...
    TrialDataUtilities.SpikeData.thresholdExtractSnippetsGroup(contData, contTime, threshold_unscaled, ...
    p.Results.samplesPrePost, 'lockoutPrePost', p.Results.lockoutPrePost, 'mode', p.Results.mode);

%% Add in newly sorted waveforms

if td.hasChannel(array)
    td = td.dropChannel(array);
end

nSpikes = sum(cellfun(@numel, spikes), 1);
for iC= 1:cd.nChannels
    debug('Adding %.1f spikes / trial for electrode %d\n', nSpikes(iC) / td.nTrialsValid, cd.electrodes(iC));
end

units = zeros(cd.nChannels, 1);
td = td.addSpikeArrayChannel(array, cd.electrodes, units, spikes, 'waveforms', waveforms, 'waveformsTime', waveTvec, ...
    'waveformsInMemoryScale', true, 'waveformsScaleFromLims', cd.scaleFromLims, 'waveformsScaleToLims', cd.scaleToLims);

