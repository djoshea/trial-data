function [tdca, sortInfo] = rethresholdChannel(tdca, chName, varargin)

p = inputParser();

p.addParameter('threshold', [], @isscalar);
p.addParameter('method', '', @(x) any(validatestring(x, {'', 'rmsThreshold', 'preserveThreshold', 'setThreshold', 'maxRmsPreserveSet'})));
p.addParameter('rmsMultiplier', -4.5, @isscalar);
p.addParameter('rmsPerTrial', false, @islogical);
p.addParameter('smoothOverTrials', 0, @isscalar);
p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
p.addParameter('mode', 'largestFirst', @ischar);

p.parse(varargin{:});

tdca = tdca.reset();

sortMethod = p.Results.method;
if isempty(sortMethod)
    if ~isempty(p.Results.threshold)
        sortMethod = 'setThreshold';
    else
        sortMethod = 'rmsThreshold';
    end
end

if ismember(sortMethod, {'rmsThreshold', 'maxRmsPreserveSet'})
    % get rms threshold
    rms1 = Broadband.computeRMSThreshold(tdca, chName, 1, ...
        'perTrial', p.Results.rmsPerTrial, 'smoothOverTrials', p.Results.smoothOverTrials);
    rmsThresh = rms1 * p.Results.rmsMultiplier;
%     debug('Mean RMS is %.0f uV, %g x RMS is %.0f uV\n', nanmean(rms1), p.Results.rmsMultiplier, nanmean(rmsThresh));
end

if ismember(sortMethod, {'preserveThreshold', 'maxRmsPreserveSet'})
    % get preserve threshold
    preserveThresh = TrialDataUtilities.SpikeData.estimateThresholdFromSpikeWaveforms(tdca, chName);
    debug('Preserve thresh is %.0f uV\n', preserveThresh)
end

setThresh = p.Results.threshold;

if strcmp(sortMethod, 'maxRmsPreserveSet')
    if p.Results.rmsPerTrial
        % one value of rms per trial
        preserveThresh = repmat(preserveThresh, size(rmsThresh));
        setThresh = repmat(setThresh, size(rmsThresh));
    end
    
    values = [rmsThresh, preserveThresh, setThresh];
    [~, which] = nanmax(abs(values), 2);
    maxThresh = TensorUtils.selectSpecificIndicesAlongDimensionEachPosition(values, 2, which);
%     list = {'rms', 'preserve', 'set'};
%     whichMax = list{which};
end

% take artifact removed, HP-filtered signal and extract spike waveforms
switch sortMethod
    case 'rmsThreshold'
        thresh = rmsThresh;
    case 'setThreshold'
        thresh = setThresh;
    case 'preserveThreshold'
        thresh = preserveThresh;
    case 'maxRmsPreserveSet'
        thresh = maxThresh;
    otherwise
        error('Unknown sortmethod');
end

% if strcmp('maxRmsPreserveSet', sortMethod)
%     debug('Using %s threshold (%s) of %.0f uV\n', sortMethod, whichMax, thresh);
% else
%     debug('Using %s threshold of %.0f uV\n', sortMethod, thresh);
% end
[tdca, sortInfo.spikeCh] = Broadband.resortChannelUsingSpecifiedThreshold(tdca, ...
    chName, thresh, 'mode', p.Results.mode, 'lockoutPrePost', p.Results.lockoutPrePost, ...
    'thresholdPerTrial', numel(thresh) > 1);

sortInfo.thresh = thresh;
sortInfo.sortMode = thresh;
sortInfo.rmsMultiplier = p.Results.rmsMultiplier;

end
