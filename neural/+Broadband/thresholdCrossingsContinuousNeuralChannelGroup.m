function [td, sortInfo] = thresholdCrossingsContinuousNeuralChannelGroup(td, groupName, varargin)

    p = inputParser();

    p.addParameter('threshold', [], @isscalar);
    p.addParameter('method', '', @(x) any(validatestring(x, {'', 'rmsThreshold', 'setThreshold'})));
    p.addParameter('rmsMultiplier', -4.5, @isscalar);
    p.addParameter('rmsPerTrial', false, @islogical);
    p.addParameter('smoothOverTrials', 0, @isscalar);
    p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
    p.addParameter('mode', 'largestFirst', @ischar);

    p.parse(varargin{:});

    td = td.reset();

    sortMethod = p.Results.method;
    if isempty(sortMethod)
        if ~isempty(p.Results.threshold)
            sortMethod = 'setThreshold';
        else
            sortMethod = 'rmsThreshold';
        end
    end

    if ismember(sortMethod, {'rmsThreshold'})
        % get rms threshold
        debug('Computing RMS per channel for group %s\n', groupName);
        rms1 = Broadband.computeGroupRMSThresholdByChannel(td, groupName, 1, ...
            'perTrial', p.Results.rmsPerTrial, 'smoothOverTrials', p.Results.smoothOverTrials);
        rmsThresh = rms1 * p.Results.rmsMultiplier;
    end

    setThresh = p.Results.threshold;

    % take artifact removed, HP-filtered signal and extract spike waveforms
    switch sortMethod
        case 'rmsThreshold'
            thresh = rmsThresh;
        case 'setThreshold'
            thresh = setThresh;
        otherwise
            error('Unknown sortmethod');
    end

    td = Broadband.resortChannelGroupUsingSpecifiedThreshold(td, ...
        groupName, thresh, 'mode', p.Results.mode, 'lockoutPrePost', p.Results.lockoutPrePost);

    sortInfo.thresh = thresh;
    sortInfo.sortMode = p.Results.mode;
    sortInfo.rmsMultiplier = p.Results.rmsMultiplier;


end