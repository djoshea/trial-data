function trialInfo = getTrialSegmentation(trialSegmentationInfo, spikeData, eventData, nsxData, varargin)
% trialInfo = getTrialSegmentation(trialSegmentationInfo, spikeData, eventData, nsxData)
%   Determines how to parse neural data into trials
%
% trialSegmentationInfo : struct describing how to segment trials
%   .mode : mandatory, determines which segmentation function to use, 
%           e.g. getTrialSegmentation_{mode}.m
% trialInfo : struct array of trial start and stop times (in ms)
%    .startTime (ms)
%    .stopTime (ms)

if isempty(trialSegmentationInfo) || ~isfield(trialSegmentationInfo, 'mode')
    error('No trial segmentation mode specified');
end

mode = trialSegmentationInfo.mode;

% convert the mode to an m file name and call it if 
fnName = sprintf('LoadNev.getTrialSegmentation_%s', mode);
fn = str2func(fnName);

trialInfo = fn(trialSegmentationInfo, spikeData, eventData, nsxData);

if isempty(trialInfo)
    return;
end

% Filter out trials that are too long
if isfield(trialSegmentationInfo, 'maxTrialLength')
    trialDurations = [trialInfo.endTime] - [trialInfo.startTime];
    invalidTrials = trialDurations > trialSegmentationInfo.maxTrialLength;
    if any(invalidTrials)
        fprintf('\tWarning: skipping %d trials which exceed %g ms\n', ...
            nnz(invalidTrials), trialSegmentationInfo.maxTrialLength);
        trialInfo = trialInfo(~invalidTrials);
    end
end
