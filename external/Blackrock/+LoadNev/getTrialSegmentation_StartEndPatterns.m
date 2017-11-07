function trialInfo = getTrialSegmentation_StartEndPatterns(trialSegmentationInfo, spikeData, eventData, nsxData, varargin)
% This takes the event stream from the cerebus file and 
% decodes the trial ID code and returns the time of the
% start and stop signal received from TEMPO
% 
% Start and end times are defined by a pattern of successive pulse values (all 16 bits) 
% which must match exactly, except NaNs may be used as a wildcard.
%
% The startPattern is scanned against the event codes to find exact matches,
% nans are considered to be wildcards that match anything
% The trial actually begins at the first element of this list
%
% The endPattern is scanned similarly, the trial actually ends
% at the last element of the list

% example:
%  startPattern = [32768 NaN NaN NaN NaN 32769];
%  endPattern = 0;


if ~isfield(trialSegmentationInfo, 'startPattern') || ~isfield(trialSegmentationInfo, 'endPattern')
    error('Must specify both startPattern and endPatterns');
end

startPattern = trialSegmentationInfo.startPattern;
endPattern = trialSegmentationInfo.endPattern;

if isequal(startPattern, endPattern)
    error('Identical start and end patterns are not supported');
end

% add a timestamp at the very end of the recording
lastTimeNSx = nanmax(arrayfun(@(nsx) nsx.time(end), nsxData));
lastTimeSpike = spikeData.timestamp(end);
if isempty(lastTimeNSx)
    lastTime = lastTimeSpike;
else    
    lastTime = nanmax(lastTimeNSx, lastTimeSpike);
end
eventData.timestamp(end+1) = lastTime;

% Find pattern matches
startInds = matchPatternAgainstStream(startPattern, eventData.code);
if isempty(endPattern)
    endInds = [startInds(2:end) numel(eventData.timestamp)];
else
    endInds = matchPatternAgainstStream(endPattern, eventData.code);
end

if isempty(startInds) || isempty(endInds)
	trialInfo = [];
	return;
end

startInds = startInds(startInds < endInds(end));
endInds = endInds(endInds > startInds(1));

if isempty(startInds) || isempty(endInds)
	trialInfo = [];
	return;
end

% for each start ind, find the most proximal subsequent end code
nextEndIndFn = @(startInd) endInds(find(endInds > startInd, 1, 'first'));
endInds = arrayfun(nextEndIndFn, startInds);

% now for each end ind, find the most proximal previous start code
prevStartIndFn = @(endInd) startInds(find(startInds < endInd, 1, 'last'));
startInds = arrayfun(prevStartIndFn, endInds);

if isempty(startInds) || isempty(endInds)
	trialInfo = [];
	return;
end

% at this point we should have an equal number of start and end codes
assert(numel(startInds) == numel(endInds), 'Unequal number of start and end codes');

% for the sake of odd bugs cropping up in the segmentation, 
% check that every time stamp is unique
% thus allowing us to go back and parse by time rather than ind
if numel(unique(eventData.timestamp)) ~= numel(eventData.timestamp)
    fprintf(2, '\t\tWarning: %d non-unique timestamps found in nev file\n', ...
        numel(eventData.timestamp) - numel(unique(eventData.timestamp)));
end

% build trialInfo struct with start/end times 
for it = 1:length(startInds)
    trialInfo(it).startTime = eventData.timestamp(startInds(it));
	trialInfo(it).endTime = eventData.timestamp(endInds(it));
end

end

function matchInds = matchPatternAgainstStream(pattern, stream)
    % start with all potential inds
    matchInds = 1:(length(stream)-length(pattern)+1);

    % loop over the elements of the pattern and keep only matchInds that match
    % up to that point
    for ip = 1:length(pattern)
        if isnan(pattern(ip))
            continue;
        end
    
        matchInds = matchInds(stream(matchInds+ip-1) == pattern(ip));
    end
end
