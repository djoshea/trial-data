function trialInfo = getTrialSegmentation_StartEndChannels(trialSegmentationInfo, spikeData, eventData, nsxData, varargin)
% This takes the event stream from the cerebus file and 
% decodes the trial ID code and returns the time of the
% start and stop signal received from TEMPO

%trialSegmentation.mode = 'StartEndChannels';
%trialSegmentation.startCh = 137;
%trialSegmentation.stopCh = 138;

startTimes = findRisingEdges(trialSegmentationInfo.startCh, nsxData);
endTimes   = findRisingEdges(trialSegmentationInfo.endCh, nsxData);

if isempty(startTimes) || isempty(endTimes)
	trialInfo = [];
	return;
end

% remove outliers
startTimes = startTimes(startTimes < endTimes(end));
endTimes   = endTimes(endTimes > startTimes(1));

if isempty(startTimes) || isempty(endTimes)
	trialInfo = [];
	return;
end

% for each start ind, find the most proximal subsequent end code
nextEndTimeFn = @(startTime) endTimes(find(endTimes > startTime, 1, 'first'));
endTimes = arrayfun(nextEndTimeFn, startTimes);

% now for each end ind, find the most proximal previous start code
prevStartTimeFn = @(endTime) startTimes(find(startTimes < endTime, 1, 'last'));
startTimes = arrayfun(prevStartTimeFn, endTimes);

if isempty(startTimes) || isempty(endTimes)
	trialInfo = [];
	return;
end

% at this point we should have an equal number of start and end codes
assert(numel(startTimes) == numel(endTimes), 'Unequal number of start and end codes');

% build trialInfo struct with start/end times 
for it = 1:length(startTimes)
    trialInfo(it).startTime = startTimes(it);
	trialInfo(it).endTime = endTimes(it);
end

end

function times = findRisingEdges(chId, nsxData)
	
	% find the channel by its id in nsxData
	for insx = 1:length(nsxData)
		channelIds = nsxData(insx).channelIds;

		ind = find(channelIds == chId, 1);
		if ~isempty(ind)
			data = nsxData(insx).data(ind, :);
			time = nsxData(insx).time;
			break;
		end
	end	

	if isempty(nsxData) || isempty(ind) || isempty(data)
		error('Could not find channel %d', chId);
	end

	% locate the rising pulse edges
	edgeIdx = find(edge(data) == 1);
	riseEdgeIdx = edgeIdx(1:2:end);

	% and look up the times at which these occur
	times = time(riseEdgeIdx);
end

