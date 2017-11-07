function trialInfo = getTrialSegmentation_StartChannel(trialSegmentationInfo, spikeData, eventData, nsxData, varargin)
% This takes the analog info from the cerebus file and 
% decodes the trial ID code and returns the time of the
% start of each trial

%trialSegmentation.mode = 'StartChannel';
%trialSegmentation.startCh = 137;

[startTimes, startIdx, tvec] = findRisingEdges(trialSegmentationInfo.startCh, nsxData);

if isempty(startTimes)
	trialInfo = [];
	return;
end

% for each start time, the subsequent pulse defines the end time of the
% previous trial
endTimes = [tvec(startIdx(2:end)-1), tvec(end)];

% build trialInfo struct with start/end times 
for it = 1:length(startTimes)
    trialInfo(it).startTime = startTimes(it);
	trialInfo(it).endTime = endTimes(it);
end

end

function [times, riseEdgeIdx, tvec] = findRisingEdges(chId, nsxData)
	
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
    
    tvec = time;
    
end

