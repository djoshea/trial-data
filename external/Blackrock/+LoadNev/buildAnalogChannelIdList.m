function channelIdList = buildAnalogChannelIdList(analogInfo, varargin)
% channelIdList = buildAnalogChannelIdList(analogInfo, varargin)
% Based on the description of the analog channels requested, builds a simple list of all channel ids
% requested to pass into nevExtractAnalog()
%
% analogInfo : 2-deep struct array, 
%   	analogInfo.channelGroupName.channelName = Channel_ID
%   		results in Q(i).channelGroupName.channelName = [ parsed data ]
%   				   Q(i).channelGroupName.time = [ associated time vector ]
%   	analogInfo.channelGroupName = [ Channel_ID List ]
%   		results in Q(i).channelGroupName.data = [data matrix]
%   				   Q(i).channelGroupName.
%   				   Q(i).channelGroupName.time = [ time vector]
%       all channels within a channelGroup must come from the same nsxFile (same element of nsxData)
%
% channelIdList: a simple scalar array of Channel_IDs

channelIdList = [];

channelGroups = fieldnames(analogInfo);
nChannelGroups = length(channelGroups);
for icg = 1:nChannelGroups
    groupName = channelGroups{icg};
    group = analogInfo.(groupName);
    if isempty(group)
        continue;
    end
    
    if isstruct(group)
        % a group of single channel named lookups
        names = fieldnames(group);
       
        channelIdsThisGroup = cellfun(@(name) group.(name), names);
    else
        % multiple channel lookup to be squashed into one array
        channelIdsThisGroup = group;
    end

    channelIdList = [channelIdList; makecol(channelIdsThisGroup)];
end

end
