function [tdca, sortInfo] = rethresholdAllChannelsInGroup(tdca, groupName, varargin)

    chList = tdca.listAnalogChannelsInGroup(groupName);
    sortInfo = cellvec(numel(chList));

    for c = 1:numel(chList)
%         debug('Rethreholding channel %s\n', chList{c});
        [tdca, sortInfo{c}] = Broadband.rethresholdChannel(tdca, chList{c}, varargin{:});

        if isempty(sortInfo{c}.spikeCh)
            continue;
        end
    end

end