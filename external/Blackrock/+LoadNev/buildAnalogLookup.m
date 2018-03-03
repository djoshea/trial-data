function analogLookup = buildAnalogLookup(analogInfo, nsxData, varargin)
% lookup the indices of the analog channels and build lookup table to used by
% addRequestedAnalogData.m
%
% analogInfo : 2-deep struct array, 
%       analogInfo.channelName = Channel_ID
%           results in Q(i).nsxSingle.channelName.data = [ parsed data]
%                      Q(i).nsxSingle.channelName.time = [ associated time vector ]
%   	analogInfo.channelGroupName.channelName = Channel_ID
%       analogInfo.channelGroupName = [ Channel_ID List ]
%   		results in Q(i).nsxGroup.channelGroupName.data = [ parsed data matrix ]
%   				   Q(i).nsxGroup.channelGroupName.time = [ associated time vector ]
%
%       all channels within a channelGroup must come from the same nsxFile (same element of nsxData)
%
%   also the field meta_nsxExtMaskByGroup.groupName = {'.ns5'} allows you to specify
%   that a given group must be found within the nsx with that extension
% 
% analogLookup(i) has fields
%    .groupName --> name of signal group or '' for signals
%    .names --> cellstr of signal name(s)
%    .nsxIndex: which nsx this comes from 
%    .indWithinNsx : nsxData(nsxIndex).data(chInd, :) are these channels' data
%    .single : false for groups, true for single
%    .lookup : for mutliple channels only, the lookup table corresponding to the elements of chInd
%               i.e. data(chInd(i), :) corresponds to analogInfo.group( .lookup (i) )
%              this is necessary because not all channels in the analogInfo.group array may be found
%    .scaleFns, .scaleLims

    p = inputParser();
%     p.addParameter('suppressNotFoundWarningsForGroups', {'lfp', 'broadband'}, @iscellstr);
    p.addParameter('keepExtraAnalogChannels', true, @islogical); % if false, channels not mentioned in analogInfo will be discarded
    p.parse(varargin{:});

%     suppressNotFoundWarningsForGroups = p.Results.suppressNotFoundWarningsForGroups;
    keepExtraAnalogChannels = p.Results.keepExtraAnalogChannels;

    % keep track of which channels are used, so that we know which are extra
    nsxDataChannelIsExtraMask = arrayfun(@(nsx) true(size(nsx.data, 1), 1), nsxData, 'UniformOutput', false);

    fprintf('\tBuilding analog channel lookup table\n');

    analogLookup = [];

    if isempty(analogInfo)
        analogInfo = struct();
    end
    
    if isfield(analogInfo, 'meta_nsxExtMaskByGroup')
        nsxExtByGroup = analogInfo.meta_nsxExtMaskByGroup;
        if isempty(nsxExtByGroup)
            nsxExtByGroup = struct();
        end
        analogInfo = rmfield(analogInfo, 'meta_nsxExtMaskByGroup');
    else
        nsxExtByGroup = struct();
    end

    channelGroups = fieldnames(analogInfo);
    nChannelGroups = length(channelGroups);
    for icg = 1:nChannelGroups
        groupName = channelGroups{icg};
        group = analogInfo.(groupName);
        if isempty(group)
            continue;
        end
        
        if isfield(nsxExtByGroup, groupName)
            extMatches = @(ext) ismember(ext, nsxExtByGroup.(groupName)) || ismember(ext(2:end), nsxExtByGroup.(groupName));
            nsxMaskThisGroup = arrayfun(@(nsx) extMatches(nsx.ext), nsxData);
        else
            nsxMaskThisGroup = true(numel(nsxData), 1);
        end
        if ~any(nsxMaskThisGroup), continue; end

        if isscalar(group) && isnumeric(group)
            id = group;
            signalName = groupName;
            % single named channel
            
            [found, nsxIndex, idxWithinNsx, nsxDataChannelIsExtraMask] = findChannelIdInNsxData(nsxData, id, nsxMaskThisGroup, nsxDataChannelIsExtraMask);

            if found
                
                lookup.groupName = '';
                lookup.names = {signalName};
                lookup.nsxIndex = nsxIndex;
                lookup.idxWithinNsx = idxWithinNsx;
                lookup.channelIds = nsxData(nsxIndex).channelIds(idxWithinNsx);
                lookup.single = true;
                lookup.lookup = 1;
                lookup.scaleFn = nsxData(nsxIndex).scaleFns{idxWithinNsx};
                lookup.scaleLims = nsxData(nsxIndex).scaleLims{idxWithinNsx};
                lookup.units = getCheckUnits(signalName, nsxData, nsxIndex, idxWithinNsx);
            else
                warning('Could not find NSx channel %s', signalName);
                continue;
            end
            
        elseif isstruct(group) || (isvector(group) && isnumeric(group))
            % named individual channels if struct
            % e.g. analogInfo.groupName.channelName1 = 1;
            %      analogInfo.groupName.channelName2 = 2;
            % a range of channel ids if vector
            % e.g. analogInfo.groupName = 1:10;
            
            if isstruct(group)
                names = fieldnames(group);
                ids = cellfun(@(name) group.(name), names);
            else
                ids = group;
                names = arrayfun(@(ind) sprintf('%s%03d', groupName, ind), 1:numel(ids), 'UniformOutput', false);
            end
            
            found = false(numel(names), 1);
            [nsxIndex, idxWithinNsx] = deal(nan(numel(names), 1));

            % loop over analog channels in this group
            for ic = 1:length(names)
                [found(ic), nsxIndex(ic), idxWithinNsx(ic), nsxDataChannelIsExtraMask] = findChannelIdInNsxData(nsxData, ids(ic), nsxMaskThisGroup, nsxDataChannelIsExtraMask);
            end
            
            if ~any(found)
                warning('Could not find any NSx channels in group %s', groupName);
                continue;
            else
                ids = ids(found);
                nsxIndex = nsxIndex(found);
                idxWithinNsx = idxWithinNsx(found);
                names = names(found);

                assert(numel(unique(nsxIndex)) == 1, 'Channels in group %s were spread over multiple nsx with different sampling rates', groupName);
                nsxIndex = nsxIndex(1);
                
                scale = nsxData(nsxIndex).scaleLims{1};
                if numel(nsxData(nsxIndex).scaleLims) > 1
                    for j = 2:numel(ids)
                        assert(isequal(scale, nsxData(nsxIndex).scaleLims{j}), 'Scale limits for channel %s.%s (id == %d) does not match first channel in group', groupName, names{j}, ids(j));
                    end
                end
                
                lookup.groupName = groupName;
                lookup.names = names;
                lookup.nsxIndex = nsxIndex;
                lookup.idxWithinNsx = idxWithinNsx;
                lookup.channelIds = nsxData(nsxIndex).channelIds(idxWithinNsx);
                lookup.single = false;
                lookup.lookup = find(found);
                lookup.scaleFn = nsxData(nsxIndex).scaleFns{1};
                lookup.scaleLims = nsxData(nsxIndex).scaleLims{1};
                lookup.units = getCheckUnits(groupName, nsxData, nsxIndex, idxWithinNsx);
            end
        end
        
        % found it, add to the lookup list
        if isempty(analogLookup)
            analogLookup = lookup;
        else
            analogLookup(end+1) = lookup; %#ok<AGROW>
        end
    end

    if keepExtraAnalogChannels
        % nsxDataChannelIsExtraMask{insx} is true for all signals that
        % weren't captured by name above. First filter those masks so that
        % each channel id is only true in the nsx with the highest sampling
        % rate
        for insx = 1:numel(nsxData)
            indsExtra = find(nsxDataChannelIsExtraMask{insx});
            for iInd = 1:numel(indsExtra)
                id = nsxData(insx).channelIds(indsExtra(iInd));
                [~, bestInsx] = findChannelIdInNsxData(nsxData, id, true(numel(nsxData), 1), nsxDataChannelIsExtraMask);
                if insx ~= bestInsx
                    nsxDataChannelIsExtraMask{insx}(indsExtra(iInd)) = false;     
                end
            end
        end
        
        % build a lookup for all the extra channels that we didn't use
        % group them by ns index since they'll have the same sampling rate
        lastNeuralChannelId = 128;
        for insx = 1:numel(nsxData)
            % need to split channels between [1 128] and [129 - Inf]
            % because front panel channels have different scaling limits
            for channelGroup = 1:2
                if channelGroup == 1
                    idxWithinNsx = find(nsxDataChannelIsExtraMask{insx} & nsxData(insx).channelIds <= lastNeuralChannelId);
                    groupSuffix = '_neural';
                else
                    idxWithinNsx = find(nsxDataChannelIsExtraMask{insx} & nsxData(insx).channelIds > lastNeuralChannelId);
                    groupSuffix = '_frontPanel';
                end
                    
                if ~any(idxWithinNsx), continue, end
            
                if nsxData(insx).ext(1) == '.'
                    lookup.groupName =  sprintf('%s%s', nsxData(insx).ext(2:end), groupSuffix);
                else
                    lookup.groupName =  sprintf('%s%s', nsxData(insx).ext, groupSuffix);
                end
                nsxIndex = insx;

                lookup.names = nsxData(nsxIndex).chLabels(idxWithinNsx);
                lookup.nsxIndex = nsxIndex;
                lookup.idxWithinNsx = idxWithinNsx;
                lookup.channelIds = nsxData(nsxIndex).channelIds(idxWithinNsx);
                lookup.single = false;
                lookup.lookup = 1:numel(idxWithinNsx);

                scale = nsxData(nsxIndex).scaleLims{1};
                for j = 2:numel(idxWithinNsx)
                    assert(isequal(scale, nsxData(nsxIndex).scaleLims{idxWithinNsx(j)}), 'Scale limits for extra channel %s does not match first channel in group', lookup.names{j});
                end

                lookup.scaleFn = nsxData(nsxIndex).scaleFns{1};
                lookup.scaleLims = nsxData(nsxIndex).scaleLims{1};
                lookup.units = getCheckUnits(lookup.groupName, nsxData, nsxIndex, idxWithinNsx);

                if isempty(analogLookup)
                    analogLookup = lookup;
                else
                    analogLookup(end+1) = lookup; %#ok<AGROW>
                end
            end
        end
    end
    
end

function [found, bestInsx, idxWithinNsx, nsxDataChannelIsExtraMask] = findChannelIdInNsxData(nsxData, id, nsxMask, nsxDataChannelIsExtraMask) % loop thru the nsx files looking for the highest sampling rate
    % nsx file with this channel
    bestInsx = NaN;
    idxWithinNsx = NaN;
    found = false;
    for insx = 1:length(nsxData)
        if ~nsxMask(insx), continue; end
        channelIdList = nsxData(insx).channelIds;
        ind = find(channelIdList == id); 
        if ~isempty(ind)
            % mark this one as used since we don't want to mark
            % it as extra
            nsxDataChannelIsExtraMask{insx}(ind) = false;

            if ~isnan(bestInsx) && nsxData(insx).samplingHz < nsxData(bestInsx).samplingHz
                % ensure we take the nsx file with the highest
                % sampling rate - this one has been sampled at a
                % higher rate already
                continue;
            end

            bestInsx = insx;
            found = true;
            idxWithinNsx = ind;
        end
    end
end

function units = getCheckUnits(groupName, nsxData, nsxIndex, idxWithinNsx)
    trim = @(s) strtok(s, char(0));

    if numel(nsxData(nsxIndex).chUnits) == 1
        units = trim(nsxData(nsxIndex).chUnits{1});
    else
        units = arrayfun(@(ind) trim(nsxData(nsxIndex).chUnits{ind}), idxWithinNsx, 'UniformOutput', false);
        if numel(unique(units)) > 1
            warning('Channels in group %s use different units: %s', groupName, strjoin(unique(units)));
        end
        units = units{1};  
    end
end
