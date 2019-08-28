function td = remapNeuralChannels(td, spikeMapFn, continuousMapFn)
%% THIS FUNCTION HAS NOW BEEN MADE OBSOLETE BY SPIKE ARRAYS

error('Not useful anymore');


% spikeMapFn and continuousMapFn receive channel descriptors and return the new
% name for the channel, and the meta field.
% [newName, meta] = spikeMapFn(spikeChannelDescriptor)
% [newName, meta] = continuousMapFn(continuousChannelDescriptor)

td.warnIfNoArgOut(nargout);

if ~isempty(spikeMapFn)
    spikeCh = td.listSpikeChannels();
    prog = ProgressBar(numel(spikeCh), 'Mapping spike channels');
    for i = 1:numel(spikeCh)
        ch = spikeCh{i};
        [newCh, meta] = spikeMapFn(td.channelDescriptorsByName.(ch)); 
        prog.update(i, 'Mapping spike channel %s to %s', ch, newCh);
        td = td.renameChannel(ch, newCh);
        td = td.setChannelMeta(newCh, meta);
    end
    prog.finish();
end

if ~isempty(continuousMapFn)
    continuousCh = td.listContinuousNeuralChannels();
    prog = ProgressBar(numel(continuousCh), 'Mapping spike channels');
    for i = 1:numel(continuousCh)
        ch = continuousCh{i};
        [newCh, meta] = continuousMapFn(td.channelDescriptorsByName.(ch)); 
        prog.update(i, 'Mapping continuous neural channel %s to %s', ch, newCh);
        td = td.renameChannel(ch, newCh);
        td = td.setChannelMeta(newCh, meta);
    end
    prog.finish();
    
    [groups, channelsInGroup] = td.listContinuousNeuralChannelGroups();
    
    % rename continuous neural channel groups to match their contents if
    % their contents are all of one type and on the same array
    for g = 1:numel(groups)
        cdsThisGroup = td.getChannelDescriptorMulti(channelsInGroup{g});
        cdsThisGroup = cat(1, cdsThisGroup{:});
        arrays = unique({cdsThisGroup.array});
        types = unique({cdsThisGroup.type});
        
        if numel(arrays) == 1 && numel(types) == 1
            % this group should be renamed to match its contents
            newName = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArray(types{1}, arrays{1});
            td = td.renameChannel(groups{g}, newName);
        end
    end
end

end