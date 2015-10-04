function td = remapNeuralChannels(td, spikeMapFn, continuousMapFn)
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
end

end