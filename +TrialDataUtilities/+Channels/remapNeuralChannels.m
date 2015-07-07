function td = remapNeuralChannels(td, spikeMapFn, lfpMapFn)
% spikeMapFn and lfpMapFn receive channel descriptors and return the new
% name for the channel, and the meta field.
% [newName, meta] = spikeMapFn(spikeChannelDescriptor)
% [newName, meta] = lfpMapFn(lfpChannelDescriptor)

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

if ~isempty(lfpMapFn)
    lfpCh = td.listLFPChannels();
    prog = ProgressBar(numel(lfpCh), 'Mapping spike channels');
    for i = 1:numel(lfpCh)
        ch = lfpCh{i};
        [newCh, meta] = lfpMapFn(td.channelDescriptorsByName.(ch)); 
        prog.update(i, 'Mapping spike channel %s to %s', ch, newCh);
        td = td.renameChannel(ch, newCh);
        td = td.setChannelMeta(newCh, meta);
    end
    prog.finish();
end

end