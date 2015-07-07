function td = addOffsetToAllSpikeTimes(td, tOffset)
% manually add a time offset to all spike times inside a trial data

spikeCh = td.listSpikeChannels();
prog = ProgressBar(numel(spikeCh), 'Adding time offset to spike channels');
for iC = 1:numel(spikeCh)
   chName = spikeCh{iC};
   prog.update(iC, 'Adding time offset to %s spikes', chName);
   times = td.getSpikeTimesUnaligned(chName);
   times = cellfun(@(t) t + tOffset, times, 'UniformOutput', false);
   td = td.setSpikeChannel(chName, times, 'isAligned', false);
end
prog.finish();

end