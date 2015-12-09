function launchMKSortSingleChannel(exportData)

    iC = 1;
    channelInfo = [];
    channelInfo.sorts = exportData.sorts(iC);
    channelInfo.preview = exportData.previews(iC);
    
    % save wave file to /tmp
    datadir = tempdir;
    waveFile = fullfile(datadir, exportData.wavefileNameShort{iC});
    
    waveforms = exportData.waveforms(iC);
    save(waveFile,  'waveforms', '-v7');
    channelInfo.path = datadir;
    channelInfo.array = exportData.previews(iC).array;
    channelInfo.thisChannel = iC;
    hSort = oneChannelSorter(channelInfo);
    
    hSort.HandleVisibility = 'on';
    handles = guidata(hSort);
%     handles.dataDir = channelInfo.datadir;
    
end