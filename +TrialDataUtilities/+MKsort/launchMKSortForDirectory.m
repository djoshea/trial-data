function launchMKSortForDirectory(d)

    h = mksort;
    h.HandleVisibility = 'on';
    handles = guidata(h);
    set(h, 'CurrentAxes', handles.axMulti);
    handles.dataDir = d;
    mksort('loadDirectory', handles);

end