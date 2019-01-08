function td = addSegmentedDatasetToTrialData(td, seg, varargin)
    p = inputParser();
    p.addParameter('arrayPrefix', 'npix', @ischar);
    p.addParameter('cluster_groups', ["good", "mua"], @(x) isstring(x) || iscategorical(x))
    p.addParameter('separateArrayPerClusterGroup', false, @islogical);
    p.addParameter('sortByClusterPosition', true, @islogical);
    p.parse(varargin{:});

    assert(td.nTrials == seg.nTrials, 'Trial counts do not match');

    clusterCOM = seg.dataset.computeCenterOfMassLocationByCluster();
    
    cluster_groups_keep = categorical(p.Results.cluster_groups);
    
    if p.Results.separateArrayPerClusterGroup
        for iG = 1:numel(cluster_groups_keep)
            cluster_mask = ismember(seg.cluster_groups, cluster_groups_keep(iG));
   
            if ~any(cluster_mask)
                continue;
            end
            
            if p.Results.sortByClusterPosition
                cluster_ind = sortClusters(cluster_mask);
            else
                cluster_ind = find(cluster_mask);
            end
            electrodes = seg.cluster_ids(cluster_ind);
            units = zeros(numel(electrodes), 1);
        
            data = seg.spike_times_ms_rel_start(:, cluster_ind);
                
            arrayName = sprintf('%s_%s', p.Results.arrayPrefix, cluster_groups_keep(iG));
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);
            
            % add metadata
            td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
            td = td.setChannelMetaKey(arrayName, 'cluster_idx', seg.cluster_ids(cluster_ind));
            td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(cluster_ind));
        end
    else
        cluster_mask = ismember(seg.cluster_groups, cluster_groups_keep);
   
        if any(cluster_mask)  
            if p.Results.sortByClusterPosition
                cluster_ind = sortClusters(cluster_mask);
            else
                cluster_ind = find(cluster_mask);
            end
            electrodes = seg.cluster_ids(cluster_ind);
            units = zeros(numel(electrodes), 1);

            data = seg.spike_times_ms_rel_start(:, cluster_ind);

            arrayName = p.Results.arrayPrefix;
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);
            
            % add metadata
            td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
            td = td.setChannelMetaKey(arrayName, 'cluster_idx', seg.cluster_ids(cluster_ind));
            td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(cluster_ind));
        end
    end

    td = td.addOrUpdateBooleanParam(sprintf('%s_hasData', p.Results.arrayPrefix), seg.trial_has_data);
    

function takeIdx = sortClusters(cluster_mask)
    takeIdx = find(cluster_mask);
    y = clusterCOM(takeIdx, 2);
    [~, sortIdx] = sort(y, 1, 'descend');
    takeIdx = takeIdx(sortIdx);
end

end