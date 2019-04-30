function td = addSegmentedDatasetToTrialData(td, seg, varargin)
    p = inputParser();
    p.addParameter('arrayPrefix', 'npix', @ischar);
    p.addParameter('cluster_groups', ["good", "mua", "noise", "unsorted"], @(x) isstring(x) || iscategorical(x))
    p.addParameter('separateArrayPerClusterGroup', false, @islogical);
    p.addParameter('sortByClusterPosition', false, @islogical);
    p.addParameter('cluster_ids', seg.cluster_ids, @isvector);
    p.parse(varargin{:});

    td = td.reset();
    assert(td.nTrials == seg.nTrials, 'Trial counts do not match');

    if p.Results.sortByClusterPosition
        debug('Computing cluster centers of mass\n');
        m = seg.dataset.getMetrics();
        clusterCOM = m.cluster_centerOfMass;
    end
    
    cluster_ids = p.Results.cluster_ids;
    clusterInds = seg.lookup_clusterIds(cluster_ids);
    
    cluster_groups_keep = categorical(p.Results.cluster_groups);
    
    if p.Results.separateArrayPerClusterGroup
        for iG = 1:numel(cluster_groups_keep)
            
            cluster_mask = ismember(seg.cluster_groups(clusterInds), cluster_groups_keep(iG));
            if ~any(cluster_mask)
                continue;
            end
            
            clusterIndsThis = clusterInds(cluster_mask);
            if p.Results.sortByClusterPosition
                clusterIndsThis = sortClusters(clusterCOM, clusterIndsThis);
            end
            electrodes = seg.cluster_ids(clusterIndsThis);
            units = zeros(numel(electrodes), 1);
        
            data = seg.spike_times_ms_rel_start(:, clusterIndsThis);
                
            arrayName = sprintf('%s_%s', p.Results.arrayPrefix, cluster_groups_keep(iG));
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);
            
            % add metadata
            td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
            td = td.setChannelMetaKey(arrayName, 'cluster_idx', seg.cluster_ids(clusterIndsThis));
            td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(clusterIndsThis));
        end
    else
        cluster_mask = ismember(seg.cluster_groups(clusterInds), cluster_groups_keep);
        if any(cluster_mask)
            clusterInds = clusterInds(cluster_mask);
            if p.Results.sortByClusterPosition
                clusterInds = sortClusters(clusterCOM, clusterInds);
            end
            electrodes = seg.cluster_ids(clusterInds);
            units = zeros(numel(electrodes), 1);

            data = seg.spike_times_ms_rel_start(:, clusterInds);

            arrayName = p.Results.arrayPrefix;
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);

            % add metadata
            td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
            td = td.setChannelMetaKey(arrayName, 'cluster_idx', seg.cluster_ids(clusterInds));
            td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(clusterInds));
        end
    end

    td = td.addOrUpdateBooleanParam(sprintf('%s_hasData', p.Results.arrayPrefix), seg.trial_has_data);  
end

function clusterInds = sortClusters(clusterCOM, clusterInds)
    y = clusterCOM(clusterInds, 2);
    [~, sortIdx] = sort(y, 1, 'descend');
    clusterInds = clusterInds(sortIdx);
end