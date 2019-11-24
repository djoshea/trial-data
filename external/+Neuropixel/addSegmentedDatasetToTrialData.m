function td = addSegmentedDatasetToTrialData(td, seg, varargin)
    p = inputParser();
    p.addParameter('arrayPrefix', 'npix', @ischar);
    p.addParameter('cluster_groups', ["good", "mua", "noise", "unsorted"], @(x) isstring(x) || iscategorical(x))
    %p.addParameter('separateArrayPerClusterGroup', false, @islogical);
    p.addParameter('sortByClusterPosition', false, @islogical);
    p.addParameter('cluster_ids', seg.cluster_ids, @isvector);
    
    p.addParameter('addHasDataParam', true, @islogical);
    p.addParameter('include_valid', true, @islogical);
    p.addParameter('include_cutoff', false, @islogical);
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
    
    include_valid = p.Results.include_valid;
    include_cutoff = p.Results.include_cutoff;
    
    cluster_mask = ismember(seg.cluster_groups(clusterInds), cluster_groups_keep);
    if any(cluster_mask)
        clusterInds = clusterInds(cluster_mask);
        if p.Results.sortByClusterPosition
            clusterInds = sortClusters(clusterCOM, clusterInds);
        end
        
        % name as array#unit# since electrode isn't meaningful
        units = seg.cluster_ids(clusterInds);
        electrodes = nan(size(units));
        if include_valid
            if include_cutoff
                data_valid = seg.spike_times_ms_rel_start(:, clusterInds);
                data_cutoff = seg.cutoff_spike_times_ms_rel_start(:, clusterInds);
                data = cellfun(@(x1, x2) sort(cat(1, x1, x2)), data_valid, data_cutoff);
            else
                data = seg.spike_times_ms_rel_start(:, clusterInds);
            end
        elseif include_cutoff
            data = seg.cutoff_spike_times_ms_rel_start(:, clusterInds);
        else
            error('Either valid (include_valid==true) or cutoff (include_cutoff) spikes must be included');
        end
        
        % data is now referenced to trial_start, not to "TimeZero", so we pass it in as isAligned==true, 
        % so that this offset is handled for us automatically
        arrayName = p.Results.arrayPrefix;
        td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', true, 'raw', true);

        % add metadata
        td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
        td = td.setChannelMetaKey(arrayName, 'cluster_ids', seg.cluster_ids(clusterInds));
        td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(clusterInds));
        td = td.setChannelMetaKey(arrayName, 'includes_valid_spikes', include_valid);
        td = td.setChannelMetaKey(arrayName, 'includes_cutoff_spikes', include_cutoff);
    end

    if p.Results.addHasDataParam
        ch_has_data = sprintf('%s_hasData', p.Results.arrayPrefix);
        td = td.addOrUpdateBooleanParam(ch_has_data, seg.trial_has_data);
        td = td.setChannelDisplayGroup(ch_has_data, 'neuropixel');
    end
    
    if isfield(td.saveFastPartitionInfo, 'spikes')
        td.saveFastPartitionInfo.spikes = union(td.saveFastPartitionInfo.spikes, string(arrayName));
    else
        td.saveFastPartitionInfo.spikes = string(arrayName);
    end
end

function clusterInds = sortClusters(clusterCOM, clusterInds)
    y = clusterCOM(clusterInds, 2);
    [~, sortIdx] = sort(y, 1, 'descend');
    clusterInds = clusterInds(sortIdx);
end