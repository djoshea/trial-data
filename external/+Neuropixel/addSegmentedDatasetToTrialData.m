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
                cluster_idx = sortClusters(cluster_mask);
            else
                cluster_idx = find(cluster_mask);
            end
            electrodes = seg.cluster_ids(cluster_idx);
            units = zeros(numel(electrodes), 1);
        
            data = seg.spike_times_ms_rel_start(:, cluster_idx);
                
            arrayName = sprintf('%s_%s', p.Results.arrayPrefix, cluster_groups_keep(iG));
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);
        end
    else
        cluster_mask = ismember(seg.cluster_groups, cluster_groups_keep);
   
        if any(cluster_mask)  
            if p.Results.sortByClusterPosition
                cluster_idx = sortClusters(cluster_mask);
            else
                cluster_idx = find(cluster_mask);
            end
            electrodes = seg.cluster_ids(cluster_idx);
            units = zeros(numel(electrodes), 1);

            data = seg.spike_times_ms_rel_start(:, cluster_idx);

            arrayName = p.Results.arrayPrefix;
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', false, 'raw', true);
        end
    end

    td = td.addOrUpdateBooleanParam(sprintf('%s_hasData', p.Results.arrayPrefix), seg.trial_has_data);

function takeIdx = sortClusters(cluster_mask)
    y = clusterCOM(cluster_mask, 2);
    [~, takeIdx] = sort(y, 1, 'descend');
end

end