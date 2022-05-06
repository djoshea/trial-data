function td = addSegmentedDatasetToTrialData(td, seg, varargin)
    p = inputParser();
    p.addParameter('arrayPrefix', 'npix', @ischar);
    p.addParameter('cluster_groups', ["good", "mua", "noise", "unsorted"], @(x) isstring(x) || iscategorical(x))
    %p.addParameter('separateArrayPerClusterGroup', false, @islogical);
    p.addParameter('sortByClusterPosition', false, @islogical);
    p.addParameter('cluster_ids', seg.cluster_ids, @isvector);
    p.addParameter('cluster_ratings', [], @(x) true);
    
    p.addParameter('markNoDataIfNoSpikes', true, @islogical); % if a trial has 0 spikes in the segmented KS dataset, consider it as having no data in the file
    p.addParameter('markNoDataFilterContiguousTrials', 10, @isscalar); % filter trials with spikes to be in contigous groups of at least this many to count (removes patchy ends from the group)
    
    p.addParameter('addHasDataParam', true, @islogical);
    p.addParameter('include_valid', true, @islogical);
    p.addParameter('include_cutoff', false, @islogical);

    p.addParameter('convertToMs', true, @islogical); % if true, spike times converted to ms, if false, left as samples to be scaled dynamically by TrialData

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
    
    cluster_ratings = p.Results.cluster_ratings;
    if ~isempty(cluster_ratings)
        assert(numel(cluster_ratings) == numel(cluster_ids));
    end

    cluster_groups_keep = categorical(p.Results.cluster_groups);
    
    include_valid = p.Results.include_valid;
    include_cutoff = p.Results.include_cutoff;
    convertToMs = p.Results.convertToMs;
    
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
                if convertToMs
                    data_valid = seg.spike_times_ms_rel_start(:, clusterInds);
                    data_cutoff = seg.cutoff_spike_times_ms_rel_start(:, clusterInds);
                else
                    data_valid = seg.spike_times_rel_start(:, clusterInds);
                    data_cutoff = seg.cutoff_spike_times_rel_start(:, clusterInds);
                end
                data = cellfun(@(x1, x2) sort(cat(1, x1, x2)), data_valid, data_cutoff);
            else
                if convertToMs
                    data = seg.spike_times_ms_rel_start(:, clusterInds);
                else
                    data = seg.spike_times_rel_start(:, clusterInds);
                end
            end
        elseif include_cutoff
            if convertToMs
                data = seg.cutoff_spike_times_ms_rel_start(:, clusterInds);
            else
                data = seg.spike_times_ms_rel_start(:, clusterInds);
            end
        else
            error('Either valid (include_valid==true) or cutoff (include_cutoff) spikes must be included');
        end

        % data is now referenced to trial_start, not to "TimeZero", so we pass it in as isAligned==true, 
        % so that this offset is handled for us automatically
        arrayName = p.Results.arrayPrefix;

        if convertToMs
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', true, 'raw', true);
        else
            scalingToMs = 1000 / seg.fsAP;
            data = cellfun(@uint32, data, 'UniformOutput', false);
            td = td.addSpikeArrayChannel(arrayName, electrodes, units, data, 'isAligned', true, 'raw', true, ...
                'timeScaling', scalingToMs, 'timeOriginalDataClass', 'uint32');
        end

        % add metadata
        td = td.setChannelMetaKey(arrayName, 'kilosort_pathLeaf', seg.dataset.pathLeaf);
        td = td.setChannelMetaKey(arrayName, 'cluster_ids', seg.cluster_ids(clusterInds));
        td = td.setChannelMetaKey(arrayName, 'cluster_groups', seg.cluster_groups(clusterInds));
        td = td.setChannelMetaKey(arrayName, 'includes_valid_spikes', include_valid);
        td = td.setChannelMetaKey(arrayName, 'includes_cutoff_spikes', include_cutoff);
        td = td.setChannelMetaKey(arrayName, 'cluster_ratings', cluster_ratings);
    end

    if p.Results.addHasDataParam
        ch_has_data = sprintf('%s_hasData', p.Results.arrayPrefix);
        
        if p.Results.markNoDataIfNoSpikes
            trial_has_data = seg.trial_has_data & seg.trial_has_nonzero_spikes;
        else
            trial_has_data = seg.trial_has_data;
        end
        
        if p.Results.markNoDataFilterContiguousTrials > 0
            % try to find a band where trial_has_data is true and trial_has_nonzero_spikes > 0.1 median
            mask_contiguous = imopen(seg.trial_has_nonzero_spikes, ones(p.Results.markNoDataFilterContiguousTrials, 1)); 
            trial_has_data = trial_has_data & mask_contiguous;
        end
        
        td = td.addOrUpdateBooleanParam(ch_has_data, trial_has_data);
        td = td.setChannelDisplayGroup(ch_has_data, 'neuropixel');
    end
end

function clusterInds = sortClusters(clusterCOM, clusterInds)
    y = clusterCOM(clusterInds, 2);
    [~, sortIdx] = sort(y, 1, 'descend');
    clusterInds = clusterInds(sortIdx);
end