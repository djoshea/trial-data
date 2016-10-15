function td = mergeSortsIntoTrialData(td, handSortPath, varargin)

    p = inputParser();
    p.addParameter('keepUnsortedUnit', false, @islogical);
    p.addParamValue('useAlignedWaves', true, @islogical); % may make repeat sorting impossible since wave data may be dropped
    p.addParamValue('sortMethod', 'mksort', @ischar); % what to specify as the sort method, for reference only
    p.parse(varargin{:});

    % Load sorts structure, figure out which channels are worth including
    loadVar = load(fullfile(handSortPath, 'sorts.mat'));
    sorts = loadVar.sorts;

    loadVar = load(fullfile(handSortPath, 'exportMeta.mat'));
    exportMeta = loadVar.exportMeta;
    offsetByTrial = exportMeta.offsetByTrial;

    % check trial counts match sort data
    assert(td.nTrials == exportMeta.nTrials, ...
        'TrialData has %d trials but sort data has %d trials', td.nTrials, exportMeta.nTrials);
    nTrials = td.nTrials;

    td = td.reset();
    
    td = td.dropChannels(td.listSpikeChannels());

    % build list of array/electrodes saved into waveform files
    waveformFiles = {};
    alpha = alphabet();
    for ch = 1:length(sorts)
      if sorts(ch).onlineSorted || sorts(ch).userSorted
        if isnan(sorts(ch).array)
          lett = '';
        else
          lett = alpha(sorts(ch).array);
        end
        waveformFiles{end+1} = makeWaveFilename(lett, sorts(ch).electrode); %#ok<AGROW>
      end
    end
    nFiles = length(waveformFiles);

    % loop over array/electrodes
    prog = ProgressBar(nFiles, 'Merging spikes into TrialData');
    for wf = 1:length(waveformFiles)
        prog.update(wf, 'Merging spikes from %s', waveformFiles{wf});

        loadVar = load(fullfile(handSortPath, waveformFiles{wf}));
        w = loadVar.waveforms;
        electrode = w.electrode;
        arrayIdx = w.array;
        channelMask = [sorts.electrode] == electrode & [sorts.array] == arrayIdx;
        sortInfo = sorts(channelMask);

        arrayName = exportMeta.arrayNames{arrayIdx};

        % for generating new channel names for these units
        spikeChNameFn = @(unit) sprintf('%s%02d_%d', arrayName, electrode, unit);

        units = unique(w.units);
        if ~p.Results.keepUnsortedUnit
            units = setdiff(units, 0);
        end

        if p.Results.useAlignedWaves && ~isempty(w.alignedWaves)
            waves = w.alignedWaves;
        else
            waves = w.waves;
        end
        nWaveSamples = size(waves, 1);
        waveTvec = w.waveformsTime;
        if numel(waveTvec) > nWaveSamples
            waveTvec = waveTvec(1:nWaveSamples);
        elseif numel(waveTvec) < nWaveSamples
            error('Waveforms time vector too short');
        end
        
        for iU = 1:numel(units)
            spikeChName = spikeChNameFn(units(iU));
            [spikesByTrial, wavesByTrial] = cellvec(nTrials);

            progT = ProgressBar(nTrials, 'Segmenting spikes for unit %s', spikeChName);
            for iT = 1:nTrials
                if mod(iT, 10) == 0
                    progT.update(iT);
                end
                mask = w.units == units(iU) & w.trialIdxBySpike == iT;
                % transpose back (nSpikes along the first dimension)
                spikesByTrial{iT} = w.spikeTimes(mask)' - offsetByTrial(iT);
                wavesByTrial{iT} = waves(:, mask)';
            end
            progT.finish();

            ratingsByTrial = convertRatingsByEpochToRatingsByTrial(units(iU), w.ratings, offsetByTrial);
            sortQuality = round(nanmean(ratingsByTrial), 1);
            
            td = td.addSpikeChannel(spikeChName, spikesByTrial, 'waveforms', wavesByTrial,...
                'waveformsTime', waveTvec, 'sortQuality', sortQuality, ...
                'sortMethod', p.Results.sortMethod, 'sortQualityEachTrial', ratingsByTrial, ...
                'updateValidOnly', false); 
        end
    end
    prog.finish();
    
    td = td.setMetaKey('spikeSorting', 'hand-sorted via mksort');
    
end

function ratingsByTrial = convertRatingsByEpochToRatingsByTrial(unit, ratings, offsetsByTrial)
    ratingsByTrial = nanvec(numel(offsetsByTrial));
    for iR = 1:numel(ratings)
        if isempty(ratings(iR).epoch)
            ratingsByTrial(:) = ratings(iR).ratings(unit);
        else
            tStart = ratings(iR).epoch(1);
            tStop = ratings(iR).epoch(2);
            trialMask = offsetsByTrial >= tStart & offsetsByTrial <= tStop;
            ratingsByTrial(trialMask) = ratings(iR).ratings(unit);
        end  
    end
end


