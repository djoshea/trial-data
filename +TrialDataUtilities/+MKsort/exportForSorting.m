function exportForSorting(td, outDir, varargin)
% artifactThresh (optional) specifies how big a point in a waveform is
% allowed to be before we consider it an artifact and throw it away. Sign
% does not matter, specify in mV. Default is Inf.
%
% snippetLen is an optional argument specifying the maximum length of
% snippets to extract. Default 32.
%
% firstSnippetPt is an optional argument specifying the first point of the
% snippet within the .nev file to save (i.e., it will save points
% firstSnippetPt : firstSnippetPt+snippetLen-1). Default 1.

p = inputParser;
% p.addParamValue('artifactThresh', Inf, @isscalar);
% p.addParamValue('conserveMemory', true, @islogical);
p.addParamValue('maxUnits', 4, @isscalar);
p.addParamValue('lockoutMs', 1.6, @isscalar);
p.parse(varargin{:});

if ~isa(td, 'TrialDataConditionAlign')
    warning('Provide TrialDataConditionAlign instance for condition info');
end

maxUnits = p.Results.maxUnits;
lockout = p.Results.lockoutMs;

if ~exist(outDir, 'dir')
    mkdirRecursive(outDir);
end

% Gather spike channels from this TrialData
spikeChannels = td.listSpikeChannels();
nS = numel(spikeChannels);
arrays = cellfun(@(ch) td.channelDescriptorsByName.(ch).array, spikeChannels, 'UniformOutput', false);
electrodes = cellfun(@(ch) td.channelDescriptorsByName.(ch).electrode, spikeChannels);
units = cellfun(@(ch) td.channelDescriptorsByName.(ch).unit, spikeChannels);

% convert array names into array indices (and later into indexes into
% alphabet because mksort wants it this way)
[uniqueArrayNames, ~, arrayIdx] = unique(arrays);
alpha = alphabet();

% find unique array electrode combinations
[uniqAE, ~, whichAE] = unique([arrayIdx, electrodes], 'rows');
% nAE = numel(uniqAE);

% save array names in the meta file so we can use them when merging
exportMeta.arrayNames = uniqueArrayNames;
exportMeta.nTrials = td.nTrials;

% reset condition info and validity so that all spikes are sorted
% regardless of the current configuration
tdReset = td.reset();
nTrials = tdReset.nTrials;

% % REMOVE THIS
% uniqAE = uniqAE(1:5, :);

chCounter = 0;
chWithSpikesCounter = 0;
waveTimeVec = [];

% here we'll stack the timestamps of each trial so that they can
% concatenate with non-overlapping timestamps
durationsByTrial = tdReset.getDurationsRaw();
cumDurations = cumsum(durationsByTrial);
offsetByTrial = [0; cumDurations(1:end-1)];

exportMeta.offsetByTrial = offsetByTrial; %#ok<STRNU> % later saved to disk

% figure out trial info for waveforms struct, this is used to display
% consistency of tuning over time. we'll condition using the present
% conditions found in the trial data, and use the align start and stop in
% order to focus on a specific time window within each trial to count
% spikes within (rather than include useless times within each trial)
if isa(td, 'TrialDataConditionAlign')
    addTrialInfo = true;
    [tStartRel, tStopRel] = td.getTimeStartStopZeroRelativeToTrialStartEachTrial();
    trialInfo.condition = td.conditionIdx';
    trialInfo.trialStartTimes = (offsetByTrial + tStartRel)';
    trialInfo.trialEndTimes = (offsetByTrial + tStopRel)';
else
    addTrialInfo = false;
end

prog = ProgressBar(nS, 'Extracting spikes for each unit');
for iAE = 1:size(uniqAE, 1)
    % find spike channels on this array, electrode
    matchingIdx = find(whichAE == iAE);
    matchingCh = spikeChannels(matchingIdx);
    nMatch = numel(matchingCh);
    
    % gather spike waveforms, times, and trial ids across units on this
    % electrode, we'll combine them later
    [wavesCell, timesCell, trialIdsCell, trialIdxCell] = deal(cellvec(nMatch));
    for iM = 1:nMatch
        chIdx = matchingIdx(iM);
        chName = spikeChannels{chIdx};
        thisArray = arrayIdx(chIdx);
        thisArrayStr = arrays{chIdx};
        thisElectrode = electrodes(chIdx);
        thisUnit = units(chIdx);
        chCounter = chCounter + 1;
        if isempty(thisArray)
            prog.update(chCounter, 'Extracting waveforms for electrode %d unit %d', thisElectrode, thisUnit);
        else
            prog.update(chCounter, 'Extracting waveforms for array ''%s'' electrode %d unit %d', thisArrayStr, thisElectrode, thisUnit);
        end

        [wavesByTrial, thisWaveTimeVec, timesByTrial] = tdReset.getSpikeWaveforms(chName);
        
        % assemble lenWave x nSpikes matrix of waveforms
        wavesCell{iM} = cat(1, wavesByTrial{:});
        % assemble nSpikes x 1 vector of times in ms
        % add the right offset to each trial
        timesByTrialOffset = cellfun(@(times, off) times + off, timesByTrial, num2cell(offsetByTrial), 'UniformOutput', false);
        [timesCell{iM}, trialIdxCell{iM}] = TensorUtils.catWhich(1, timesByTrialOffset{:});
        
        % compare waveTimeVecs to make sure they stay constant
        if isempty(waveTimeVec)
            waveTimeVec = thisWaveTimeVec;
        else
            assert(isequal(waveTimeVec, thisWaveTimeVec), 'Wave time vector for %s is not consistent with other spike channels', chName);
        end
        
        if addTrialInfo
            % build a cell which identifies the trial for each spike, or
            % NaN if the spike falls outside the valid window
            trialIdBySpikeCell = cellvec(nTrials);
            for iT = 1:nTrials
                if td.valid(iT)
                    idx = repmat(iT, numel(timesByTrial{iT}), 1);
                    idx(timesByTrial{iT} < tStartRel(iT) | timesByTrial{iT} > tStopRel(iT)) = NaN;
                else
                    idx = nanvec(numel(timesByTrial{iT}));
                end
                trialIdBySpikeCell{iT} = idx;
            end
            trialIdsCell{iM} = cat(1, trialIdBySpikeCell{:});
        end
    end
    
    % concatenate the waves and units on this electrode
    wavesCat = cat(1, wavesCell{:});
    [timesCat, unitIdxCat] = TensorUtils.catWhich(1, timesCell{:});
    matchingUnits = units(matchingIdx);
    unitsCat = matchingUnits(unitIdxCat);
    trialIdxCat = cat(1, trialIdxCell{:});
    
    waveforms = createWaveformsStruct(thisElectrode, wavesCat', unitsCat', timesCat', thisArray, maxUnits); 

    % add trial number for each unit to simplify the merging later on
    waveforms.trialIdxBySpike = trialIdxCat';
    
    % keep track of the waveforms time vector
    waveforms.waveformsTime = waveTimeVec;
    
    % add trial info for tuning analysis
    if addTrialInfo
        waveforms.trialInfo = trialInfo;
        waveforms.trialInfo.trial = cat(1, trialIdsCell{:})'; % mksort expects row vec
    end
    
    if ~isempty(waveforms.waves)
      %% Save channel
      waveFilename = makeWaveFilename(alpha(thisArray), thisElectrode);
      save(fullfile(outDir, waveFilename), 'waveforms', '-v7');
      
      %% Generate this channel's worth of previews and sorts
      chWithSpikesCounter = chWithSpikesCounter + 1;
      [previews(chWithSpikesCounter), sorts(chWithSpikesCounter)] = generateSortsAndPreviews(waveforms, maxUnits, lockout); %#ok<AGROW,NASGU>
    end
end
prog.finish();

%% Save previews and sorts
debug('Saving sorts, previews, and exportMeta files\n');
save(fullfile(outDir, 'previews.mat'), 'previews', '-v6');
save(fullfile(outDir, 'sorts.mat'), 'sorts', '-v6');

save(fullfile(outDir, 'exportMeta.mat'), 'exportMeta', '-v6');

return;

% for arr = 1:length(nevNames)
%   for fileNum = 1:length(nevNames{arr})
%  
%     %% Figure out lockout period
%     % This is the number of waveform samples divided by the sampling rate,
%     % converted into ms. This may be overridden at the end if lockout
%     % violations are found. Violations imply either that RRR was used, which
%     % has a shorter lockout, or that a buggy version of the Cerebus software
%     % was used, which can add extraneous events.
%     snippetLockout = 1000 * ((header.datasize - 8) / waveBytesPerSample) / header.SampleRes;
%     
%     
%     %% Check for lockout violations
%     [lockout, lockoutViolations, flagrantLockoutViolations] = ...
%       checkForLockoutViolations(spikeTimes, allElectrodeNums, snippetLockout, lockoutViolations, flagrantLockoutViolations, RRRlockout);


% 
% %% Stitch together chunks into whole channels and generate previews
% fprintf('Stitching chunks together, generating previews\n');
% fprintf('Channel: ');
% 
% % Loop through arrays
% usedChs = 0;
% for array = 1:size(allChunkFileInfo, 1)
%   % Array letter if applicable
%   if size(allChunkFileInfo, 1) > 1
%     lett = alpha(array);
%     fprintf('%s', lett);
%   else
%     lett = '';
%   end
%   
%   % Loop through electrodes
%   for electrode = 1:size(allChunkFileInfo, 2)
%     %% Figure out whether this was a used channel
%     % Figure out how many chunks for this array/electrode
%     chunkInfo = allChunkFileInfo{array, electrode};
%     nChunks = size(chunkInfo, 1);
%     
%     % If this channel was disabled, skip it
%     if nChunks == 0
%       continue;
%     end
%     
%     fprintf('%d ', electrode);
%     
%     %% Load all chunks for this channel
%     % elecWaveforms will hold all the chunks to be combined in one step at
%     % the end
%     elecWaveforms = cell(1, nChunks);
%     
%     % Load chunks, delete chunks files
%     for chunk = 1:nChunks
%       chunkFilename = makeChunkFilename(chunkDir, lett, electrode, chunkInfo(chunk, 1), chunkInfo(chunk, 2));
%       loadVar = load(chunkFilename);
%       elecWaveforms{chunk} = loadVar.waveforms;
%       
%       delete(chunkFilename);
%     end
%     
%     %% Create full waveforms structure for this electrode
%     waveforms = elecWaveforms{1};
%     waveforms.waves = cellfun(@(x)x.waves, elecWaveforms, 'UniformOutput', false);
%     waveforms.waves = [waveforms.waves{:}];
%     waveforms.units = cellfun(@(x)x.units, elecWaveforms, 'UniformOutput', false);
%     waveforms.units = [waveforms.units{:}];
%     waveforms.spikeTimes = cellfun(@(x)x.spikeTimes, elecWaveforms, 'UniformOutput', false);
%     waveforms.spikeTimes = [waveforms.spikeTimes{:}];
%     waveforms.ratings.epoch = [1 length(waveforms.spikeTimes)];
%     waveforms.sourceFiles = chunkFilenamesToSourceFiles(elecWaveforms);
%     
%     %% Remove artifacts
%     % Note that artifactThresh is specified in mV, need to convert to uV
%     if ~isinf(artifactThresh)
%       maxes = max(abs(waveforms.waves));
%       artifacts = (maxes > abs(artifactThresh) * 1000);
%       
%       waveforms = removeBadSpikes(waveforms, artifacts);
%     end
%     
%     %% If there were refractory violations caused by a Cerebus bug, fix
%     % We need to do a two-pass repair: first we'll remove all waveforms
%     % that don't cross threshold at the right time, then we'll remove the
%     % small number of potentially remaining lockout violators.
%     if flagrantLockoutViolations
%       %% Find threshold value and sample
%       [thresh, threshi] = findThreshold(waveforms.waves);
%       
%       %% Remove spikes that don't cross threshold at the right time
%       if ~isempty(waveforms.waves)
%         nonCrossers = (waveforms.waves(threshi-1, :) <= thresh | waveforms.waves(threshi, :) > thresh);
%         fracNonCrossers(end+1) = sum(nonCrossers) / size(waveforms.waves, 2); %#ok<AGROW>
%         
%         waveforms = removeBadSpikes(waveforms, nonCrossers);
%       end
%       
%       %% Find remaining lockout violations
%       if ~isempty(waveforms.waves)
%         diffs = diff(waveforms.spikeTimes);
%         violations = [false (diffs > 0 & diffs < lockout)];
%         
%         waveforms = removeBadSpikes(waveforms, violations);
%         
%         remLockoutViol(end+1) = sum(violations); %#ok<AGROW>
%       else
%         remLockoutViol(end+1) = 0; %#ok<AGROW>
%       end
%     end
%     
%     
%     % Check that the channel didn't end up with absolutely nothing in it
%     if ~isempty(waveforms.waves)
%       %% Save channel
%       waveFilename = makeWaveFilename(lett, electrode);
%       save(fullfile(outDir, waveFilename), 'waveforms', '-v7');
%       
%       %% Generate this channel's worth of previews and sorts
%       usedChs = usedChs + 1;
%       [previews(usedChs), sorts(usedChs)] = generateSortsAndPreviews(waveforms, maxUnits, lockout); %#ok<AGROW,NASGU>
%     end
%       
%     if mod(electrode, 20) == 0, fprintf('\n'); end
%   end
%   fprintf('\n');
% end
% 
% 
% %% Attempt to delete temp chunk directory
% 
% if ~strcmp(chunkDir, outDir)
%   status = rmdir(chunkDir);
%   
%   if ~status
%     fprintf('Warning: could not delete the temporary extraction directory.\n');
%   end
% end
% 
% 
% %% If relevant, report .nev bug repairs
% if flagrantLockoutViolations
%   fprintf('\nPercent of waveforms that failed to cross threshold at the right time, by channel:\n');
%   blockPrintNumbers(100 * fracNonCrossers, 20);
%   fprintf('Mean: %0.3f\n', mean(100 * fracNonCrossers));
%   fprintf('Total number of remaining waveforms involved in a residual lockout violation, by channel:\n');
%   blockPrintNumbers(remLockoutViol, 20);
%   fprintf('Mean: %0.3f%%\n', mean(remLockoutViol));
% end
%   
% 
% %% Save previews and sorts
% save(fullfile(outDir, 'previews.mat'), 'previews', '-v6');
% save(fullfile(outDir, 'sorts.mat'), 'sorts', '-v6');
% 
% fprintf('Done.\n');
% 
% 
% 
% 

% %%
% function chunkFilename = makeChunkFilename(outDir, lett, electrode, fileNum, chunki)
% name = makeWaveFilename(lett, electrode);
% stem = name(1:end-4);
% chunkFilename = fullfile(outDir, sprintf('%s_%d_%d.mat', stem, fileNum, chunki));
% 

% 
% %%
% function [lockout, lockoutViolations, flagrantLockoutViolations] = ...
%   checkForLockoutViolations(spikeTimes, allElectrodeNums, snippetLockout, lockoutViolations, flagrantLockoutViolations, RRRlockout)
% % There are three possibilities. One is that the lockout period is
% % determined normally by the snippet length. A second is that data was
% % recorded with a buggy version of the Cerebus software, which will cause
% % flagrant lockout violations (e.g., 100 us ISIs). A third is that the data
% % was recorded with RRR, which has a shorter lockout than normal Cerebus.
% 
% % ELB_NOTE: Lockout -- Deadtime
% 
% if ~flagrantLockoutViolations
%   % Check for violations of the lockout period calculated by snippet
%   % length. If we find any, assume the data was collected with RRR,
%   % which uses a lockout shorter than the snippet length.
%   lockout = snippetLockout;
%   
%   electrodes = unique(allElectrodeNums);
%   
%   for electrode = electrodes
%     diffs = diff(spikeTimes(allElectrodeNums == electrode));
%     
%     if any(diffs >= 0 & diffs < RRRlockout)
%       fprintf('Found violations of lockout period; assuming lockout period is snippet length\n');
%       flagrantLockoutViolations = 1;
%       lockout = snippetLockout;
%       break;
%     elseif lockoutViolations || any(diffs >= 0 & diffs < snippetLockout)
%       lockoutViolations = 1;
%       lockout = RRRlockout;
%     else
%       lockout = snippetLockout;
%     end
%   end
% end
% 
% if flagrantLockoutViolations
%   lockout = snippetLockout;
% end
% 
% 
% 
% %%
% function [thresh, threshi] = findThreshold(waves)
% % Find threshold value and time. We can't use the value in the spike
% % headers, because Blackrock changed between 0 and 1 indexing for the
% % number of pre-threshold samples between version 4 and 6 of the
% % software/firmware (and possibly introduced another weird offset).
% 
% % Note - This only returns the negative threshold.  There may be different
% % negative and positive thresholds
% 
% % The hack we'll use to find the threshold-crossing sample is to find the
% % sample at which the most waveforms have just gone negatively
% diffs = diff(waves);
% nNegs = sum(diffs < 0, 2);
% [junk, threshi] = max(nNegs); %#ok<ASGLU>
% threshi = threshi + 1;        % Need to correct for off-by-1 of diff
% 
% % We'll try different strategies depending on how many waveforms are
% % present
% nWaves = size(waves, 2);
% if nWaves < 20
%   thresh = -1e-9;
%   quantiles = NaN;
% else
%   % To find the actual threshold value, we'll take the percentiles of
%   % post-threshold values from 95th to 99.9th (to exclude the few
%   % violators). We'll try to find two percentiles 0.2% apart that agree,
%   % and take the highest ones that do. May have to use coarser-grained
%   % quantiles if the number of waveforms is small.
%   if nWaves < 100
%     quantiles = 0.999:-1/nWaves:0.5;
%     criterion = 5;
%   elseif nWaves < 1000
%     quantiles = 0.999:-1/nWaves:0.5;
%     criterion = 1;
%   else
%     quantiles = (0.999:-0.001:0.950);
%     criterion = 1e-9;
%   end
%   wavesQ = quantile_mksort(waves(threshi, :), quantiles);
%   thresh = NaN;
%   for q = 1:length(quantiles) - 2
%     if abs(wavesQ(q) - wavesQ(q + 2)) < criterion
%       thresh = wavesQ(q);
%       break;
%     end
%   end
% end
% 
% if ~isnan(thresh)
%   if ~isnan(quantiles)
%     fprintf('Using threshold at %0.1f percentile\n', quantiles(q) * 100);
%   else
%     fprintf('So few threshold crossings present that the cleaning threshold was chosen arbitrarily at -1e-9\n');
%   end
% else
%   fprintf('WARNING: Could not find a valid threshold between 95th and 99.9th percentile of post-threshold values\n');
%   fprintf('Since a consistent threshold could not be found, arbitrarily using -1e-9. Invalid waveforms likely remain.\n');
%   thresh = -1e-9;
% end
% 
% 
% 
% %%
% function sourceFiles = chunkFilenamesToSourceFiles(elecWaveforms)
% % Take the elecWaveforms cell array, and use the spike arrays and
% % fileSource info contained therein to produce a cell array of filenames
% % and an nFile x 2 table of spike indices telling which file spikes came
% % from.
% 
% sourceFiles.filenames{1} = elecWaveforms{1}.sourceFiles;
% sourceFiles.spikesByFile = [1 0];
% 
% filei = 1;
% 
% for ch = 1:length(elecWaveforms)
%   chunk = elecWaveforms{ch};
%   
%   if ~strcmp(chunk.sourceFiles, sourceFiles.filenames{filei})
%     filei = filei + 1;
%     sourceFiles.filenames{filei} = chunk.sourceFiles;
%     sourceFiles.spikesByFile(filei, :) = sourceFiles.spikesByFile(filei - 1, 2) + [1, length(chunk.spikeTimes)];
%   else
%     sourceFiles.spikesByFile(filei, 2) = sourceFiles.spikesByFile(filei, 2) + length(chunk.spikeTimes);
%   end
% end
% 
% % Now, check for rows where no spikes were present. They will have a last
% % spike less than their first spike, but subsequent rows will be ok.
% % Replace those bad rows with NaNs.
% badRows = find(diff(sourceFiles.spikesByFile, 1, 2) < 0);
% if ~isempty(badRows)
%   sourceFiles.spikesByFile(badRows, :) = [NaN NaN];
% end
% 
% 
% %%
% function waveforms = removeBadSpikes(waveforms, badSpikes)
% 
% waveforms.waves(:, badSpikes) = [];
% waveforms.units(:, badSpikes) = [];
% waveforms.spikeTimes(:, badSpikes) = [];
% waveforms.ratings.epoch = [1 length(waveforms.spikeTimes)];
% waveforms.sourceFiles = correctSourceFilesTable(waveforms.sourceFiles, badSpikes);
% 
% 
% 
% 
% %%
% function sourceFiles = correctSourceFilesTable(sourceFiles, badSpikes)
% % Take a sourceFiles struct and a logical array of bad spikes that will be
% % removed, and fix the table values.
% 
% badAtGivenPos = cumsum(badSpikes);
% for filei = 1:size(sourceFiles.spikesByFile, 1)
%   val1 = sourceFiles.spikesByFile(filei, 1);
%   if ~isnan(val1)
%     % Have to add back the badSpikes value at first element, since
%     % otherwise we may decrement it (and we never should, or else it will
%     % be on top of the last spike of the previous file)
%     sourceFiles.spikesByFile(filei, 1) = val1 - badAtGivenPos(val1) + badSpikes(val1);
%     val2 = sourceFiles.spikesByFile(filei, 2);
%     sourceFiles.spikesByFile(filei, 2) = val2 - badAtGivenPos(val2);
%   end
% end
