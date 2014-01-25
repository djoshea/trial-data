classdef SpikeRaster < handle & matlab.mixin.Copyable

    properties(SetAccess=protected)
        spikes % cell array over trials containing spike times (must be sorted!)
        spikesExcludingPad
        waveforms   % cell array over trials containing waveforms
        meta        % array over trials with arbitrary metadata

        conditionInfo   % string or integer indicating which condition 

        alignInfo % instance of alignInfo from which this was built
        tUnits = 'ms';
        tUnitsPerSec = 1000;
 
        psthMinTrials = 5; % time bins with too few trials are marked as NaN in psth

        % maps (R, unit) --> spikeTimesCell of size nTrials x 1
        getSpikesFn = @SpikeRaster.defaultGetSpikesFn;

        tBinWidth
        
        % these manually specify the time window and override tMin and tMax
        % when present. They were added to deal with issues with changing 
        % time windows when resampling trials occurs
        tMinManual
        tMaxManual
        
        padWindow % padding added to the windows in each trial to facilitate artifact free smoothing
    end

    properties
        name = '';
        info % arbitrary field for user storage
        
        preserveTimeWindowWhenResampling = true;
        
        % if true, tMin and tMax will be set to the widest 
        % setting this calls clearCache()
        useWidestCommonValidTimeWindow = false;
        
        % see set.spikeFilter 
        spikeFilter 
    end

    properties(SetAccess=protected)
        valid
    end

    properties(SetAccess=protected, Hidden) % cached data, cleared using clearCaches
        cachedSmoothedRatesByTrial
    end

    properties(Dependent)
        nTrials
        nTrialsValid
        nSpikesPerTrial
        nSpikesTotal
        hasWaveforms
        hasMeta

        minTrialsByCondition % minimum trial count across all conditions
        minTrialsByNonEmptyCondition % same as above only ignore conditions with no trials
        nTrialsByCondition
        nConditions
        idxByCondition
        conditionNames
        conditionAppearance
        
        conditionByTrial

        % aligned start and stop times for each trial
        tMinByTrial
        tMaxByTrial

        % widest valid time window per-condition
        tMinByCondition
        tMaxByCondition

        % widest valid time window across ALL trials, overridden by setting
        % tMinManual and tMaxManual above
        tMin
        tMax

        tBins % used for all binned counts, psths, etc.
        nTBins

        tLengthByTrial
    end

    methods % constructor 
        function obj = SpikeRaster(R, unit, varargin)
            alignDescriptor = [];
            alignInfo = [];
            conditionInfo = [];
            conditionDescriptor = [];
            name = '';
            tBinWidth = 1;
            spikeFilter = [];
            def = structargs(varargin);

            obj.tBinWidth = def.tBinWidth;

            obj.spikeFilter = def.spikeFilter;
            if isempty(obj.spikeFilter)
                obj.spikeFilter = GaussianSpikeFilter();
            end

            alignDescriptor = def.alignDescriptor;
            alignInfo = def.alignInfo;

            % copy this to avoid interference by modifying the condition info externally
            if isempty(def.conditionInfo) && isempty(def.conditionDescriptor)
                error('Please provide conditionDescriptor or conditionInfo argument to describe how to structure comparisons'); 
            end
            
            if isempty(def.conditionInfo)
                def.conditionInfo = ConditionInfo.fromConditionDescriptor(def.conditionDescriptor, R);
            end
            
            obj.conditionInfo = def.conditionInfo;
            
            if ~isempty(def.name)
                obj.name = def.name;
            else
                try 
                    desc = getRStructDescription(R);
                    obj.name = sprintf('%s Unit %s', desc, unit);
                catch
                    obj.name = sprintf('Unit %s', unit);
                end
            end

            rawSpikes = obj.getSpikesFn(R, unit);
            % grab extra data to accommodate SpikeFilter
            obj.padWindow = [obj.spikeFilter.preWindow obj.spikeFilter.postWindow];
            
            if isempty(alignDescriptor) && isempty(alignInfo)
                error('Please provide alignDescriptor or alignInfo argument to describe how to align and segment spikes'); 
            end
            if isempty(alignInfo)
                if ~isa(alignDescriptor, 'AlignInfo');
                    obj.alignInfo = AlignInfo.fromAlignDescriptor(alignDescriptor);
                else
                    obj.alignInfo = alignDescriptor;
                end
                obj.alignInfo = obj.alignInfo.bind(R);
            else
                assert(isa(alignInfo, 'AlignInfo'), 'alignInfo argument must be an AlignInfo instance');
                obj.alignInfo = alignInfo;
            end
            obj.alignInfo = obj.alignInfo.pad(obj.padWindow);
            
            %obj.alignTimeInfo = obj.alignDescriptor.getTimeInfo(R, obj.padWindow);
           % [obj.spikes] = obj.alignDescriptor.getAlignedTimes(obj.alignTimeInfo, rawSpikes);
            obj.spikes = obj.alignInfo.getAlignedTimes(rawSpikes, true);
            
            obj.spikesExcludingPad = obj.getSpikesExcludingPad();

            emptyMask = cellfun(@isempty, obj.spikes);
            obj.spikes(emptyMask) = {[]};
            
            % update the validity mask in the time info in order for draw time axis to work correctly
%             ciValid = makecol(obj.conditionInfo.valid);
%             if ~any(ciValid)
%                 warning('No valid trials found via ConditionInfo');
%             end
%             tiValid = makecol([obj.alignInfo.valid]);
%             if ~any(tiValid)
%                 warning('No valid trials found via AlignDescriptor');
%             end
%             valid = ciValid & tiValid;
%             if ~any(valid)
%                 warning('No valid trials found via ConditionInfo and AlignDescriptor');
%             end
%             
%             obj.valid = valid;
            
            % mark as invalid trials which do not meet the alignment descriptor's inclusion criteria
            % this ensures they will not be included when we call conditionInfo.getGroups()
            %obj.alignInfo = obj.alignInfo.markInvalid(~valid);
            %obj.conditionInfo = obj.conditionInfo.markInvalid(~obj.valid);
        end
    end
    
    methods(Static)
        function timesCell = defaultGetSpikesFn(data, unitName)
            if(~isa(data, 'TrialData'))
                error('Please set .getSpikesFn to provide access to raw spike times');
            end
            
            timesCell = data.getRawSpikeTimes(unitName);
        end
    end

    methods(Access=protected) % copyElement for Copyable
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            cpObj.conditionInfo = obj.conditionInfo.copy();
            cpObj.spikeFilter = obj.spikeFilter.copy();
        end
    end

    methods % dependent properties
        function n = get.nTrials(obj)
            n = length(obj.spikes);
        end
        
        function n = get.nTrialsValid(obj)
            n = nnz(obj.valid);
        end

        function spikes = getSpikesExcludingPad(obj)
            spikesOrig = obj.spikes;
            spikes = cell(obj.nTrials, 1);
            
            for i = 1:obj.nTrials
                ti = obj.alignInfo.timeInfo(i);
                if ti.valid
                    spikes{i} = spikesOrig{i}(spikesOrig{i} >= (ti.start - ti.zero) & spikesOrig{i} <= (ti.stop - ti.zero));
                end
            end
        end 
        
        function counts = get.nSpikesPerTrial(obj)
            counts = cellfun(@length, obj.spikesExcludingPad);
        end
    
        function n = get.nSpikesTotal(obj)
           n = sum(obj.nSpikesPerTrial);
        end

        function tf = get.hasWaveforms(obj)
            tf = ~isempty(obj.waveforms);
        end

        function tf = get.hasMeta(obj)
            tf = ~isempty(obj.meta);
        end

        function set.useWidestCommonValidTimeWindow(obj, val)
            if obj.useWidestCommonValidTimeWindow ~= val
                obj.useWidestCommonValidTimeWindow = val;
                obj.clearCaches();
            end
        end
        
        function set.spikeFilter(obj, val)
            % auto-clear the psths
            if ~isequal(obj.spikeFilter, val)
                obj.spikeFilter = val;
                obj.clearCaches();
            end
        end

        function tMinByTrial = get.tMinByTrial(obj)
            tMinByTrial = makecol([obj.alignInfo.timeInfo.start] - [obj.alignInfo.timeInfo.zero]);
        end

        function tMaxByTrial = get.tMaxByTrial(obj)
            tMaxByTrial = makecol([obj.alignInfo.timeInfo.stop] - [obj.alignInfo.timeInfo.zero]);
        end

        function tLengthByTrial = get.tLengthByTrial(obj)
            tLengthByTrial = obj.tMaxByTrial - obj.tMinByTrial;
        end

        function tMinByCondition = get.tMinByCondition(obj)
            tMinByCondition = nan(obj.conditionInfo.conditionsSize);
            for iC = 1:obj.nConditions
                if obj.nTrialsByCondition(iC) > 0
                    tMinByCondition(iC) = nanmax(obj.tMinByTrial(obj.idxByCondition{iC}));
                else
                    tMinByCondition(iC) = NaN;
                end
            end
        end
    
        function tMaxByCondition = get.tMaxByCondition(obj)
            tMaxByCondition = nan(obj.conditionInfo.conditionsSize);
            for iC = 1:obj.nConditions
                if obj.nTrialsByCondition(iC) > 0
                    tMaxByCondition(iC) = nanmin(obj.tMaxByTrial(obj.idxByCondition{iC}));
                else
                    tMaxByCondition(iC) = NaN;
                end
            end
        end

        function tMin = get.tMin(obj)
            if isempty(obj.tMinManual)
                if obj.useWidestCommonValidTimeWindow
                    tMin = nanmax(obj.tMinByTrial(obj.valid));
                else
                    tMin= nanmin(obj.tMinByTrial(obj.valid));
                end
            else
                tMin = obj.tMinManual;
            end
            if isempty(tMin)
                tMin = NaN;
            end
        end
    
        function tMax = get.tMax(obj)
            if isempty(obj.tMaxManual)
                if obj.useWidestCommonValidTimeWindow
                    tMax = nanmin(obj.tMaxByTrial(obj.valid));
                else
                    tMax = nanmax(obj.tMaxByTrial(obj.valid));
                end
            else
                tMax = obj.tMaxManual;
            end
            if isempty(tMax)
                tMax = NaN;
            end
        end

        function tBins = get.tBins(obj)
            tBins = makecol(obj.tMin : obj.tBinWidth : obj.tMax);
        end

        function nTBins = get.nTBins(obj)
            nTBins = length(obj.tBins);
        end

        function nConditions = get.nConditions(obj)
            nConditions = numel(obj.idxByCondition);
        end

        function nTrialsByCondition = get.nTrialsByCondition(obj)
            nTrialsByCondition = cellfun(@length, obj.idxByCondition);
        end
        
        function minTrialsByCondition = get.minTrialsByCondition(obj)
            minTrialsByCondition = min(obj.nTrialsByCondition(:));
        end
        
        function minTrialsByCondition = get.minTrialsByNonEmptyCondition(obj)
            nTrials = obj.nTrialsByCondition(:);
            nTrials = nTrials(nTrials > 0);
            if isempty(nTrials)
                minTrialsByCondition = 0;
            else
                minTrialsByCondition = min(nTrials);
            end
        end

        function idxByCondition = get.idxByCondition(obj)
            idxByCondition = obj.conditionInfo.listByCondition;
        end
        
        function conditionByTrial = get.conditionByTrial(obj)
            conditionByTrial = obj.conditionInfo.conditionIdx;
        end

        function conditionNames = get.conditionNames(obj)
            conditionNames = obj.conditionInfo.names;
        end

        function conditionAppearance = get.conditionAppearance(obj)
            conditionAppearance = obj.conditionInfo.appearances;
        end

        function spikesByCond = getSpikesForCondition(obj, condIdx)
%             if nargin < 2
%                condIdx = true(obj.conditionInfo.conditionsSize);
%             end
%             if islogical(condIdx)
%                 fullCondInds = containingLinearInds(obj.conditionInfo.conditionsSize);
%                 condIdx = fullCondIdx(condIdx);
%             end
%                    
            % return a cell of cells, each entry containing spikes(idxForCond) for each condition
            
            
            spikesByCond = obj.spikesExcludingPad(obj.idxByCondition{condIdx});
        end

    end

    methods
        function clearCaches(obj)
            obj.cachedSmoothedRatesByTrial = [];
        end

        % internally reorder the trials using selector mask 
        function resampleTrials(obj, mask, varargin)
            p = inputParser;
            p.addParamValue('resampleConditionInfo', true, @islogical)
            p.parse(varargin{:});
            
            if obj.preserveTimeWindowWhenResampling
                % the time window will stay fixed at what it is now
                % to prevent the resampling from growing the validity
                % window and making things confusing for consumers of the
                % PSTHs.
                obj.tMinManual = obj.tMin;
                obj.tMaxManual = obj.tMax;
            end

            obj.spikes = obj.spikes(mask);
            obj.spikesExcludingPad = obj.spikesExcludingPad(mask);
            if ~isempty(obj.waveforms)
                obj.waveforms = obj.waveforms(mask);
            end
            if ~isempty(obj.meta)
                obj.meta = obj.meta(mask);
            end

            if p.Results.resampleConditionInfo
                obj.conditionInfo.resampleTrials(mask);
            end

            if ~isempty(obj.alignInfo)
                obj.alignInfo = obj.alignInfo.selectTrials(mask);
            end

            if ~isempty(obj.cachedSmoothedRatesByTrial)
                obj.cachedSmoothedRatesByTrial = obj.cachedSmoothedRatesByTrial(mask,:);
            end
            
            obj.valid = true(nnz(mask), 1);
        end

        function sr = filteredByAttribute(obj, attr, value)
            sr = obj.copy();
            sr.clearCaches();
            sr.conditionInfo = sr.conditionInfo.filteredByAttribute(attr, value);
        end
    end

    methods % filtering trials / conditions
        function filterValidTrials(obj)
            % drop all trials which are invalid according to obj.valid
            obj.resampleTrials(obj.valid);
        end
    end

    methods % Alignment time info 
        % build a cell array in the size of conditionInfo.getGroups containing the
        % aligned time info struct for all trials in that condition
        function alignTimeInfo = getAlignTimeInfoByCondition(obj, varargin)
            idxByCondition = obj.conditionInfo.listByCondition;
            alignTimeInfo = cell(size(idxByCondition));

            for i = 1:numel(alignTimeInfo)
                alignTimeInfo{i} = makecol(obj.alignInfo.timeInfo(idxByCondition{i}));
            end
        end

        function alignTimeInfo = getAlignTimeInfoForCondition(obj, conditionInd, varargin)
            idxByCondition = obj.conditionInfo.listByCondition;
            alignTimeInfo = makecol(obj.alignInfo.timeInfo(idxByCondition{conditionInd}));
        end
    end

    methods % binned spike counts and smoothed rates

        function [counts timeBins] = getBinnedCountsByTrial(obj, varargin)
            % for each trial, counts the number of spikes within each time bin
            % returns a matrix nTrials x nTimes containing the number of spikes observed
            % or NaN if outside the trials time window
            assignargs(varargin);

            timeBins = obj.tBins;
            binWidth = obj.tBinWidth;
            timeBinEdges = [timeBins; obj.tMax];
            nTimeBins = length(timeBins);

            tMinByTrial = obj.tMinByTrial;
            tMaxByTrial = obj.tMaxByTrial;
            counts = nan(obj.nTrials, nTimeBins);

            for iTrial = 1:obj.nTrials
                spikes = obj.spikes{iTrial};
                validBinsThisTrial = find(timeBins >= tMinByTrial(iTrial) & ...
                                timeBins + binWidth <= tMaxByTrial(iTrial));
                if isempty(spikes)
                    countsThisTrial = zeros(1, length(timeBinEdges)-1);
                else
                    countsThisTrial =  histc(spikes, timeBinEdges);
                end
                counts(iTrial, validBinsThisTrial) = countsThisTrial(validBinsThisTrial);
            end
        end

        function [countsByCond timeBins] = getBinnedCountsByCondition(obj, varargin)
            idxByCondition = [];
            assignargs(varargin);
            
           [counts timeBins] = obj.getBinnedCountsByTrial(varargin{:}); 

            if isempty(idxByCondition)
                idxByCondition = obj.conditionInfo.getGroups();
            end
            nGroups = numel(idxByCondition);
            countsByCond = cell(size(idxByCondition));
            for iCond = 1:nGroups
                countsByCond{iCond} = counts(idxByCondition{iCond}, :); 
            end
       end

        function [smoothedCounts timeBins] = getSmoothedRatesByTrial(obj, varargin)
            assignargs(varargin);

            if isempty(obj.cachedSmoothedRatesByTrial)
                %counts = obj.getBinnedCountsByTrial();
                % TODO probably need to make sure this handles NaNs at the tails correctly
                %obj.cachedSmoothedRatesByTrial = FilterSpikes(obj.filterWidth, counts, 'causal', obj.filterCausal) * ...
                %    obj.tUnitsPerSec; 
                obj.cachedSmoothedRatesByTrial = obj.spikeFilter.filterSpikeTrains(obj.spikes, [obj.tMin obj.tMax]);
                
                % go through and mark as NaN any time outside each trial's
                % valid window, in case tMin / tMax are larger than that
                timeBins = obj.tBins;
                tMinByTrial = obj.tMinByTrial;
                tMaxByTrial = obj.tMaxByTrial;
                for iTrial = 1:obj.nTrials
                    mask = false(length(timeBins), 1);
                    mask(timeBins < tMinByTrial(iTrial)) = true;
                    mask(timeBins > tMaxByTrial(iTrial)) = true;
                    obj.cachedSmoothedRatesByTrial(iTrial, mask) = NaN;
                end
            end

            smoothedCounts = obj.cachedSmoothedRatesByTrial;
            timeBins = obj.tBins;
        end

        function [smoothedRatesByCond timeBins] = getSmoothedRatesByCondition(obj, varargin)
            idxByCond = [];
            assignargs(varargin);

            [smoothedRates timeBins] = obj.getSmoothedRatesByTrial(varargin{:});

            if isempty(idxByCond)
                idxByCond = obj.conditionInfo.getGroups();
            end

            nGroups = numel(idxByCond);
            smoothedRatesByCond = cell(obj.conditionInfo.conditionsSize);
            for iCond = 1:nGroups
                smoothedRatesByCond{iCond} = smoothedRates(idxByCond{iCond}, :); 
            end
        end

        function [rates timeBins] = getRatesByTrial(obj, varargin)
            smooth = true;
            assignargs(varargin);

            if smooth
                [rates timeBins] = obj.getSmoothedRatesByTrial();
            else
                [counts timeBins] = obj.getBinnedCountsByTrial();
                rates = counts ./ (obj.tBinWidth / obj.tUnitsPerSec);
            end
            
            timeBins = makecol(timeBins);
        end

        function [psthByCondition semByCondition timeBins] = getPSTHByCondition(obj, varargin)
            idxByCondition = [];
            smooth = true;
            computeSem = true;
            def = assignargs(varargin);

            [counts timeBins] = obj.getRatesByTrial(def);

            if isempty(idxByCondition)
                idxByCondition = obj.conditionInfo.listByCondition;
            end
            nConds = numel(idxByCondition);
            psthByCondition = zeros(nConds, length(timeBins));
            if computeSem
                semByCondition = zeros(nConds, length(timeBins));
            else
                semByCondition = [];
            end

            for iCond = 1:nConds    
                idx = idxByCondition{iCond};
                % TODO probably need to make sure this deals with nans correctly
                psthByCondition(iCond, :) = nanMeanMinCount(counts(idx, :), 1, obj.psthMinTrials); 
                if computeSem
                    semByCondition(iCond, :) = nansem(counts(idx, :), 1);
                    semByCondition(iCond, isnan(psthByCondition(iCond,:))) = NaN;
                end
            end
        end

    end

    methods % PSTH distance along comparison axis of conditionInfo

        function [distanceByTime meanDistance timeBins] = getPSTHDistance(obj, compareAcross, varargin)
            % compare PSTHs along some axis of the condition info and compare the distance between them 
            % over time, in units of Hz
            % distanceByTime nComparisons x nTimeBins
            % distanceMean nComparisons x 1 
            
            idxCompare = [];
            namesCompare = [];
            rectify = false;
            smooth = true;
            assignargs(varargin);
           
            % default is from condition info, but allow specification for the sake of doing bootstraps
            if isempty(idxCompare) || isempty(namesCompare);
                [idxCompare namesCompare appearanceCompare] = ...
                    obj.conditionInfo.getComparisonAxis(compareAcross);
            end
            
            nCompare = numel(idxCompare);
            comparisonNames = cell(nCompare, 1);
            meanDistance = nan(nCompare, 1);
            timeBins = obj.tBins;
            distanceByTime = nan(nCompare, length(timeBins));
            
            for iCompare = 1:nCompare
                idxThisComparison = idxCompare{iCompare};
                if isempty(idxThisComparison)
                    continue;
                end
                nTrials = cellfun(@length, idxThisComparison);
                
                if all(nTrials > 0)
                    assert(length(idxThisComparison) == 2, 'Comparison axis must span exactly two elements');

                    [psthBoth, ~, timeBins] = obj.getPSTHByCondition('idxByCondition', idxThisComparison, 'smooth', smooth); 
                    delta = psthBoth(2,:) - psthBoth(1,:);

                    if rectify
                        delta = abs(delta);
                    end

                    timePerBin = obj.tBinWidth / obj.tUnitsPerSec;
                    useDelta = removenan(delta);
                    if ~isempty(useDelta)
                        meanDistance(iCompare) = mean(useDelta);
                    else
                        meanDistance(iCompare) = NaN;
                    end

                    if ~exist('distanceByTime', 'var')
                        distanceByTime = nan(nCompare, length(timeBins));
                    end
                    distanceByTime(iCompare, :) = delta;

                   % comparisonNames{iCompare} = sprintf('%s vs %s', ...
                   %     namesCompare{iCompare}{2}, namesCompare{iCompare}{1});
                end       
            end
        end

        function [distanceByTimeFloor meanDistanceFloor distanceByTimeMat meanDistanceMat] = ...
                getPSTHDistanceShuffled(obj, compareAcross, varargin)
            nShuffles = 100;
            p = 0.05;
            def = assignargs(varargin);

            % we pass this to getPSTH to avoid unnecessary repeated calculation
            [smoothedCounts timeBins] = obj.getSmoothedRatesByTrial();

            % these are cells over shuffles of cells over conditions of cells over the comparison axis (2 elements)
            [idxCompareShuffles namesCompare] = ...
                obj.conditionInfo.getComparisonAxisShuffled(compareAcross, 'nShuffles', nShuffles);
           
            [distanceByTime meanDistanceShuffles] = deal(cell(nShuffles,1));
            for iShuffle = 1:nShuffles 
                % pull out the cells over conditions of cells over the comparison axis here
                idxCompare = idxCompareShuffles{iShuffle};

                [distanceByTime{iShuffle} meanDistanceShuffles{iShuffle} timeBins] = ...
                    obj.getPSTHDistance(compareAcross, def, 'idxCompare', idxCompare, ...
                    'namesCompare', namesCompare);
            end

            % distanceByTime is a nShuffles x cell array of arrays of size nCompare x nTimeBins
            % distanceByTimeMat is a nCompare x nTimeBins x nShuffle array 
            distanceByTimeMat = cell2mat(shiftdim(distanceByTime,-2));

            % and we compute (1-p) percentile of the data along the shuffle axis as the floor per comparison
            distanceByTimeFloor = prctile(abs(distanceByTimeMat), 100*(1-p), 3);

            % meanDistanceMat is nCompare x nShuffles
            meanDistanceMat = cell2mat(shiftdim(meanDistanceShuffles, -1));
            meanDistanceFloor = prctile(abs(meanDistanceMat), 100*(1-p), 2);
        end

        % similar to shuffled except resamples trial labels rather than shuffling them
        % returns upper and lower confidence intervals in distanceByTimeIntervals
        % specify param 'rectify', true if you want absolute values
        %
        % nCompare is the number of unique values for the other attributes
        %
        % distanceByTimeIntervals : nCompare x nTimeBins x 2 (low/high) confidence intervals
        % meanDistanceInterval : nCompare x 2 mean confidence intervals
        % distanceByTimeMat : nCompare x nTimeBins x nResample (default 200) raw bootstrap resampled distance traces
        % meanDistanceMat : nCompare x nResample raw mean resampled 
        function [distanceByTimeIntervals meanDistanceInterval distanceByTimeMat meanDistanceMat] = ...
                getPSTHDistanceResampled(obj, compareAcross, varargin)
            nResample = 200;
            p = 0.05;
            def = assignargs(varargin);

            % these are cells over resamples of cells over conditions of cells over the comparison axis (2 elements)
            [idxCompareResample namesCompare] = ...
                obj.conditionInfo.getComparisonAxisResampled(compareAcross, 'nResample', nResample);
           
            [distanceByTime meanDistanceResampled] = deal(cell(nResample,1));
            for iResample = 1:nResample 
                % pull out the cells over conditions of cells over the comparison axis here
                idxCompare = idxCompareResample{iResample};

                [distanceByTime{iResample} meanDistanceResampled{iResample} timeBins] = ...
                    obj.getPSTHDistance(compareAcross, def, 'idxCompare', idxCompare, ...
                    'namesCompare', namesCompare);
            end

            % distanceByTime is a nResample x cell array of arrays of size nCompare x nTimeBins
            % distanceByTimeMat is a nCompare x nTimeBins x nResample array 
            distanceByTimeMat = cell2mat(shiftdim(distanceByTime,-2));

            % and we compute p/2 percentile edges of the data along the shuffle axis as the floor per comparison
            szMat = size(distanceByTimeMat);
            distanceByTimeIntervals = nan(szMat(1), szMat(2), 2);
            distanceByTimeIntervals(:, :, 1) = prctile(distanceByTimeMat, 100*(p/2), 3);
            distanceByTimeIntervals(:, :, 2) = prctile(distanceByTimeMat, 100*(1-p/2), 3);

            % meanDistanceMat is nCompare x nShuffles
            meanDistanceMat = cell2mat(shiftdim(meanDistanceResampled, -1));
            meanDistanceInterval = nan(size(meanDistanceMat, 1), 2);
            meanDistanceInterval(:,1) = prctile(meanDistanceMat, 100*(p/2), 2);
            meanDistanceInterval(:,2) = prctile(meanDistanceMat, 100*(1-p/2), 2);
        end

        % similar to getPSTHDistanceResampled except resamples all trial labels from
        % the first element of the comparison axis, such that we end up with a measure
        % of the measurement noise if all trials were drawn from the first trial type
        % along the comparison axis rather than from both
        %
        % returns upper and lower confidence intervals in distanceByTimeIntervals
        % specify param 'rectify', true if you want absolute values
        %
        % nCompare is the number of unique values for the other attributes
        %
        % distanceByTimeIntervals : nCompare x nTimeBins x 2 (low/high) confidence intervals
        % meanDistanceInterval : nCompare x 2 mean confidence intervals
        % distanceByTimeMat : nCompare x nTimeBins x nResample (default 200) raw bootstrap resampled distance traces
        % meanDistanceMat : nCompare x nResample raw mean resampled 
        function [distanceByTimeIntervals meanDistanceInterval distanceByTimeMat meanDistanceMat] = ...
                getPSTHDistanceResampledFromSame(obj, compareAcross, varargin)
            nResample = 200;
            p = 0.05;
            def = assignargs(varargin);

            % these are cells over resamples of cells over conditions of cells over the comparison axis (2 elements)
            [idxCompareResample namesCompare] = ...
                obj.conditionInfo.getComparisonAxisResampledFromSame(compareAcross, 'resampleFrom', 1, 'nResample', nResample);
           
            [distanceByTime meanDistanceResampled] = deal(cell(nResample,1));
            for iResample = 1:nResample 
                % pull out the cells over conditions of cells over the comparison axis here
                idxCompare = idxCompareResample{iResample};

                [distanceByTime{iResample} meanDistanceResampled{iResample} timeBins] = ...
                    obj.getPSTHDistance(compareAcross, def, 'idxCompare', idxCompare, ...
                    'namesCompare', namesCompare);
            end

            % distanceByTime is a nResample x cell array of arrays of size nCompare x nTimeBins
            % distanceByTimeMat is a nCompare x nTimeBins x nResample array 
            distanceByTimeMat = cell2mat(shiftdim(distanceByTime,-2));

            % and we compute p/2 percentile edges of the data along the shuffle axis as the floor per comparison
            szMat = size(distanceByTimeMat);
            distanceByTimeIntervals = nan(szMat(1), szMat(2), 2);
            distanceByTimeIntervals(:, :, 1) = prctile(distanceByTimeMat, 100*(p/2), 3);
            distanceByTimeIntervals(:, :, 2) = prctile(distanceByTimeMat, 100*(1-p/2), 3);

            % meanDistanceMat is nCompare x nShuffles
            meanDistanceMat = cell2mat(shiftdim(meanDistanceResampled, -1));
            meanDistanceInterval = nan(size(meanDistanceMat, 1), 2);
            meanDistanceInterval(:,1) = prctile(meanDistanceMat, 100*(p/2), 2);
            meanDistanceInterval(:,2) = prctile(meanDistanceMat, 100*(1-p/2), 2);
        end

        function figh = drawPSTHDistance(obj, compareAcross, varargin)
            drawTimeAxis = true;
            axh = [];
            color = 0.4 *ones(1,3);
            significantColor = 'r';
            rectify = false;
            splitPlots = true;
            
            colorBySignificance = false; % if true, gray below contour, red above, else use appearanceFn colors
            thickenBySignificance = false; % if true, show thicker lines where significantly different from zero
            alphaNonSignificant = 0.5;
            showShuffleContour = false;
            showShuffled = false;

            showResampledContour = false;
            showResampled = false;

            showResampledFromSameContour = false;
            showResampledFromSame = false;

            shuffleColor = 0.8 * ones(1,3);
            resampleColor = [1 0.6 0.6];
            resampleFromSameColor = 0.6 * ones(1,3);
            def = assignargs(varargin);

            [distanceByTime meanDistance timeBins namesCompareCommon] = ...
                obj.getPSTHDistance(compareAcross, def);
            
            if colorBySignificance || showShuffleContour || showShuffled || thickenBySignificance
                [distanceByTimeFloor meanDistanceFloor distanceByTimeMat] = ...
                    obj.getPSTHDistanceShuffled(compareAcross, def);
            end
            
            if showResampledContour || showResampled
                [distanceByTimeIntervals meanDistanceIntervals distanceByTimeReshuffleMat meanDistanceReshuffleMat] = ...
                    obj.getPSTHDistanceResampled(compareAcross, def);
            end
            
            if showResampledFromSameContour || showResampledFromSame
                [distanceByTimeIntervalsFromSame meanDistanceIntervalsFromSame distanceByTimeReshuffleFromSameMat meanDistanceReshuffleFromSameMat] = ...
                    obj.getPSTHDistanceResampledFromSame(compareAcross, def);
            end

            [idxCompare namesCompare appearanceCompare] = ...
                obj.conditionInfo.getComparisonAxis(compareAcross);
            
            if isempty(axh)
                figh = fig(1050);
                axh = gca();
            end

            if colorBySignificance
                significantByTime = distanceByTime;
                significantByTime(abs(distanceByTime) < abs(distanceByTimeFloor)) = NaN;
            end

            nCompare = size(distanceByTime, 1);
            for iCompare = 1:nCompare
                if splitPlots
                    axh = subplot(nCompare, 1, iCompare);
                end

                if iCompare == 1
                    title(obj.name);
                end

                hold on
                if splitPlots && showShuffled
                    plot(timeBins, squeeze(distanceByTimeMat(iCompare, :, :))', 'Color', shuffleColor, 'Parent', axh);
                end

                if splitPlots && showShuffleContour
                    plot(timeBins, distanceByTimeFloor(iCompare, :), 'Color', shuffleColor, 'Parent', axh);
                    if ~rectify
                        plot(timeBins, -distanceByTimeFloor(iCompare, :), 'Color', shuffleColor, 'Parent', axh);
                    end
                end

                if showResampledContour
                    plot(timeBins, squeeze(distanceByTimeIntervals(iCompare, :, :)), ...
                        'Color', resampleColor, 'Parent', axh);
                end
                if showResampledFromSameContour
                    plot(timeBins, squeeze(distanceByTimeIntervalsFromSame(iCompare, :, :)), ...
                        'Color', resampleFromSameColor, 'Parent', axh);
                end

                if showResampled
                    plot(timeBins, squeeze(distanceByTimeReshuffleMat(iCompare, :, :)), ...
                        'Color', resampleColor, 'Parent', axh);
                end
                if showResampledFromSame
                    plot(timeBins, squeeze(distanceByTimeReshuffleFromSameMat(iCompare, :, :)), ...
                        'Color', resampleFromSameColor, 'Parent', axh);
                end

                if ~rectify
                    plot(timeBins, 0*timeBins, 'k-'); 
                end

                if colorBySignificance 
                    nonSigColor = color;
                    nonSigWidth = appearanceCompare{iCompare}(1).lineWidth;
                    sigColor = significantColor;
                    sigWidth = nonSigWidth;
                else
                    condColor = appearanceCompare{iCompare}(1).color;
                    if ischar(condColor)
                        nonSigColor = condColor;
                    else
                        nonSigColor = alphaNonSignificant*condColor + (1-alphaNonSignificant)*[1 1 1];
                    end
                    nonSigWidth = 2;
                    sigColor = condColor;
                    sigWidth = 3*nonSigWidth;
                end

                plot(timeBins, distanceByTime(iCompare, :), '-', 'Color', nonSigColor, 'LineWidth', nonSigWidth,'Parent', axh);
                
                if colorBySignificance || thickenBySignificance
                    plot(timeBins, significantByTime(iCompare, :), '-', 'LineWidth', sigWidth, 'Color', sigColor, 'Parent', axh);
                end

                hold off
                box off
                if splitPlots
                    ylabel(sprintf('%s', namesCompareCommon{iCompare}));
                else
                    ylabel('Delta (Hz)');
                end
                
                xlim([min(timeBins) max(timeBins)]);
                
                if drawTimeAxis
                    obj.drawTimeAxis('axh', axh);
                end
            end

        end
           
        function [fighStem fighHistogram] = drawPSTHDistanceMean(obj, compareAcross, varargin)
            p = 0.05;
            rectify = false;
            plotHistogram = true;
            def = assignargs(varargin);

            [distanceByTime meanDistance timeBins namesCompareCommon] = obj.getPSTHDistance(compareAcross, def); 
            [distanceByTimeFloor meanDistanceFloor distanceByTimeMat meanDistanceShuffles] = ...
                obj.getPSTHDistanceShuffled(compareAcross, def);

            nCompare = size(distanceByTime, 1);
            maxDist = max([max(abs(meanDistance)), max(abs(meanDistanceShuffles(:)))]);
            
            if plotHistogram
                fighHistogram = fig(156);
                figname('Distance shuffle control');

                for iCompare = 1:nCompare
                    axh = subplot(nCompare, 1, iCompare);
                    shuffleVals = meanDistanceShuffles(iCompare, :);
                    mdist = meanDistance(iCompare);
                    mfloor = meanDistanceFloor(iCompare);
                    if abs(mdist) > mfloor
                        arrowColor = 'r';
                    else
                        arrowColor = 'k';
                    end
                
                    if rectify
                        bins = linspace(0, maxDist, 20);
                    else
                        bins = linspace(-maxDist, maxDist, 20);
                    end

                    if rectify
                        lineAt = mfloor;
                    else
                        lineAt = [-mfloor mfloor];
                    end

                    plotHist(shuffleVals, bins, 'axh', axh, ... 
                        'arrowAt', mdist, 'arrowLabel', sprintf('%.3g Hz', mdist), 'arrowColor', arrowColor, ...
                        'lineAt', lineAt, 'lineColor', 0.6*ones(3,1));
                     

                    ylabel(namesCompareCommon{iCompare});
                    if iCompare == 1 
                        title(obj.name);
                    end
                    if iCompare == nCompare
                        xlabel('Delta (Hz)');
                    end

                   % makePrettyAxis;
                end
                figsize(8, 6);
            end

            fighStem = fig(157);
            figname('Delta PSTH');
            xl = [0 nCompare+1];
            plot(xl, [0 0], '--', 'Color', 0.6*ones(3,1));
            hold on
            for iCompare = 1:nCompare
                mdist = meanDistance(iCompare);
                mfloor = meanDistanceFloor(iCompare);

                if abs(mdist) > mfloor 
                    color = 'r';
                else
                    color = 'k';
                end

                plot(iCompare, mdist, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'MarkerSize', 6);
                hold on
                plot([iCompare iCompare], [mdist mdist-sign(mdist)*mfloor], '-', 'Color', color);
            end

            set(gca, 'XLim', xl, 'XTick', 1:nCompare, 'XTickLabel', namesCompareCommon);
            ylabel('Delta (Hz)');
            title(obj.name);
            makePrettyAxis;

        end
    end

    methods % quantification of firing rates
        function [meanHz semHz] = getMeanHz(obj)
            meanHzByTrials = obj.getMeanHzByTrial();
            meanHz = nanmean(meanHzByTrials);
            semHz = nansem(meanHzByTrials);
        end


        function [meanHzByTrial] = getMeanHzByTrial(obj)
            meanHzByTrial = nan(obj.nTrials, 1);

            nSpikesPerTrial = obj.nSpikesPerTrial;
            tLengthSecPerTrial = obj.tLengthByTrial ./ obj.tUnitsPerSec;
            meanHzByTrial = nSpikesPerTrial ./ tLengthSecPerTrial;
        end

        function [meanHzByCondition semHzByCondition] = getMeanHzByCondition(obj)
            idxByCond = obj.conditionInfo.getGroups();
            nConds = length(idxByCond);
            meanHzByCondition = nan(nConds,1);
            semHzByCondition = nan(nConds,1);

            meanHzByTrial = obj.getMeanHzByTrial();
            
            for iCond = 1:nConds
                idx = idxByCond{iCond};              
                nTrials = length(idx);

                meanHzByCondition(iCond) = nanmean(meanHzByTrial(idx));
                semHzByCondition(iCond) =  nansem(meanHzByTrial(idx));
            end
        end
        
        function [meanHzByCondition semHzByCondition] = plotMeanHzByCondition(obj, varargin)
            [meanHzByCondition semHzByCondition] = obj.getMeanHzByCondition();
            cla;
            
            mask = obj.conditionInfo.countByCondition > 0;
            
            xval = 1:nnz(mask);
            errorbar(xval, meanHzByCondition(mask), semHzByCondition(mask), semHzByCondition(mask), ...
                'k-');
            
            set(gca, 'XTick', xval, 'XTickLabel', obj.conditionInfo.names(mask));
            %makePrettyAxis;           
            xlabel('Condition');
            ylabel('Mean Firing Rate (sp/s)')
            box off
        end
        

        function [maxHzByCondition minHzByCondition] = getExtremaHzByCondition(obj)
            psthByCond = obj.getPSTHByCondition();
            maxHzByCondition = max(psthByCond, [], 2);
            minHzByCondition = min(psthByCond, [], 2);
        end
        
        function [maxHz minHz] = getExtremaHz(obj)
            % get extrema of each PSTH
            [maxHzByCondition minHzByCondition] = getExtremaHzByCondition(obj);
            
            maxHz = max(maxHzByCondition);
            minHz = min(minHzByCondition);
        end
            
        function [tuningDepth time] = getTuningDepthOverTime(obj)
            [psthByCond time] = obj.getPSTHByCondition();
            tuningDepth = max(psthByCond, [], 1) - min(psthByCond, [], 1);
        end

        function [maxTuningDepth] = getMaxTuningDepth(obj)
            maxTuningDepth = max(obj.getTuningDepthOverTime());
        end

        function [meanTuningDepth, semTuningDepth] = getMeanTuningDepth(obj)
            tuningDepth = obj.getTuningDepthOverTime();
            meanTuningDepth = nanmean(tuningDepth);
            semTuningDepth = nansem(tuningDepth);
        end

        function [stdByCond] = getTemporalVariationByCondition(obj)
            % compute the std deviation
            [psthByCond time] = obj.getPSTHByCondition();
            stdByCond = nanstd(psthByCond, [], 2);
        end

        function maxVariation = getMaxTemporalVariation(obj);
            maxVariation = max(obj.getTemporalVariationByCondition()); 
        end
    end

    methods % plotting Raster and PSTHs
        function drawTimeAxis(obj, varargin)
            p = inputParser();
            p.addParamValue('trialIdx', true(obj.nTrials, 1), @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            ti = obj.alignInfo.timeInfo(p.Results.trialIdx);
            
           % obj.alignInfo.drawTimeAxis(ti, 'setXLim', true);
        end

        function hLine = drawSpikesForCondition(obj, iCond, varargin)
            axh = [];
            tickHeight = 0.99;
            rowHeight = 1;
            xOffset = 0;
            yOffset = 0; % offset for top of spikes
            tLims = [];
            assignargs(varargin); 

            if obj.conditionInfo.countByCondition(iCond) == 0
                return;
            end

            if isempty(axh)
                axh = gca;
                cla(axh);
            end

            % build line commands
            spikes = obj.getSpikesForCondition(iCond);
            nTrials = length(spikes);

            XByTrial = cell(1, nTrials);
            YByTrial = cell(1, nTrials);
            for iE = 1:nTrials 
                XByTrial{iE} = repmat(makerow(spikes{iE}), 2, 1);
                YByTrial{iE} = repmat([-rowHeight*(iE-1); -rowHeight*(iE-1)-tickHeight], 1, length(spikes{iE}));
            end
            
            X = cell2mat(XByTrial) + xOffset;
            Y = cell2mat(YByTrial) + yOffset;
            
            % filter within time limits?
            if ~isempty(X)
                if ~isempty(tLims)
                    withinTLimsMask = X(1,:) >= tLims(1) & X(1,:) <= tLims(2);
                    X = X(:, withinTLimsMask);
                    Y = Y(:, withinTLimsMask);
                end

                %hLine = line(X, Y, 'Parent', axh, 'Color', 'k', 'LineWidth', conditionAppearance(iCond).lineWidth);
                hLine = line(X, Y, 'Parent', axh, 'Color', 'k', 'LineWidth', 0.2);
            else
                hLine = NaN;
            end
        end

        function hLine = drawSpikes(obj, varargin)
            rowHeight = 1;
            yOffset = 0;
            yConditionGap = 0.5; % as percent multiplier of total trial count, excluding interval gaps
            drawPrettyAxis = true;
            drawConditionNames = true;
            drawTimeAxis = true;
            drawTitle = true;
            drawIntervals = false;
            intervalStartGap = 0.25; % as percent multiplier of total trial count
            intervalHeight = 1.0; % as percent multiplier of total trial count
            drawTimeScaleBar = false;
            timeScaleBarWidth = 100;
            timeScaleBarHeight = 1.5; % in % of trials
            axh = [];
            assignargs(varargin);
            
            if isempty(axh)
                clf;
                axh = gca;
            end
            
            % convert to actual units
            yConditionGap = yConditionGap * obj.nTrials / 100;
            intervalHeight = intervalHeight * obj.nTrials / 100;
            intervalStartGap = intervalStartGap * obj.nTrials / 100;
            
            if drawIntervals
                yConditionGap = yConditionGap + intervalHeight + intervalStartGap;
                infoByCond = obj.alignInfo.getIntervalInfoByCondition(obj.conditionInfo);
            end
            
            %tLims = obj.alignInfo.getTimeAxisLims();

            % work from top to bottom, figure out starting position
            yStartByCondition = nan(obj.nConditions,1);
            yEndByCondition = nan(obj.nConditions, 1);
            yStart = yOffset + obj.nTrials*rowHeight + (obj.conditionInfo.nConditionsNonEmpty-1)*yConditionGap;
            yEnd = yStart;
            for iCond = 1:obj.nConditions;
                if obj.nTrialsByCondition(iCond) > 0
                    yStartByCondition(iCond) = yEnd;
                    yEndByCondition(iCond) = yEnd - obj.nTrialsByCondition(iCond)*rowHeight;
                    yEnd = yEndByCondition(iCond) - yConditionGap;
                end
            end
            yEnd = yEnd + yConditionGap;
            
            if drawIntervals
                yStart = yStart + intervalHeight + intervalStartGap;
            end
            
            for iCond = 1:obj.nConditions;
                if obj.nTrialsByCondition(iCond) > 0
                    obj.drawSpikesForCondition(iCond, 'yOffset', yStartByCondition(iCond), ...
                        'axh', axh, 'rowHeight', rowHeight);
                    hold on
                    if drawIntervals
                        % place above the first row of spikes
                        yIntStart = yStartByCondition(iCond) + intervalStartGap;
                        yIntStop = yStartByCondition(iCond) + intervalStartGap + intervalHeight;
                        
                        for iInt = 1:length(infoByCond(iCond).interval)
                            intervals = infoByCond(iCond).interval{iInt};
                            % intervals is nPeriods x 2
                            % poly* is 4 x nPeriods, top-left, bot-left, bot-right, bot-left
                            polyX = [intervals(:,1)'; intervals(:,1)'; intervals(:,2)'; intervals(:,2)'];
                            polyY = repmat([yIntStart; yIntStop; yIntStop; yIntStart], 1, size(polyX, 2));
                            h = patch(polyX, polyY, obj.alignInfo.intervalColors{iInt});
                            set(h, 'EdgeColor', 'none');
                            hold on
                        end
                    end
                end
            end

            hold off
            
            if drawTitle
                title(obj.name);
            end
            
            %xlim(tLims);
            ylim([yEnd yStart]);
            
            if drawTimeScaleBar
                timeScaleBarHeight = timeScaleBarHeight*(yStart/ 100);
                yh = yEnd-yConditionGap - timeScaleBarHeight;
                yl = yEnd-yConditionGap - 2*timeScaleBarHeight;
                xlims = get(gca, 'XLim');
                xh = xlims(2);
                xl = xlims(2) - timeScaleBarWidth;
                h = patch([xl;xl;xh;xh], [yl;yh;yh;yl], 'k');
                set(h, 'EdgeColor', 'none');
                
                str = sprintf('%d ms', timeScaleBarWidth);
                text(xh, yh, str, 'VerticalAlign', 'bottom', 'HorizontalAlign', 'right');
                ylim([yl yStart]);
            end

            if drawPrettyAxis
                if drawTimeAxis
                    obj.drawTimeAxis('axh', axh, 'drawY', false);
                else
                    %set(axh, 'XTick', [], 'XTickLabel', {});
                end
                
                for iCond = 1:obj.nConditions
                    if obj.nTrialsByCondition(iCond) == 0
                        continue;
                    end
                    if drawConditionNames
                        condName = obj.conditionNames{iCond};
                    else
                        condName = '';
                    end
%                     drawAxis([yStartByCondition(iCond) mean([yStartByCondition(iCond) yEndByCondition(iCond)]) yEndByCondition(iCond)], ...
%                         'axisOrientation', 'v', ...
%                         'tickLabels', {'', condName, ''}, ...
%                         'color', obj.conditionAppearance(iCond).color, ...
%                         'labelColor', obj.conditionAppearance(iCond).color);
                end
            end
        end

        function figh = drawPSTH(obj, varargin)
            xOffset = 0;
            yOffset = 0;
            showLegend = true;
            legendPosition = 'Best';
            axh = [];
            drawTimeAxis = true;
            idxByCondition = [];
            names = [];
            appearance = [];
            drawLegend = true;
            yLabel = 'spikes / sec';
            drawSem = false;
            drawTitle = true;
            yLim = [];
            assignargs(varargin);

            if obj.nTrialsValid == 0
                warning('No valid trials found for this SpikeRaster');
            end
            
            if isempty(axh)
                figh = gcf;
                clf
                axh = gca;
            else
                figh = get(axh, 'Parent');
            end

            if isempty(appearance) || isempty(names) || isempty(idxByCondition) 
                idxByCondition = obj.conditionInfo.listByCondition;
                appearance = obj.conditionAppearance;
                names = obj.conditionNames;
            end
            nConditions = numel(idxByCondition);

            [psthByCond semByCond time] = obj.getPSTHByCondition('idxByCondition', idxByCondition);
  
            if all(isnan(psthByCond(:)))
                warning('Not enough trials per condition to build any PSTHs');
            end
            
            h = nan(nConditions,1);
            
            for iCond = 1:nConditions
                if isempty(idxByCondition{iCond})
                    continue;
                end
                lineWidth = appearance(iCond).lineWidth;
                lineWidthMean = lineWidth;

                % plot central mean with flanking standard error if asked
                h(iCond) = plot(time, psthByCond(iCond, :), 'Parent', axh, ...
                    'Color', appearance(iCond).color, 'LineWidth', lineWidthMean); 
                hold on
                if drawSem
                    lineWidthSem = 1; %lineWidth/2;
                    semTraces = [psthByCond(iCond,:)+semByCond(iCond,:); psthByCond(iCond,:)-semByCond(iCond,:)]';
                    semTraces(semTraces < 0) = NaN;
                    plot(time, semTraces, ...
                       'Color', appearance(iCond).color, 'LineWidth', lineWidthSem);
                end
            end
            
            if ~isnan(obj.tMin) && ~isnan(obj.tMax)
                xlim([obj.tMin obj.tMax]);
            end
            hold off
            xlabel(sprintf('time (%s)', obj.tUnits));
            if ~isempty(yLabel)
                ylabel(yLabel);
            end
  
            if drawLegend && nConditions > 1
                legend(h(~isnan(h)), names(~isnan(h)), 'Location', legendPosition);
                legend('boxoff');
            end

            box(axh, 'off');
            if drawTitle
                title(obj.name);
            end

            if ~isempty(yLim)
                ylim(yLim);
            end

            if drawTimeAxis
                allIdx = cell2mat(makecol(idxByCondition(:)));
                obj.drawTimeAxis('trialIdx', allIdx);
            end

        end

        function figh = drawSpikesAndPSTH(obj, varargin)
            rowHeight = 1;
            yConditionGap = 10;
            axh = [];
            drawConditionNames = true;
            yLabelPSTH = 'spikes/sec';
            psthYLim = [];
            assignargs(varargin);

            if isempty(axh)
                clf
                figh = figure();
                axh = gca();
            else
                figh = get(axh, 'Parent');
            end

            obj.drawPSTH(varargin{:}, 'axh', axh, 'yLabel', yLabelPSTH, 'drawLegend', false, 'drawTitle', false, 'yLim', psthYLim);
            ylims = ylim();
            yOffset = ylims(2) + yConditionGap;
            obj.drawSpikes(varargin{:}, 'axh', axh, 'drawConditionNames', drawConditionNames, 'drawPrettyAxis', true, 'drawTimeAxis', false, 'yOffset', yOffset);

            yFinal = yOffset + sum(obj.nTrialsByCondition(:))*rowHeight + yConditionGap*(obj.conditionInfo.nConditionsNonEmpty-1);
            ylim([0 yFinal]);
        end

        function drawPSTHCompare(obj, attrNames, varargin)
            if nargin < 2
                attrNames = {};
            end
             
            [idxCompare, namesCompare, appearanceCompare, ...
                conditionIdxCell, conditionDescriptorOuter, ...
                conditionDescriptorInnerCell] = obj.conditionInfo.getComparisonAxis(attrNames);
            
            % filter everything by comparisons with at least one trial inside
            nonEmptyMask = cellfun(@(idxCell) any(~cellfun(@isempty, idxCell(:))), idxCompare);
            conditionDescriptorInnerCell = conditionDescriptorInnerCell(nonEmptyMask);
            idxCompare = idxCompare(nonEmptyMask);
            namesCompare = namesCompare(nonEmptyMask);
            nCompare = numel(idxCompare);

            clf;
            p = panel();
            p.de.margin = 0;
            p.pack(nCompare,1);
            
            for iCompare = 1:nCompare
                idx = idxCompare{iCompare};
                namesInner = conditionDescriptorInnerCell{iCompare}.names;
                appearanceInner = conditionDescriptorInnerCell{iCompare}.appearances;

                p(iCompare, 1).select();
                obj.drawPSTH(varargin{:}, 'axh', gca, 'drawTitle', iCompare==1, ...
                    'idxByCondition', idx, 'names', namesInner, 'appearance', appearanceInner, ...
                    'yLabel', namesCompare{iCompare}); 
                axis tight;
            end
        end
        
        function drawPSTHCompareWithin(obj, attrNames, varargin)
            otherAttr = setdiff(obj.conditionInfo.groupByList, attrNames);
            obj.drawPSTHCompare(otherAttr, varargin{:});
        end
    end

    methods % resample and shuffle controls
        function srCopy = buildShuffledAlong(sr, compareAlong)
            srCopy = sr.copy();
            [srCopy.conditionInfo newIdxList] = sr.conditionInfo.buildShuffledAlong(compareAlong);
            srCopy.resampleTrials(newIdxList, 'resampleConditionInfo', false);
        end

        function srCopy = buildResampled(sr)
            srCopy = sr.copy();
            [srCopy.conditionInfo newIdxList] = sr.conditionInfo.buildResampled();
            srCopy.resampleTrials(newIdxList, 'resampleConditionInfo', false);
        end

        function srCopy = buildResampledFromSingleAttributeValue(sr, attr, value)
            srCopy = sr.copy();
            [srCopy.conditionInfo newIdxList] = sr.conditionInfo.buildResampledFromSingleAttributeValue(attr, value);
            srCopy.resampleTrials(newIdxList, 'resampleConditionInfo', false);
        end

        % legacy functions
        function srShuffle = buildShuffledAlongComparisonAxis(sr, compareAcross, varargin)
            p = inputParser();
            p.addRequired('compareAcross', @ischar);
            p.addParamValue('nShuffles', 1, @isscalar);
            p.parse(compareAcross, varargin{:});
            nShuffles = p.Results.nShuffles;

            ciCell = sr.conditionInfo.buildShuffledAlongComparisonAxis(compareAcross, ...
                'nShuffles', nShuffles);

            srShuffle = cell(nShuffles, 1);
            for i = 1:nShuffles
                srShuffle{i} = sr.copy();
                srShuffle{i}.filterValidTrials();
                srShuffle{i}.conditionInfo = ciCell{i};
                srShuffle{i}.clearCaches();
            end

            if nShuffles == 1
                srShuffle = srShuffle{1};
            end
        end
        
        function srResample = buildResampledWithinCondition(sr, varargin)
            p = inputParser();
            p.addParamValue('nResample', 1, @isscalar);
            p.parse(varargin{:});
            nResample= p.Results.nResample;

            srResample = cell(nResample, 1);
            idxCell = sr.conditionInfo.buildResampledWithinCondition('nResample', nResample);
            for i = 1:nResample
                srResample{i} = sr.copy();
                srResample{i}.resampleTrials(idxCell{i});
            end
            
            if nResample == 1
                srResample = srResample{1};
            end
        end
    end
end
