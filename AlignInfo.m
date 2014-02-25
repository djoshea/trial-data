classdef AlignInfo < AlignDescriptor
% AlignInfo is a subclass of AlignDescriptor that may be bound to a set of trials
% where the events have actual timestamps by trial.

    properties(SetAccess=protected)
        % AlignInfoOnDemandCache instance
        odc
    end

    properties(SetAccess=protected)
        % this function maps (R, eventList) --> eventTimes array nTrials x nEvents
        getEventTimesFn = @AlignInfo.defaultGetEventTimesFn;

        % struct with fieldnames == events.
        % each value is nTrials x 1 numeric or cell array 
        % containing the absolute times of that event
        % for the trials, as returned by getEventTimesFn
        eventData
        
        % struct with fieldnames == events
        % each is nTrials x 1 array of event counts
        eventCounts

        manualInvalid

        applied = false; 
    end
    
    % properties which are computed on demand and stored in odc
    properties(Transient, Dependent, SetAccess=protected)
        % struct array of nTrials x 1 containing the absolute times of the
        % start, stop, zero events, as well as intervals and marks 
        timeInfo
        
        % valid is a dependent property formed by merging computedValid
        % with manualInvalid
        computedValid
        
        % nTrials x nMarks cell array of mark times
        markData
        
        % nTrials x nMarks array of mark occurrence counts
        markCounts
        
        % nIntervals cell array of nTrials x nMaxOccurrences interval start times
        intervalStartData
        
        % nIntervals cell array of nTrials x nMaxOccurrences interval stop times
        intervalStopData
        
        % nTrials x nIntervals array of interval occurrence counts
        intervalCounts
    end
    
    properties(Dependent)
        valid
        nTrials
    end
    
    methods
        function v = get.timeInfo(ad)
            v = ad.odc.timeInfo;            
            if isempty(v)
                ad.buildTimeInfo();
                v = ad.odc.timeInfo;
            end
        end
        
        function ad = set.timeInfo(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.timeInfo = v;
        end
        
        function v = get.computedValid(ad)
            v = ad.odc.computedValid;            
            if isempty(v)
                ad.buildTimeInfo();
                v = ad.odc.computedValid;
            end
        end
        
        function ad = set.computedValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.computedValid = v;
        end
        
        function v = get.markData(ad)
            v = ad.odc.markData;            
            if isempty(v)
                ad.buildMarkData();
                v = ad.odc.markData;
            end
        end
        
        function ad = set.markData(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.markData = v;
        end
        
        function v = get.markCounts(ad)
            v = ad.odc.markCounts;            
            if isempty(v)
                ad.buildMarkData();
                v = ad.odc.markCounts;
            end
        end
        
        function ad = set.markCounts(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.markCounts = v;
        end

        function v = get.intervalStartData(ad)
            v = ad.odc.intervalStartData;            
            if isempty(v)
                ad.buildIntervalData();
                v = ad.odc.intervalStartData;
            end
        end
        
        function ad = set.intervalStartData(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalStartData = v;
        end
        
        function v = get.intervalStopData(ad)
            v = ad.odc.intervalStopData;            
            if isempty(v)
                ad.buildIntervalData();
                v = ad.odc.intervalStopData;
            end
        end
        
        function ad = set.intervalStopData(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalStopData = v;
        end

        function v = get.intervalCounts(ad)
            v = ad.odc.intervalCounts;            
            if isempty(v)
                ad.buildIntervalData();
                v = ad.odc.intervalCounts;
            end
        end
        
        function ad = set.intervalCounts(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalCounts = v;
        end
    end
    
    methods(Static) % construct from another align descriptor, used primarily by AlignInfo
        function adNew = fromAlignDescriptor(ad, adNew)
            if nargin < 2
                % allow subclasses to provide their own instance to
                % populate
                adNew = AlignInfo();
            end
            adNew = AlignDescriptor.fromAlignDescriptor(ad, adNew);
            adNew.applied = false;
        end
    end
    
    methods % build a AlignDescriptor back from this AlignInfo
        % build a static ConditionDescriptor for the current groupByList
        function ad = getAlignDescriptor(ad)
            ad = AlignDescriptor.fromAlignDescriptor(ad);
        end
    end

    methods % Construct and bind to trial data or event struct
        function ad = AlignInfo(varargin)
            ad = ad@AlignDescriptor(varargin{:});
            ad.odc = AlignInfoOnDemandCache();
        end

        % bind this AlignInfo to a set of trials
        function ad = applyToTrialData(ad, R)
            ad.warnIfNoArgOut(nargout);
            [ad.eventData, ad.eventCounts] = ad.requestEventInfo(R);
            ad.applied = true;
            ad.manualInvalid = falsevec(ad.nTrials);
            ad = ad.update();
        end

        function printOneLineDescription(ad)
            desc = ad.getStartStopZeroPadDescription();
            if ad.applied
                validStr = sprintf('(%d valid)', nnz(ad.computedValid));
            else
                validStr = '(not applied)';
            end
            if ~ad.nameDefault
                tcprintf('inline', '{yellow}%s: {bright blue}%s : {none}%s %s\n', class(ad), ad.name, desc, validStr);
            else
                % name will just match desc, so don't print it twice
                tcprintf('inline', '{yellow}%s: {none}%s %s\n', class(ad), desc, validStr);
            end
        end

        function ad = update(ad)
            ad.warnIfNoArgOut(nargout);
            
            if ~ad.applied
                % nothing to udpate if we haven't applied to trial data yet
                return;
            end
            
            ad.odc = ad.odc.copy();
            ad.odc.flush();
        end
        
        function ad = updateMark(ad)
            ad.warnIfNoArgOut(nargout);
            
            if ~ad.applied
                % nothing to udpate if we haven't applied to trial data yet
                return;
            end
            
            ad.odc = ad.odc.copy();
            ad.odc.flushMarkData();
        end
        
        function ad = updateInterval(ad)
            ad.warnIfNoArgOut(nargout);
            if ~ad.applied
                % nothing to udpate if we haven't applied to trial data yet
                return;
            end
            
            ad.odc = ad.odc.copy();
            ad.odc.flushIntervalData();
        end
        
        % internal use function that simply grabs the event times relative to the trial
        % start and start/stop times as well
        % eventInfo(iTrial).event_eventName is time of eventName in trial
        % iTrial relative to eventInfo(iTrial).event_start
        function [eventData, eventCounts] = requestEventInfo(ad, td)
            eventList = ad.getEventList();
            
            % grab the event times for all needed events
            [eventData, eventCounts] = ad.getEventTimesFn(td, eventList);
        end
        
        function assertHasEvent(ad, eventName)
            if ad.applied
                assert(isfield(ad.eventData, eventName), ...
                    'AlignInfo has no event named %s', eventName);
            end
        end

        % internal utility functions for accessing specific event times
        
        % n must be a scalar, times is a numeric array
        function times = getEventNthTimeVector(ad, event, n, offset, roundRelativeTo)
            ad.assertHasEvent(event);
           
            timesMat = ad.eventData.(event);
            counts = ad.eventCounts.(event);
            
            if strcmp(n, 'end')
                % TODO implement 'end-1' type index
                inds = sub2ind(size(timesMat), 1:ad.nTrials, counts);
                times = timesMat(inds);
            else
                times = timesMat(:, n);
            end
           
            times = times + offset;
            
            if ~isempty(roundRelativeTo) && ad.roundTimes
                % shift event times to maintain an integer multiple of
                % timeDelta offset from .zero
                times = ad.roundTimesRelativeTo(times, roundRelativeTo);
            end
        end
        
        % n may be an index or a selector (e.g. '2:end')
        function [timeMatrix, counts] = getEventIndexedTimeMatrix(ad, event, n, offset, roundRelativeTo)
            % similar to getEventIndexedTimeVector but returns a nTrials x nMaxOccurrences matrix,
            % with a variable number of event occurrences, and n may be be a
            % string of the form '1:2', '1:end', ':', etc
            %
            % trials with fewer occurrences will be padded with NaNs to
            % complete the matrix form
            
            ad.assertHasEvent(event);
            
            if ~ischar(n) || strcmp(n, 'end')
                % defer to simpler function when scalar or 'end'
                timeMatrix = ad.getEventNthTimeVector(event, n, offset, roundRelativeTo);
                counts = double(~isnan(timeMatrix));
                
            else
                % must have a colon, parse into tokens, end keyword okay
                % before and/or after colon
                pat = '(?<end1>end)?(?<ind1>-?\d*)?:(?<end2>end)?(?<ind2>-?\d*)?';
                info = regexp(n, pat, 'names', 'once');

                if isempty(info)
                    error('Unable to parse event index %s(%s)', event, n);
                end
                
                % info.ind1 will contain the scalar offset, 
                % info.end1 will be 'end' if end is part of the expression 
                % same for .ind2, .end2

                % convert ind1, ind2 to doubles
                end1 = ~isempty(info.end1);
                if isempty(info.ind1)
                    if end1             
                        ind1 = 0; % default is end+0
                    else
                        ind1 = 1; 
                    end
                else
                    ind1 = str2double(info.ind1);
                end
                
                if isempty(info.ind2)
                    % default to :end, even if not specified as 'end'
                    end2 = true;
                    ind2 = 0;
                else
                    end2 = ~isempty(info.end2);
                    ind2 = str2double(info.ind2);
                end
                
                countsFull = ad.eventCounts.(event);
                nTrials = numel(countsFull);             
                timesFull = ad.eventData.(event);
   
                
                % figure out start and end columns to index into nTrials x
                % nMaxOccurrence matrix timesFull
                if end1
                    colStart = countsFull + ind1;
                else
                    colStart = repmat(ind1, nTrials, 1);
                end
                if end2
                    colStop = countsFull + ind2;
                else
                    colStop = repmat(ind2, nTrials, 1);
                end
                    
                % compute the counts
                counts = colStop - colStart + 1;
                
                % and populate the matrix row by row, leaving NaNs where
                % fewer events occurred
                timeMatrix = nan(nTrials, max(counts));
                for i = 1:nTrials
                    timeMatrix(i, 1:counts(i)) = timesFull(i, colStart(i):colStop(i));
                end
            end
            
            timeMatrix = timeMatrix + offset;
            
            if ~isempty(roundRelativeTo)
                % shift event times to maintain an integer multiple of
                % timeDelta offset from .zero
                timeMatrix = ad.roundTimesRelativeTo(timeMatrix, roundRelativeTo);
            end
        end
        
        % same as above but empty cells are filled with NaN
%         function timeCell = getEventIndexedTimeCellFillEmptyWithNaN(ad, event, n, offset, roundRelativeTo)
%             timeCell  = ad.getEventIndexedTimeCell(event, n, offset, roundRelativeTo);
%             emptyMask = cellfun(@isempty, timeCell);
%             [timeCell{emptyMask}] = deal(NaN);
%         end
        
        function times = roundTimesRelativeTo(ad, times, ref)
            % rounds time points such that they align with a time reference
            % i.e. each time lies an exact integer multiple of minTimeDelta
            % away from ref
            
            if ~ad.roundTimes
                % no rounding enabled
                return;
            end
            
            delta = ad.minTimeDelta;
            roundFn = @(times, ref) round((times - ref) / delta) * delta + ref;
            if isnumeric(times)
                times = roundFn(times, ref);
            else
                times = cellfun(roundFn, times, num2cell(ref), 'UniformOutput', false);
            end
        end
        
        function times = filterTimesWithin(ad, times, startByTrial, stopByTrial) %#ok<INUSL>
            % filters times (numeric or cell) with size along dim1 == nTrials
            % to include only times within the start to stop window 
            if isnumeric(times)
                % numeric --> replace invalid with NaN
                times(times < startByTrial | times > stopByTrial) = NaN;
            else
                % cell matrix --> remove invalid times
                filterFn = @(vec, start, stop) vec(vec >= start & vec <= stop);
                times = map(filterFn, times, startByTrial, stopByTrial);
            end
        end
        
        function times = filterTimesWithinStartStop(ad, times, includePadding)
            if includePadding
                [startData, stopData] = ad.getStartStopZeroByTrial();
            else
                [startData, stopData] = ad.getStartStopZeroByTrialWithPadding();
            end
            
            times = ad.filterTimesWithin(times, startData, stopData);
        end
        
        % get the aligned start/stop/zero/mark time windows for each trial, 
        % respecting all truncation and invalidation instructions
        function buildTimeInfo(ad)
            % returns a struct array with the actual time window and time of zero for trial i as
            %   timeInfo(i).start, .stop, .zero
            %
            % timeInfo(i).valid and valid(i) indicate whether trial i satisfied inclusion criteria
            % as specified by the various means of trial invalidation 
           
            if ~ad.applied
                return;
            end
            
            padPre = ad.padPre;
            padPost = ad.padPost;

            % temporary, the dependent property .nTrials only works after this
            % is called
            nTrials = ad.nTrials;

            valid = truevec(nTrials);
            t.invalidCause = cellvec(nTrials);

            % get zero alignment event without rounding
            t.zero = ad.getEventNthTimeVector(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset, []);
            noZero = isnan(t.zero);
            [t.invalidCause{noZero & valid}] = deal(sprintf('Missing zero event %s', ad.zeroUnabbreviatedLabel));
            valid(noZero) = false;
            
            t.trialStart = ad.getEventNthTimeVector('TrialStart', 1, 0, t.zero); 
            t.trialStop = ad.getEventNthTimeVector('TrialEnd', 1, 0, t.zero);

            % get start event
            t.start = ad.getEventNthTimeVector(ad.startEvent, ad.startEventIndex, ad.startOffset, t.zero);
            noStart = isnan(t.start);
            [t.invalidCause{noStart}] = deal(sprintf('Missing start event %s', ad.startUnabbreviatedLabel));
            valid(noStart & valid) = false;

            % get stop event
            t.stop = ad.getEventNthTimeVector(ad.stopEvent, ad.stopEventIndex, ad.stopOffset, t.zero);
            noStop = isnan(t.stop);
            [t.invalidCause{noStop & valid}] = deal(sprintf('Missing stop event %s', ad.stopUnabbreviatedLabel));
            valid(noStop) = false;
        
            % get pad window
            t.startPad = t.start - padPre;
            t.stopPad = t.stop + padPost;
            
            % truncate trial end (including padding) based on truncateAfter events
            t.isTruncatedStop = falsevec(nTrials);
            for i = 1:length(ad.truncateAfterEvents)
                times = ad.getEventNthTimeVector(ad.truncateAfterEvents{i}, ...
                    ad.truncateAfterEventsIndex{i}, ad.truncateAfterOffsets(i), t.zero);
                t.isTruncatedStop = t.isTruncatedStop | times < t.stopPad;
                t.stopPad = nanmin(t.stop, times);
                t.stop = t.stopPad - padPost;
            end
            
            % truncate trial start (including padding) based on truncateBefore events
            t.isTruncatedStart = falsevec(nTrials);
            for i = 1:length(ad.truncateBeforeEvents)
                times = ad.getEventNthTimeVector(ad.truncateBeforeEvents{i}, ...
                    ad.truncateBeforeEventsIndex{i}, ad.truncateBeforeOffsets(i), t.zero);
                t.isTruncatedStart = t.isTruncatedStart | times > t.startPad;
                t.startPad = nanmax(t.startPad, times);
                t.start = t.startPad + padPre;
            end

            % mark trials as invalid if startPad:stopPad includes any invalidateEvents
            for i = 1:length(ad.invalidateEvents)
                timesMatrix = ad.getEventIndexedTimeMatrix(ad.invalidateEvents{i}, ...
                    ad.invalidateEventsIndex{i}, ad.invalidateOffsets(i), t.zero);
                
                maskInvalid = any(bsxfun(@ge, timesMatrix, t.startPad) & ...
                                  bsxfun(@le, timesMatrix, t.stopPad), 2); 
                
                valid(maskInvalid) = false;
                [t.invalidCause{maskInvalid}] = deal(sprintf('Overlapped invalidating event %s', ad.invalidateUnabbreviatedLabels{i}));
            end

            % handle windows which extend outside of trial
            if strcmp(ad.outsideOfTrialMode, ad.TRUNCATE) || ...
               strcmp(ad.outsideOfTrialMode, ad.INVALIDATE)

                if strcmp(ad.outsideOfTrialMode, ad.TRUNCATE) 
                    % truncate so that padded window fits within the trial
                    t.startPad(t.startPad < t.trialStart) = t.trialStart(t.startPad < t.trialStart);
                    t.stopPad(t.stopPad > t.trialStop) = t.trialStop(t.stopPad > t.trialStop);

                    t.start = t.startPad + padPre;
                    t.stop = t.stopPad - padPost;
                else
                    mask = t.startPad < t.trialStart;
                    valid(mask) = false;
                    [t.invalidCause{mask}] = deal('Start plus padding occurs before trial start');
                    
                    mask = t.stopPad > t.trialStop;
                    valid(mask) = false;
                    [t.invalidCause{mask}] = deal('Stop plus padding occurs after trial end');
                end
            end

            % handle minimum duration window
            mask = t.stop - t.start < ad.minDuration;
            valid(mask) = false;
            [t.invalidCause{mask}] = deal(sprintf('Trial duration is less than minDuration %g', ad.minDuration));
        
            % clear out values for invalid trials to avoid hard to catch bugs
            % IT IS CRITICAL THAT THIS ONLY CONSIDER computedValid, not
            % manual valid, otherwise changing the validity requires
            % timeInfo to be built again
            t.startPad(~valid) = NaN;
            t.stopPad(~valid) = NaN;
            t.start(~valid) = NaN;
            t.stop(~valid) = NaN;
            t.zero(~valid) = NaN;
            
            % store in odc
            c = ad.odc;
            c.timeInfo = t;
            c.computedValid = valid;
        end
        
        function buildMarkData(ad)
            % compute the mark times and store in odc.markInfo
            
            markData = cell(ad.nMarks, 1);
            markCounts = nan(ad.nTrials, ad.nMarks);
            for i = 1:length(ad.markEvents)
                [markData{i}, markCounts(:, i)] = ...
                    ad.getEventIndexedTimeMatrix(ad.markEvents{i}, ...
                    ad.markEventsIndex{i}, ad.markOffsets(i), ad.getZeroByTrial());
                
                markData{i}(~ad.computedValid, :) = NaN;
            end

            markCounts(~ad.computedValid, :) = NaN;
            
            c = ad.odc;
            c.markData = markData;
            c.markCounts = markCounts;
        end
        
        function buildIntervalData(ad)
            % compute the interval times by trial and store in
            % odc.intervalStartData and odc.intervalStopData
            
            % contents of each cell will be nTrials x nMaxOccurrence
            [intervalStartData, intervalStopData] = deal(cell(ad.nIntervals, 1));
            
            intervalStartCounts = nan(ad.nTrials, ad.nIntervals);
            intervalStopCounts = nan(ad.nTrials, ad.nIntervals);
            
            % collect the interval times
            for iInt = 1:ad.nIntervals
                [intervalStartData{iInt}, intervalStartCounts(:, iInt)] = ad.getEventIndexedTimeMatrix(...
                    ad.intervalEventsStart{iInt}, ad.intervalEventsIndexStart{iInt}, ...
                    ad.intervalOffsetsStart(iInt), ad.getZeroByTrial());
                [intervalStopData{iInt}, intervalStopCounts(:, iInt)] = ad.getEventIndexedTimeMatrix(...
                    ad.intervalEventsStop{iInt}, ad.intervalEventsIndexStop{iInt},  ...
                    ad.intervalOffsetsStop(iInt), ad.getZeroByTrial());
                
                intervalStartData{iInt}(~ad.computedValid, :) = NaN;
                intervalStopData{iInt}(~ad.computedValid, :) = NaN;
            end
            
            mismatch = intervalStartCounts ~= intervalStopCounts;
            if any(mismatch(:))
                warning('Encountered mismatched interval start / stop events');
            end
            intervalCounts = min(intervalStartCounts, intervalStopCounts);
            intervalCounts(~ad.computedValid, :) = NaN;
            
            c = ad.odc;
            c.intervalStartData = intervalStartData;
            c.intervalStopData = intervalStopData;
            c.intervalCounts = intervalCounts;
        end
    end

    methods % Dependent properties and post-applied alignment data request
        function assertApplied(ad)
            assert(ad.applied, 'Call .applyToTrialData first');
        end
        
        function nt = get.nTrials(ad)
            if ~ad.applied
                nt = 0;
            else
                nt = numel(ad.eventCounts.(ad.zeroEvent));
            end
       end
       
       function lengths = getValidDurationByTrial(ad) 
            ad.assertApplied();
            ti = ad.timeInfo;
            lengths = makecol([ti.stop] - [ti.start] + 1);
       end
        
       function [start, stop, zero] = getStartStopZeroByTrial(ad)
           ad.assertApplied();
           
           start = makecol([ad.timeInfo.start]);
           stop = makecol([ad.timeInfo.stop]);
           zero = makecol([ad.timeInfo.zero]);
       end
       
       function [startPad, stopPad, zero] = getStartStopZeroByTrialWithPadding(ad)
           ad.assertApplied();

           startPad = makecol([ad.timeInfo.startPad]);
           stopPad = makecol([ad.timeInfo.stopPad]);
           zero = makecol([ad.timeInfo.zero]);
       end

       % note that this should return NaNs for invalid trials
       function zero = getZeroByTrial(ad)
            ad.assertApplied();
            zero = makecol([ad.timeInfo.zero]);
       end
    end
    
    methods % Trial validity
        function ad = markInvalid(ad, invalid)
            ad.warnIfNoArgOut(nargout);
            ad.manualInvalid(invalid) = false;
        end

        function ad = setInvalid(ad, mask)
            ad.warnIfNoArgOut(nargout);
            ad.manualInvalid = mask;
        end
        
        function valid = get.valid(ad)
            if isempty(ad.timeInfo)
                valid = [];
            else
                valid = ad.computedValid & ~ad.manualInvalid;
            end
        end 
        
        function ad = selectTrials(ad, mask)
            ad.warnIfNoArgOut(nargout);
            
            flds = fieldnames(ad.eventData);
            for iF = 1:numel(flds)
                fld = flds{iF};
                ad.eventData.(fld) = ad.eventData.(fld)(mask, :);
                ad.eventCounts.(fld) = ad.eventCounts.(fld)(mask);
            end
            
            flds = fieldnames(ad.timeInfo);
            for iF = 1:numel(flds)
                fld = flds{iF};
                ad.timeInfo.(fld) = ad.timeInfo.(fld)(mask);
            end
            
            for iM = 1:ad.nMarks
                ad.markData{iM} = ad.markData{iM}(mask, :);
                ad.markCounts = ad.markCounts(mask, :);
            end
            
            for iI = 1:ad.nIntervals
                ad.intervalStartData{iI} = ad.intervalStartData{iI}(mask, :);
                ad.intervalStopData{iI} = ad.intervalStopData{iI}(mask, :);
                ad.intervalCounts = ad.intervalCounts(mask, :);
            end
            
            ad.manualInvalid = ad.manualInvalid(mask);
            ad.computedValid = ad.computedValid(mask);
        end
    end
    
    methods % Time-Aligning data transformations
        function [alignedTimes, rawTimesMask] = getAlignedTimesMatrix(ad, rawTimesMatrix, includePadding)
            if nargin < 3
                includePadding = false;
            end

            assert(size(rawTimesMatrix, 1) == ad.nTrials, 'Size(1) must match nTrials');
            
            % filter the spikes within the window and recenter on zero
            if includePadding
                start = ad.timeInfo.startPad;
                stop = ad.timeInfo.stopPad;
            else
                start = ad.timeInfo.start;
                stop = ad.timeInfo.stop;
            end
            zero = ad.timeInfo.zero;
            
            rawTimesMask = bsxfun(@ge, rawTimesMatrix, start) & bsxfun(@le, rawTimesMatrix, stop);
            alignedTimes = bsxfun(@minus, rawTimesMatrix, zero);  

            rawTimesMask(~ad.valid, :) = false;
            alignedTimes(~ad.valid, :) = NaN;
        end
        % use the alignment to shift the times in rawTimesCell to be zero relative
        % and filter by time window determined by getTimeInfo for each trial
        % if includePadding is true, will additional times found in the padWindow, see .setPadWindow
        function [alignedTimes, rawTimesMask] = getAlignedTimesCell(ad, rawTimesCell, includePadding)
            if nargin < 3
                includePadding = false;
            end
            
            if isempty(rawTimesCell)
                alignedTimes = rawTimesCell;
                rawTimesMask = rawTimesCell;
                return;
            end

            assert(numel(rawTimesCell) == ad.nTrials, 'Size must match nTrials');
            
            % filter the spikes within the window and recenter on zero
            if includePadding
                start = ad.timeInfo.startPad;
                stop = ad.timeInfo.stopPad;
            else
                start = ad.timeInfo.start;
                stop = ad.timeInfo.stop;
            end
            zero = ad.timeInfo.zero;
            
            [alignedTimes, rawTimesMask] = deal(cell(ad.nTrials, 1));
            valid = ad.valid;
            for i = 1:ad.nTrials
                if valid(i)
                    raw = rawTimesCell{i};
                    rawTimesMask{i} = raw >= start(i) & raw <= stop(i);
                    alignedTimes{i} = raw(rawTimesMask{i}) - zero(i);
                end
            end
        end
        
        function [alignedData, alignedTime] = getAlignedTimeseries(ad, dataCell, timeCell, includePadding, varargin)
            [alignedTime, rawTimesMask] = ad.getAlignedTimesCell(timeCell, includePadding);
            alignedData = cellfun(@(data, mask) data(mask, :), dataCell, rawTimesMask, ...
                'UniformOutput', false, 'ErrorHandler', @(varargin) []);
        end
        
        function [startRel, stopRel] = getStartStopRelativeToZeroByTrial(ad)
           ad.assertApplied();
           [start, stop, zero] = ad.getStartStopZeroByTrial();
           startRel = start - zero;
           stopRel = stop - zero;
        end
       
        function [markData, markDataMask] = getAlignedMarkData(ad)
            [markData, markDataMask] = ...
                cellfun(@ad.getAlignedTimesMatrix, ad.markData, 'UniformOutput', false);
        end
        
        function [intervalStartData, intervalStopData, intervalStartMask, intervalStopMask] = ...
                getAlignedIntervalData(ad)
            [intervalStartData, intervalStartMask] = ...
                cellfun(@ad.getAlignedTimesMatrix, ad.intervalStartData, 'UniformOutput', false);
            [intervalStopData, intervalStopMask] = ...
                cellfun(@ad.getAlignedTimesMatrix, ad.intervalStopData, 'UniformOutput', false);
        end
    end
    
    methods % Drawing on data
        % annotate data time-series with markers according to the labels indicated
        % by this align descriptor
        %
        % N is the number of traces to be annotated
        % T is number of time points
        % D is data dimensionality e.g. 1 or 2 or 3)
        %
        % the sizes of timeInfo and data may be one of the following:
        %   one-trial per data trace:
        %     timeInfo is N x 1 struct vec
        %     timeData is N x T matrix or N x 1 cell of T_i vectors
        %     data is N x T x D matrix or N x 1 cell of T_i x D matrices
        %
        %   many-trials per data trace:
        %     timeInfo is N x 1 cell array of ? x 1 struct vecs 
        %     timeData is N x T matrix or N x 1 cell of T_i vectors
        %     data is N x T x D matrix or N x 1 cell of T_i x D matrices 
        %     for this, the median will be computed for each data trace in timeInfo{:} and plotted accordingly
        %     on each of the N groups of timeInfos on data(m, :, :)
        %
        function drawOnData(as, timeData, data, varargin)
            p = inputParser();
            p.addParamValue('drawLegend', false, @islogical);
            p.addParamValue('drawRange', false, @islogical);
            p.parse(varargin{:});
            
            

            hold on

            N = length(timeInfo);
            assert(isvector(timeInfo) && (isstruct(timeInfo) || iscell(timeInfo)), ...
                'timeInfo must be struct vector or cell vector');
            assert(iscell(data) || N == size(data, 1), 'Length of timeInfo must match size(data, 1)');
            assert(~iscell(data) || (isvector(data) && N == length(data)), 'Data length must match timeInfo');
            assert(iscell(timeData) || N == size(timeData, 1), 'Length of timeInfo must match size(timeData, 1)');
            assert(~iscell(timeData) || (isvector(timeData) && N == length(timeData)), 'TimeData length must match timeInfo');
            assert(iscell(timeData) == iscell(data), 'TimeData and Data must both be cells or both matrices');

            hleg = nan(size(timeInfo));
            legstr = cell(size(timeInfo));
            for i = 1:length(timeInfo)
                if iscell(timeInfo)
                    % each time info is for a single data
                    ti = timeInfo{i};
                else
                    ti = timeInfo(i);
                end
                labelInfo = ad.getLabelInfo(ti);
                
                if iscell(data)
                    tvec = timeData{i};
                    dmat = data{i};
                else
                    tvec = squeeze(timeData(i, :));
                    dmat = squeeze(data(i, :, :));
                end
                
                if ~isempty(dmat)
                    drawOnSingle(ti, tvec, dmat, labelInfo);
                end
            end
            
            if p.Results.drawLegend
                idx = 1;
                for iLabel = 1:length(labelInfo)
                    info = labelInfo(iLabel).info;
                    if ~labelInfo(iLabel).markData
                        continue;
                    end
                    hleg(idx) = plot(NaN, NaN, info.marker, 'MarkerFaceColor', info.color, ...
                        'MarkerEdgeColor', info.color, 'MarkerSize', info.size);
                    legstr{idx} = labelInfo(iLabel).name;
                    idx = idx + 1;
                end
                
                legend(hleg, legstr, 'Location', 'NorthEast');
                legend boxoff;
            end
            
            function drawOnSingle(timeInfo, timeVec, dmat, labelInfo)
                % timeInfo is a struct array or single struct
                % dmat is T x D matrix
                nDim = size(dmat, 2);

                for iLabel = 1:length(labelInfo)
                    if ~labelInfo(iLabel).markData
                        continue;
                    end
                    info = labelInfo(iLabel).info;
                    ind = find(floor(labelInfo(iLabel).time) == floor(timeVec), 1);
                    if isempty(ind), continue, end
                    dvec = dmat(ind, :);
                    extraArgs = {info.marker, 'MarkerFaceColor', info.color, ...
                            'MarkerEdgeColor', info.color, 'MarkerSize', info.size};
                    if nDim == 1
                        plot(timeVec(ind), dvec(1), extraArgs{:});
                    elseif nDim == 2
                        plot(dvec(1), dvec(2), extraArgs{:});
                    elseif nDim == 3
                        plot3(dvec(1), dvec(2), dvec(3), extraArgs{:});
                    end
                end
                labelInfo = struct('name', {}, 'time', {}, 'align', {}, 'info', {}, ...
                    'markOndmat', {}, 'fixed', {});
            end
        end
    end
    
    methods(Static) % Default accessor methods
        function [eventData, eventCounts] = defaultGetEventTimesFn(td, eventNameList)
            % default function for accessing event times
            % Returns a cell array of size nTrials x nEvents
            % with each containing all events for event i, trial j
            % unless either 'takeFirst', or 'takeLast' is true, in which case not
            % regardless, if no events are found, the cell will contain a NaN rather than be empty
            
            if isa(td, 'TrialData')
                eventData = td.eventData;
                eventCounts = td.eventCounts;
            else
                error('Unsupported trial data type, please specify .getEventTimesFn');
            end
        end
    end
end

