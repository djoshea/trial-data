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
        
        timeInfoValid
        
        % valid is a dependent property formed by merging computedValid
        % with manualInvalid
        computedValid
        
        % nTrials x nMarks cell array of mark times
        markData
        
        % nTrials x nMarks cell array of mark times; NaN for invalid trials
        markDataValid
        
        % nTrials x nMarks array of mark occurrence counts
        markCounts
        
        % nTrials x nMarks array of mark occurrence counts; NaN for invalid trials
        markCountsValid
        
        % nMarks vector of maximum number of occurrences for each mark
        markMaxCounts
        
        % nIntervals cell array of nTrials x nMaxOccurrences interval start times
        intervalStartData
        
        % nIntervals cell array of nTrials x nMaxOccurrences interval start times
        % NaN for invalid trials
        intervalStartDataValid
        
        % nIntervals cell array of nTrials x nMaxOccurrences interval stop times
        intervalStopData
        
        intervalStopDataValid
        
        % nTrials x nIntervals array of interval occurrence counts
        intervalCounts
        
        % nTrials x nIntervals array of interval occurrence counts, NaN for
        % invalid trials
        intervalCountsValid
        
        % nMarks vector of maximum number of occurrences for each mark
        intervalMaxCounts
    end
    
    properties(Dependent)
        valid
        nTrials
        invalidCause
    end
    
    methods % Property get/set stored in odc
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
        
        function v = get.timeInfoValid(ad)
            v = ad.odc.timeInfoValid;
            if isempty(v)
                ad.buildTimeInfoValid();
                v = ad.odc.timeInfoValid;
            end
        end
        
        function ad = set.timeInfoValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.timeInfoValid = v;
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
        
        function v = get.markDataValid(ad)
            v = ad.odc.markDataValid;
            if isempty(v)
                ad.buildMarkDataValid();
                v = ad.odc.markDataValid;
            end
        end
        
        function ad = set.markDataValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.markDataValid = v;
        end
        
        function v = get.markCountsValid(ad)
            v = ad.odc.markCountsValid;
            if isempty(v)
                ad.buildMarkDataValid();
                v = ad.odc.markCountsValid;
            end
        end
        
        function ad = set.markCountsValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.markCountsValid = v;
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
        
        function v = get.intervalStartDataValid(ad)
            v = ad.odc.intervalStartDataValid;
            if isempty(v)
                ad.buildIntervalDataValid();
                v = ad.odc.intervalStartDataValid;
            end
        end
        
        function ad = set.intervalStartDataValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalStartDataValid = v;
        end
        
        function v = get.intervalStopDataValid(ad)
            v = ad.odc.intervalStopDataValid;
            if isempty(v)
                ad.buildIntervalDataValid();
                v = ad.odc.intervalStopDataValid;
            end
        end
        
        function ad = set.intervalStopDataValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalStopDataValid = v;
        end
        
        function v = get.intervalCountsValid(ad)
            v = ad.odc.intervalCountsValid;
            if isempty(v)
                ad.buildIntervalDataValid();
                v = ad.odc.intervalCountsValid;
            end
        end
        
        function ad = set.intervalCountsValid(ad, v)
            ad.odc = ad.odc.copy();
            ad.odc.intervalCountsValid = v;
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
        
        function ad = updateManualInvalid(ad)
            ad.warnIfNoArgOut(nargout);
            ad.odc = ad.odc.copy();
            ad.odc.flushManualValid();
        end
        
        function ad = postUpdateMark(ad)
            ad.warnIfNoArgOut(nargout);
            
            if ~ad.applied
                % nothing to udpate if we haven't applied to trial data yet
                return;
            end
            
            ad.odc = ad.odc.copy();
            ad.odc.flushMarkData();
        end
        
        function ad = postUpdateInterval(ad)
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
        % eventData.eventName is nTrials x nOccurrencesMax
        function [eventData, eventCounts] = requestEventInfo(ad, td)
            eventList = ad.getEventList();
            
            if td.nTrials == 0
                % create empty struct with the requisite fields, mostly
                % useful for empty constructor
                eventData = struct();
                eventCounts = struct();
                fields = union({'TrialStart', 'TrialEnd'}, eventList);
                for iF = 1:numel(fields)
                    eventData.(fields{iF}) = zeros(0, 0);
                    eventCounts.(fields{iF}) = zeros(0, 1);
                end
            else
                % grab the event times for all needed events
                [eventData, eventCounts] = ad.getEventTimesFn(td, eventList);
                
            end
            
            if ~isfield(eventData, 'TrialStart')
                error('Missing event field TrialStart');
                %                 eventData.TrialStart = nanvec(td.nTrials);
                %                 eventCounts.TrialStart = zerosvec(td.nTrials);
            end
            if ~isfield(eventData, 'TrialEnd')
                error('Missing event field TrialEnd');
                %                 eventData.TrialEnd = nanvec(td.nTrials);
                %                 eventCounts.TrialEnd = zerosvec(td.nTrials);
            end
            
            if isfield(eventData, 'TimeZero')
                error('Event "TimeZero" is reserved');
            end
            
            for iE = 1:numel(eventList)
                if ~isfield(eventData, eventList{iE})
                    warning('Missing event field %s', eventList{iE});
                    eventData.(eventList{iE}) = nanvec(td.nTrials); eventCounts.(eventList{iE}) = zerosvec(td.nTrials);
                end
            end
            
            % create the "TimeZero" event as 0 on each trial
            eventData.TimeZero = zerosvec(td.nTrials);
            eventCounts.TimeZero = onesvec(td.nTrials);
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
            
            if ad.nTrials == 0
                times = zeros(0, 1);
                return;
            end
            
            timesMat = ad.eventData.(event);
            counts = ad.eventCounts.(event);
            
            if strcmp(n, 'end')
                % TODO implement 'end-1' type index
                inds = sub2ind(size(timesMat), (1:ad.nTrials)', max(1, counts));
                times = timesMat(inds);
                times(counts < 1) = NaN;
            else
                if ischar(n)
                    n = str2double(n);
                elseif isempty(n)
                    n = 1;
                end
                
                if n > size(timesMat, 2)
                    times = nanvec(ad.nTrials);
                else
                    times = timesMat(:, n);
                end
            end
            
            times = times + double(offset);
            
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
            
            if isempty(n)
                n = ':';
            end
            
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
            roundFn = @(times, ref) bsxfun(@plus, round(bsxfun(@minus, times, ref) / delta) * delta, ref);
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
            nTrials = ad.nTrials; %#ok<*PROP>
            
            valid = truevec(nTrials);
            t.invalidCause = cellvec(nTrials);
            t.invalidCause(:) = {''};
            
            % get zero alignment event without rounding
            t.zero = ad.getEventNthTimeVector(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset, []);
            noZero = isnan(t.zero);
            [t.invalidCause{noZero & valid}] = deal(sprintf('Missing zero event %s', ad.zeroUnabbreviatedLabel));
            valid(noZero & valid) = false;
            
            t.trialStart = ad.getEventNthTimeVector('TrialStart', 1, 0, t.zero);
            t.trialStop = ad.getEventNthTimeVector('TrialEnd', 1, 0, t.zero);
            
            % get start event
            t.start = ad.getEventNthTimeVector(ad.startEvent, ad.startEventIndex, ad.startOffset, t.zero);
            noStart = isnan(t.start);
            [t.invalidCause{noStart & valid}] = deal(sprintf('Missing start event %s', ad.startUnabbreviatedLabel));
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
            
            % mark trials invalid if start > stop
            maskInvalid = valid & t.start > t.stop;
            valid(maskInvalid) = false;
            [t.invalidCause{maskInvalid}] = deal(sprintf('Start event after stop event after truncation'));
            
            
            % mark trials as invalid if startPad:stopPad includes any invalidateEvents
            for i = 1:length(ad.invalidateEvents)
                timesMatrix = ad.getEventIndexedTimeMatrix(ad.invalidateEvents{i}, ...
                    ad.invalidateEventsIndex{i}, ad.invalidateOffsets(i), t.zero);
                
                maskInvalid = any(bsxfun(@ge, timesMatrix, t.startPad) & ...
                    bsxfun(@le, timesMatrix, t.stopPad), 2);
                
                valid(maskInvalid) = false;
                [t.invalidCause{maskInvalid}] = deal(sprintf('Overlapped invalidating event %s', ad.invalidateUnabbreviatedLabels{i}));
            end
            
            % regardless of outsideOfTrialMode, also store truncated times so that padded window fits within the trial
            t.startTruncated = max(t.trialStart, t.start);
            t.stopTruncated = min(t.trialStop, t.stop);
            
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
            mask = max(0, t.stop - t.start) < ad.minDuration; % ignore start after stop
            valid(mask) = false;
            [t.invalidCause{mask}] = deal(sprintf('Trial duration is less than minDuration %g', ad.minDuration));
            
            % store in odc
            c = ad.odc;
            c.timeInfo = t;
            c.timeInfoValid = [];
            c.computedValid = valid;
        end
        
        function buildTimeInfoValid(ad)
            % generates timeInfo, but inserts NaNs in manual valid to avoid
            % bugs down the line. timeInfo is rarely updated, timeInfoValid
            
            t = ad.timeInfo;
            valid = ad.valid;
            
            [t.invalidCause{ad.manualInvalid}] = deal(sprintf('Marked invalid manually'));
            
            % clear out values for invalid trials to avoid hard to catch bugs
            % Here we consider computedValid AND manualValid
            t.startPad(~valid) = NaN;
            t.stopPad(~valid) = NaN;
            t.startTruncated(~valid) = NaN;
            t.stopTruncated(~valid) = NaN;
            t.start(~valid) = NaN;
            t.stop(~valid) = NaN;
            t.zero(~valid) = NaN;
            
            c = ad.odc;
            c.timeInfoValid = t;
        end
        
        function buildMarkData(ad)
            % compute the mark times and store in odc.markInfo
            % DONT MARK AS NAN here, this happens in buildMarkDataValid
            markData = cell(ad.nMarks, 1);
            markCounts = nan(ad.nTrials, ad.nMarks);
            for i = 1:length(ad.markEvents)
                [markData{i}, markCounts(:, i)] = ...
                    ad.getEventIndexedTimeMatrix(ad.markEvents{i}, ...
                    ad.markEventsIndex{i}, ad.markOffsets(i), ad.timeInfo.zero);
            end
            
            c = ad.odc;
            c.markData = markData;
            c.markCounts = markCounts;
        end
        
        function buildMarkDataValid(ad)
            markDataValid = ad.markData;
            markCountsValid = ad.markCounts;
            for i = 1:length(ad.markEvents)
                markDataValid{i}(~ad.valid, :) = NaN;
            end
            markCountsValid(~ad.computedValid, :) = NaN;
            
            c = ad.odc;
            c.markDataValid = markDataValid;
            c.markCountsValid = markCountsValid;
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
                    ad.intervalOffsetsStart(iInt), ad.timeInfo.zero);
                [intervalStopData{iInt}, intervalStopCounts(:, iInt)] = ad.getEventIndexedTimeMatrix(...
                    ad.intervalEventsStop{iInt}, ad.intervalEventsIndexStop{iInt},  ...
                    ad.intervalOffsetsStop(iInt), ad.timeInfo.zero);
                
                %intervalStartData{iInt}(~ad.computedValid, :) = NaN;
                %intervalStopData{iInt}(~ad.computedValid, :) = NaN;
            end
            
            mismatch = intervalStartCounts ~= intervalStopCounts;
            if any(mismatch(:))
                warning('Encountered mismatched interval start / stop events');
            end
            intervalCounts = min(intervalStartCounts, intervalStopCounts);
            %intervalCounts(~ad.computedValid, :) = NaN;
            
            c = ad.odc;
            c.intervalStartData = intervalStartData;
            c.intervalStopData = intervalStopData;
            c.intervalCounts = intervalCounts;
        end
        
        function buildIntervalDataValid(ad)
            intervalStartDataValid = ad.intervalStartData;
            intervalStopDataValid = ad.intervalStopData;
            intervalCountsValid = ad.intervalCounts;
            
            for iInt = 1:ad.nIntervals
                intervalStartDataValid{iInt}(~ad.valid, :) = NaN;
                intervalStopDataValid{iInt}(~ad.valid, :) = NaN;
            end
            intervalCountsValid(~ad.valid, :) = NaN;
            
            c = ad.odc;
            c.intervalStartDataValid = intervalStartDataValid;
            c.intervalStopDataValid = intervalStopDataValid;
            c.intervalCountsValid = intervalCountsValid;
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
        
        function n = get.markMaxCounts(ad)
            % markCounts is nTrials x nMarks
            n = max(ad.markCounts, [], 1)';
        end
        
        function n = get.intervalMaxCounts(ad)
            % intervalCounts is nTrials x nMarks
            n = max(ad.intervalCounts, [], 1)';
        end
        
        function lengths = getValidDurationByTrial(ad)
            ad.assertApplied();
            ti = ad.timeInfoValid;
            lengths = makecol([ti.stop] - [ti.start]); % was +1 here, not sure why this is necessary important 20160418
            lengths(lengths < 0) = 0;
        end
        
        function [TrialStart, TrialStop] = getUnalignedTrialStartStopByTrial(ad)
            ad.assertApplied();
            TrialStart = makecol([ad.timeInfoValid.trialStart]);
            TrialStop = makecol([ad.timeInfoValid.trialStop]);
        end
        
        function [start, stop, zero] = getStartStopZeroByTrial(ad)
            ad.assertApplied();
            
            start = makecol([ad.timeInfoValid.start]);
            stop = makecol([ad.timeInfoValid.stop]);
            zero = makecol([ad.timeInfoValid.zero]);
        end
        
        function [startPad, stopPad, zero] = getStartStopZeroByTrialWithPadding(ad)
            ad.assertApplied();
            
            startPad = makecol([ad.timeInfoValid.startPad]);
            stopPad = makecol([ad.timeInfoValid.stopPad]);
            zero = makecol([ad.timeInfoValid.zero]);
        end
        
        function [start, stop, zero] = getStartStopZeroByTrialTruncated(ad)
            ad.assertApplied();
            
            start = makecol([ad.timeInfoValid.startTruncated]);
            stop = makecol([ad.timeInfoValid.stopTruncated]);
            zero = makecol([ad.timeInfoValid.zero]);
        end
        
        % note that this should return NaNs for invalid trials
        function zero = getZeroByTrial(ad)
            ad.assertApplied();
            zero = makecol([ad.timeInfoValid.zero]);
        end
        
        function [startRel, stopRel, zeroRel] = getStartStopZeroRelativeToTrialStartByTrial(ad)
            ad.assertApplied();
            
            startRel = makecol([ad.timeInfoValid.start]) - makecol(ad.timeInfoValid.trialStart);
            stopRel = makecol([ad.timeInfoValid.stop]) - makecol(ad.timeInfoValid.trialStart);
            zeroRel = makecol([ad.timeInfoValid.zero]) - makecol(ad.timeInfoValid.trialStart);
        end
        
        function [startRel, stopRel, zeroRel] = getStartStopZeroRelativeToTrialStartByTrialWithPadding(ad)
            ad.assertApplied();
            
            [start, stop, zero] = ad.getStartStopZeroByTrialWithPadding();
            startRel = start - makecol(ad.timeInfoValid.trialStart);
            stopRel = stop - makecol(ad.timeInfoValid.trialStart);
            zeroRel = zero - makecol(ad.timeInfoValid.trialStart);
        end
        
        function [startRel, stopRel] = getStartStopRelativeToZeroByTrial(ad)
            ad.assertApplied();
            [start, stop, zero] = ad.getStartStopZeroByTrial();
            startRel = start - zero;
            stopRel = stop - zero;
        end
        
        function [startRel, stopRel] = getStartStopRelativeToZeroByTrialWithPadding(ad)
            ad.assertApplied();
            [start, stop, zero] = ad.getStartStopZeroByTrialWithPadding();
            startRel = start - zero;
            stopRel = stop - zero;
        end
        
        function [startRel, stopRel] = getStartStopRelativeToZeroByTrialTruncated(ad)
            ad.assertApplied();
            [start, stop, zero] = ad.getStartStopZeroByTrialTruncated();
            startRel = start - zero;
            stopRel = stop - zero;
        end
        
    end
    
    methods % Trial validity
        function ad = markInvalid(ad, invalid)
            ad.warnIfNoArgOut(nargout);
            ad.manualInvalid(invalid) = false;
            ad = ad.updateManualInvalid();
            % IT IS PERFORMANCE-CRITICAL THAT WE ONLY UPDATE TIMEINFOVALID HERE,
            % do not call .update() here
        end
        
        function ad = setManualValidTo(ad, mask)
            ad.warnIfNoArgOut(nargout);
            assert(isvector(mask) && islogical(mask) && numel(mask) == ad.nTrials);
            ad.manualInvalid = ~mask;
            ad = ad.updateManualInvalid();
            % IT IS PERFORMANCE-CRITICAL THAT WE ONLY UPDATE TIMEINFOVALID HERE,
            % do not call .update() here
        end
        
        function valid = get.valid(ad)
            if isempty(ad.timeInfo)
                valid = [];
            else
                valid = ad.computedValid & ~ad.manualInvalid;
            end
        end
        
        function cause = get.invalidCause(ad)
            if isempty(ad.timeInfo)
                cause = {};
            else
                cause = ad.timeInfo.invalidCause;
            end
        end
        
        function ad = selectTrials(ad, mask)
            ad.warnIfNoArgOut(nargout);
            
            flds = fieldnames(ad.eventData);
            for iF = 1:numel(flds)
                fld = flds{iF};
                ad.eventData.(fld) = ad.eventData.(fld)(mask, :);
                ad.eventCounts.(fld) = ad.eventCounts.(fld)(mask, :);
            end
            
            ad.timeInfo = [];
            ad.timeInfoValid = [];
            %             flds = fieldnames(ad.timeInfo);
            %             for iF = 1:numel(flds)
            %                 fld = flds{iF};
            %                 ad.timeInfo.(fld) = ad.timeInfo.(fld)(mask);
            %                 ad.timeInfoValid.(fld) = ad.timeInfoValid.(fld)(mask);
            %             end
            %
%             if ~isempty(ad.odc.markData)
%                 ad.markCounts = ad.markCounts(mask, :);
%                 nMax = max(ad.markCounts, [], 1)';
%                 for iM = 1:ad.nMarks
%                     ad.markData{iM} = ad.markData{iM}(mask, 1:nMax(iM));
%                 end
%             end
%             
%             if ~isempty(ad.odc.intervalStartData)
%                 ad.intervalCounts = ad.intervalCounts(mask, :);
%                 nMax = max(ad.intervalCounts, [], 1)';
%                 for iI = 1:ad.nIntervals
%                     ad.intervalStartData{iI} = ad.intervalStartData{iI}(mask, 1:nMax(iI));
%                     ad.intervalStopData{iI} = ad.intervalStopData{iI}(mask, 1:nMax(iI));
%                 end
%             end

            ad.manualInvalid = ad.manualInvalid(mask);
            ad = ad.update();
%             ad.computedValid = ad.computedValid(mask);
        end
    end
    
    methods % Time-Aligning data transformations
        function [rawTimesMask, indFirst, indLast, sampleDelta] = getAlignedTimesMask(ad, rawTimes, varargin)
            p = inputParser();
            p.addOptional('includePadding', false, @islogical);
            p.addParameter('isAligned', false, @islogical); % true means times already are relative to zero
            p.addParameter('singleTimepointTolerance', NaN, @isscalar); % acceptable tolerance to take sample when start == stop, NaN --> 2 * sampling rate
            p.addParameter('edgeTolerance', NaN, @isscalar); % acceptable tolerance to be within the interval, NaN --> 1e-6 * sampling rate
            p.addParameter('raw', false, @islogical); % ignore current .valid and operate on all trials that are align valid
            p.addParameter('preserveNativeScaling', 0, @isscalar); % if nonzero, this will be timeScaling to ms for the input data
            p.parse(varargin{:});
            request_raw = p.Results.raw;
            
            includePadding = p.Results.includePadding;
            singleTimepointTolerance = p.Results.singleTimepointTolerance;
            edgeTolerance = p.Results.edgeTolerance;
            preserveNativeScaling = p.Results.preserveNativeScaling;
            if preserveNativeScaling > 0
                timeScaling = double(preserveNativeScaling);
                nativeClass = TrialDataUtilities.Data.getCellElementClass(rawTimes);

                % with preserve native scaling, we can be a bit more precise about the edgeTolerance here
                if isnan(edgeTolerance)
                    edgeTolerance = double(ones(1, nativeClass)) * 1e-3;
                end
            else
                timeScaling = 1;
            end
            
            if nargout > 3 || isnan(singleTimepointTolerance) || isnan(edgeTolerance)
                sampleDelta = TrialData.computeDeltaFromTimes(rawTimes);
                if isnan(sampleDelta)
                    sampleDelta = 1;
                end
            end
            if isnan(singleTimepointTolerance)
                singleTimepointTolerance = 2 * sampleDelta;
            end
            if isnan(edgeTolerance)
                edgeTolerance = 1e-6 * sampleDelta;
            end
            
            if p.Results.raw
                timeInfo = ad.timeInfo;
            else
                timeInfo = ad.timeInfoValid;
            end
            
            % filter the times within the window
            if includePadding
                start = timeInfo.startPad;
                stop = timeInfo.stopPad;
            else
                start = timeInfo.start;
                stop = timeInfo.stop;
            end
            if p.Results.isAligned
                % add in time zero to unalign the incoming times
                offsets = timeInfo.zero;
            else
                offsets = zeros(ad.nTrials, 1);
            end
            
            valid = ad.valid; %#ok<*PROPLC>
            
            if iscell(rawTimes)
                % nTrials x J cell of vectors
                assert(size(rawTimes, 1) == ad.nTrials, 'rawTimes cell vector length must match nTrials');
                J = size(rawTimes, 2);
                rawTimesMask = cell(size(rawTimes));
            else
                % nTrials x T matrix of time
                assert(ismatrix(rawTimes) && size(rawTimes, 1) == ad.nTrials, 'Size(1) must match nTrials');
                J = 1;
                rawTimesMask = false(size(rawTimes));
            end
            
            [indFirst, indLast] = deal(nan(ad.nTrials, J));
            for i = 1:ad.nTrials
                if request_raw
                    valid_this = ~any(isnan([start(i), stop(i), offsets(i)])); 
                else
                    valid_this = valid(i);
                end
                if valid_this
                    for j = 1:J
                        if iscell(rawTimes)
                            raw = rawTimes{i, j};
                        else
                            raw = rawTimes(i, :)';
                        end

                        raw = double(raw)*timeScaling + offsets(i);
                        rawMask = false(numel(raw), 1);
                        
                        if (stop(i) - start(i)) < edgeTolerance && (stop(i) - start(i)) > -edgeTolerance % second clause is just in case somehow stop ends up pre-start
                            % only one timepoint requested, get the closest
                            % point if it's nearby
                            % provided that the sampling time is within the
                            % trial start:stop boundaries
                            if singleTimepointTolerance > 0 && start(i) + edgeTolerance >= timeInfo.trialStart(i) && ...
                                    stop(i) - edgeTolerance <= timeInfo.trialStop(i)
                                
                                [closestTime, closestIdx] = min(abs(raw - start(i)));
                                rawMask(closestIdx) = closestTime < singleTimepointTolerance;
                            end
                        else
                            % for time intervals, take only that interval
                            rawMask = raw >= start(i) - edgeTolerance & raw <= stop(i) + edgeTolerance;
                        end
                        
                        if iscell(rawTimes)
                            rawTimesMask{i, j} = rawMask;
                        else
                            rawTimesMask(i, :) = rawMask;
                        end
                        if any(rawMask)
                            indFirst(i, j) = find(rawMask, 1, 'first');
                            indLast(i, j) = find(rawMask, 1, 'last');
                        end
                    end
                else
                    % generate empty masks that match size
                    for j = 1:J
                        if iscell(rawTimes)
                            raw = rawTimes{i, j};
                        else
                            raw = rawTimes(i, :)';
                        end
                        rawTimesMask{i, j} = false(numel(raw), 1);
                    end
                end
            end
        end
        
        function [alignedTimes, rawTimesMask] = getAlignedTimesMatrix(ad, rawTimesMatrix, includePadding, beforeStartReplaceStart, afterStopReplaceStop)
            if nargin < 3
                includePadding = false;
            end
            if nargin < 4
                beforeStartReplaceStart = false;
            end
            if nargin < 5
                afterStopReplaceStop = false;
            end
            
            assert(size(rawTimesMatrix, 1) == ad.nTrials, 'Size(1) must match nTrials');
            
            % filter the spikes within the window and recenter on zero
            if includePadding
                start = ad.timeInfoValid.startPad;
                stop = ad.timeInfoValid.stopPad;
            else
                start = ad.timeInfoValid.start;
                stop = ad.timeInfoValid.stop;
            end
            zero = ad.timeInfoValid.zero;
            
            beforeStartMask = bsxfun(@lt, rawTimesMatrix, start);
            afterStopMask = bsxfun(@gt, rawTimesMatrix, stop);
            if beforeStartReplaceStart
                mask = beforeStartMask;
                startMat = repmat(start, size(rawTimesMatrix, 2));
                rawTimesMatrix(mask) = startMat(mask);
                beforeStartMask(:) = false;
            end
            if afterStopReplaceStop
                mask = afterStopMask;
                stopMat = repmat(stop, size(rawTimesMatrix, 2));
                rawTimesMatrix(mask) = stopMat(mask);
                afterStopMask(:) = false;
            end
            
            rawTimesMask = ~beforeStartMask & ~afterStopMask;
            alignedTimes = bsxfun(@minus, rawTimesMatrix, zero);
            
            rawTimesMask(~ad.valid, :) = false;
            alignedTimes(~rawTimesMask) = NaN;
        end
        
        % use the alignment to shift the times in rawTimesCell to be zero relative
        % and filter by time window determined by getTimeInfo for each trial
        % if includePadding is true, will additional times found in the padWindow, see .setPadWindow
        %
        % when an empty interval is sampled, the argument
        % simpleTimepointTolerance is important. When sampling an analog
        % simple at a specific timepoint, set this to some small positive
        % value and the closest sample in time will be taken provided it is
        % within singleTimepointTolerance of the provided value
        function [alignedTimes, rawTimesMask] = getAlignedTimesCell(ad, rawTimesCell, varargin)
            p = inputParser();
            p.addOptional('includePadding', false, @islogical);
            p.addParameter('isAligned', false, @islogical); % true means times already are relative to zero
            p.addParameter('singleTimepointTolerance', NaN, @isscalar); % acceptable tolerance to take sample when start == stop, NaN --> 2 * sampling rate
            p.addParameter('edgeTolerance', NaN, @isscalar); % acceptable tolerance to be within the interval, NaN --> 1e-6 * sampling rate
            p.addParameter('preserveNativeScaling', 0, @isscalar); % if nonzero, this will be the 
            p.parse(varargin{:});
            preserveNativeScaling = p.Results.preserveNativeScaling;

            if preserveNativeScaling > 0
                timeScaling = preserveNativeScaling;
            else
                timeScaling = 1;
            end

            if isempty(rawTimesCell)
                alignedTimes = rawTimesCell;
                rawTimesMask = rawTimesCell;
                return;
            end
            
            if ~p.Results.isAligned
                % subtract off time zero to align the incoming times
                offsets = ad.timeInfoValid.zero;
            else
                offsets = zeros(ad.nTrials, 1);
            end
            nativeClass = TrialDataUtilities.Data.getCellElementClass(rawTimesCell);
            signedClass = TrialDataUtilities.Data.getSignedDataType(nativeClass);
            
            assert(size(rawTimesCell, 1) == ad.nTrials, 'Size must match nTrials');
            
            J = size(rawTimesCell(:, :), 2);
            
            rawTimesMask = ad.getAlignedTimesMask(rawTimesCell, p.Results);
            
            alignedTimes = cell(size(rawTimesCell));

            for i = 1:ad.nTrials
                for j = 1:J
                    % subtract off the zero if needed
                    alignedTimes{i,j} = cast(rawTimesCell{i, j}(rawTimesMask{i, j}), signedClass) - ...
                        cast(offsets(i) / timeScaling, signedClass);
                end
            end
        end
        
        function intervalCell = getAlignedIntervalCell(ad, intervalCell, includePadding)
            % take data in interval cell (nTrials x 1) with nIntervals x 2
            % start stop matrices inside
            if nargin < 2
                includePadding = false;
            end
            assert(size(intervalCell, 1) == ad.nTrials);
            
            if includePadding
                start = ad.timeInfoValid.startPad;
                stop = ad.timeInfoValid.stopPad;
            else
                start = ad.timeInfoValid.start;
                stop = ad.timeInfoValid.stop;
            end
            zero = ad.timeInfoValid.zero;
            
            nEachTrial = size(intervalCell(:, :), 2);
            for iT = 1:ad.nTrials
                for j = 1:nEachTrial
                    % filter the spikes within the window and recenter on zero
                    rawTimesMatrix = intervalCell{iT, j};
                    
                    removeMask = falsevec(size(rawTimesMatrix, 1));
                    for iR = 1:size(rawTimesMatrix, 1)
                        if any(isnan(rawTimesMatrix(iR, :)))
                            removeMask(iR) = true;
                            
                        elseif rawTimesMatrix(iR, 2) < start(iT) || rawTimesMatrix(iR, 1) > stop(iT)
                            % not within start:stop, remove
                            removeMask(iR) = true;
                            
                        else
                            % constrain the edges within start:stop
                            rawTimesMatrix(iR, 1) = max(start(iT), rawTimesMatrix(iR, 1));
                            rawTimesMatrix(iR, 2) = min(stop(iT),  rawTimesMatrix(iR, 2));
                        end
                    end
                    
                    intervalCell{iT, j} = rawTimesMatrix(~removeMask, :) - zero(iT);
                end
            end
        end
        
        function [alignedData, alignedTime] = getAlignedTimeseries(ad, dataCell, timeCell, includePadding, varargin)
            % align timeseries where trials are along dimension 1
            [alignedTime, rawTimesMask] = ad.getAlignedTimesCell(timeCell, includePadding, varargin{:});
            
            %             nonEmpty = ~cellfun(@isempty, rawTimesMask);
            %             alignedData = cell(size(dataCell));
            %             alignedData(nonEmpty) = cellfun(@(data, mask) data(mask(1:size(data, 1)), :, :, :), ...
            %                 makecol(dataCell(nonEmpty)), makecol(rawTimesMask(nonEmpty)), ...
            %                 'UniformOutput', false);
            
            alignedData = cellfun(@(data, mask) data(mask(1:size(data, 1)), :, :, :), ...
                makecol(dataCell), makecol(rawTimesMask), ...
                'UniformOutput', false);
        end
        
        function [markData, markDataMask] = getAlignedMarkData(ad, includePadding)
            if nargin < 2
                includePadding = false;
            end
            [markData, markDataMask] = ...
                cellfun(@(m) ad.getAlignedTimesMatrix(m, includePadding), ad.markDataValid, 'UniformOutput', false);
        end
        
        function [intervalStartData, intervalStopData, intervalStartMask, intervalStopMask] = ...
                getAlignedIntervalData(ad, includePadding)
            if nargin < 2
                includePadding = false;
            end
            [intervalStartData, intervalStartMask] = ...
                cellfun(@(m) ad.getAlignedTimesMatrix(m, includePadding, true, false), ad.intervalStartDataValid, 'UniformOutput', false);
            [intervalStopData, intervalStopMask] = ...
                cellfun(@(m) ad.getAlignedTimesMatrix(m, includePadding, false, true), ad.intervalStopDataValid, 'UniformOutput', false);
        end
        
    end
    
    methods % Time-Aligning data transformations to individual marks and intervals
        % use the alignment to shift the times in rawTimesCell to be zero relative
        % to EACH mark and filter by time window determined by getTimeInfo for each trial
        % if includePadding is true, will additional times found in the padWindow, see .setPadWindow
        function [alignedTimes, rawTimesMask] = getMarkAlignedTimesCell(ad, rawTimesCell, markIdx, window, includePadding, varargin)
            if nargin < 5
                includePadding = false;
            end
            
            if ischar(markIdx)
                markIdx = ad.findMarkByString(markIdx);
            end
            assert(markIdx >= 1 && markIdx <= ad.nMarks, 'Mark idx out of range');
            
            if isempty(rawTimesCell)
                alignedTimes = rawTimesCell;
                rawTimesMask = rawTimesCell;
                return;
            end
            
            % grab the aligned mark data
            assert(numel(rawTimesCell) == ad.nTrials, 'Size must match nTrials');
            markDataAll = ad.getAlignedMarkData(includePadding);
            alignedMarkMat = markDataAll{markIdx};
            
            % align the raw
            [alignedTimesTrial, rawTimesMaskTrial] = ad.getAlignedTimesCell(rawTimesCell, includePadding, varargin{:});
            
            % filter the spikes within the window and recenter on zero
            startMat = alignedMarkMat + window(1);
            stopMat = alignedMarkMat + window(2);
            zeroMat = alignedMarkMat;
            
            nMarkOccur = ad.markMaxCounts(markIdx);
            [alignedTimes, rawTimesMask] = deal(cell(ad.nTrials, nMarkOccur));
            %valid = ad.valid;
            for i = 1:ad.nTrials
                selectedFromTrial = alignedTimesTrial{i}; % this includes only pre-selected times in trial
                maskWholeTrial = rawTimesMaskTrial{i}; % this is a mask over all times in trial
                for m = 1:nMarkOccur
                    % this is a sub-selection mask into selectedFromTrial
                    selected = selectedFromTrial >= startMat(i, m) & selectedFromTrial <= stopMat(i, m);
                    idxExclude = find(maskWholeTrial);
                    mask = maskWholeTrial;
                    mask(idxExclude(~selected)) = false;
                    rawTimesMask{i, m} = mask;
                    alignedTimes{i, m} = selectedFromTrial(selected) - zeroMat(i, m);
                end
            end
        end
        
    end
    
    methods % Drawing on single trial data
        function [hMarks, hIntervals] = drawOnDataByTrial(ad, varargin)
            % annotate data time-series with markers according to the labels indicated
            % by this align descriptor. Similar to AlignSummary's version
            % except operates on individual trials
            %
            % N = nTrials is the number of traces to be annotated
            % D is data dimensionality e.g. 1 or 2 or 3)
            %
            % the sizes of timeData and data may be one of the following:
            %   one-trial per data trace:
            %     timeData is T vector or N x 1 cell of T_i vectors
            %     data is N x T x D matrix or N x 1 cell of T_i x D matrices
            p = inputParser();
            p.addParameter('time', [], @(x) isnumeric(x) || iscell(x));
            p.addParameter('data', [], @(x) isnumeric(x) || iscell(x));
            p.addParameter('axh', gca, @ishandle);
            p.addParameter('tOffsetZero', 0, @isscalar);
            p.addParameter('showMarks', true, @islogical);
            p.addParameter('showIntervals', true, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markOutline', true, @islogical);
            p.addParameter('markOutlineAlpha', 1, @isscalar);
            p.addParameter('markSize', 8, @isscalar);
            p.addParameter('intervalThickness', 3, @isscalar);
            p.addParameter('intervalAlpha', 1, @isscalar);
            p.addParameter('trialIdx', 1:ad.nTrials, @isnumeric);
            p.addParameter('showInLegend', true, @islogical);
            p.addParameter('useTranslucentMark3d', false, @islogical);
            p.addParameter('clipping', 'on', @ischar);
            p.parse(varargin{:});
            
            data = p.Results.data;
            time = p.Results.time;
            
            axh = p.Results.axh;
            tOffsetZero = p.Results.tOffsetZero;
            trialIdx = p.Results.trialIdx;
            hold(axh, 'on');
            
            N = numel(trialIdx);
            assert(N == size(data, 1), 'size(data, 1) must match nTrials');
            %assert(N == size(time, 1), 'size(timeData, 1) must match nTrials');
            
            % grab information about the mark times relative to the trial
            nOccurByMark = ad.markMaxCounts;
            markData = ad.getAlignedMarkData();
            
            if iscell(data)
                emptyMask = cellfun(@isempty, data);
                if all(emptyMask), return; end
                nonEmpty = find(~emptyMask, 1);
                D = size(data{nonEmpty}, 2);
                markIntervalClipping = 'off';
            else
                D = size(data, 3);
                markIntervalClipping = 'on';
            end
            
            % plot intervals
            if p.Results.showIntervals
                hIntervals = cell(ad.nIntervals, 1);
                nOccurByInterval = ad.intervalMaxCounts;
                [intStartData, intStopData] = ad.getAlignedIntervalData();
                for iInterval = 1:ad.nIntervals
                    if ~ad.intervalShowOnData(iInterval), continue, end
                    
                    % gather mark locations
                    % nOccur x nTrials cell of T x D data in interval
                    intLoc = cell(nOccurByInterval(iInterval), N);
                    
                    for iTrial = 1:N
                        % filter by the time window specified (for this trial)
                        tStart = intStartData{iInterval}(trialIdx(iTrial), :)';
                        tStop = intStopData{iInterval}(trialIdx(iTrial), :)';
                        
                        % tvec should T vector, dmat should be T x D
                        if iscell(time)
                            tvec = time{iTrial};
                        else
                            tvec = time;
                        end
                        if iscell(data)
                            dmat = data{iTrial};
                        else
                            % data is N x T x D matrix
                            dmat = TensorUtils.squeezeDims(data(iTrial, :, :), 1);
                        end
                        
                        % constrain the time window to the interval being
                        % plotted as defined by tvec
                        valid = ~isnan(tStart) & ~isnan(tStop);
                        valid(tStart > max(tvec)) = false;
                        valid(tStop < min(tvec)) = false;
                        
                        if ~any(valid), continue; end
                        
                        tStart(tStart < min(tvec)) = min(tvec);
                        tStop(tStop > max(tvec)) = max(tvec);
                        
                        % and slice the interval location
                        % tStart, tStop is nOccur x 1
                        % dError will be nOccur cell with T x D values
                        dInterval = TrialDataUtilities.Plotting.DrawOnData.sliceIntervalLocations(tvec, dmat, tStart, tStop);
                        intLoc(:, iTrial) = dInterval;
                    end
                    
                    % add the time offset if plotting against time
                    if D == 1
                        for i = 1:numel(intLoc)
                            if isempty(intLoc{i}), continue; end
                            intLoc{i}(:, 1, :) = intLoc{i}(:, 1, :) + tOffsetZero;
                        end
                    end
                    
                    app = ad.intervalAppear{iInterval};
                    
                    hIntervals{iInterval} = TrialDataUtilities.Plotting.DrawOnData.plotInterval(axh, intLoc, D, ...
                        app, p.Results.intervalThickness, p.Results.intervalAlpha, 'clipping', p.Results.clipping);
                    if p.Results.showInLegend
                        TrialDataUtilities.Plotting.showFirstInLegend(hIntervals{iInterval}, ad.intervalLabels{iInterval});
                    else
                        TrialDataUtilities.Plotting.hideInLegend(hIntervals{iInterval});
                    end
                end
            else
                hIntervals = {};
            end
            
            if p.Results.showMarks
                % plot marks
                hMarks = cell(ad.nMarks, 1);
                for iMark = 1:ad.nMarks
                    if ~ad.markShowOnData(iMark), continue, end
                    
                    % gather mark locations
                    % nOccur x D x N
                    markLoc = nan(nOccurByMark(iMark), max(2, D), N);
                    
                    for t = 1:N
                        % get the mark times on this trial
                        tMark = markData{iMark}(trialIdx(t), :);
                        
                        if iscell(time)
                            tvec = time{t};
                        else
                            tvec = time;
                        end
                        if iscell(data)
                            dmat = data{t};
                        else
                            % data is N x T x D matrix
                            dmat = TensorUtils.squeezeDims(data(t, :, :), 1);
                        end
                        % tvec should T vector, dmat should be T x D
                        
                        % filter by the time window specified (for this trial)
                        maskInvalid = tMark < min(tvec) | tMark > max(tvec);
                        tMark(maskInvalid) = NaN;
                        
                        if all(isnan(tMark))
                            % none found in this time window for this condition
                            continue;
                        end
                        if isempty(dmat)
                            continue;
                        end
                        
                        % dMean will be nOccurThisTrial x max(2,D)
                        % since time will become dMean(:, 1, :) if D == 1
                        dMark = TrialDataUtilities.Plotting.DrawOnData.interpMarkLocation(tvec, dmat, tMark);
                        
                        markLoc(1:size(dMark, 1), :, t) = dMark;
                    end
                    
                    % add the time offset to time column if plotting against time
                    if D == 1
                        markLoc(:, 1, :) = markLoc(:, 1, :) + tOffsetZero;
                    end
                    
                    app = ad.markAppear{iMark};
                    
                    % plot mark and provide legend hint
                    if ~isempty(markLoc)
                        hMarks{iMark} = TrialDataUtilities.Plotting.DrawOnData.plotMark(axh, markLoc, app, ...
                            p.Results.markSize, 'alpha', p.Results.markAlpha, 'useTranslucentMark3d', p.Results.useTranslucentMark3d, ...
                            'outline', p.Results.markOutline, 'outlineAlpha', p.Results.markOutlineAlpha,  'clipping', p.Results.clipping);
                        
                        if p.Results.showInLegend
                            TrialDataUtilities.Plotting.showInLegend(hMarks{iMark}(1), ad.markLabels{iMark});
                            TrialDataUtilities.Plotting.hideInLegend(hMarks{iMark}(2:end));
                        else
                            TrialDataUtilities.Plotting.hideInLegend(hMarks{iMark});
                        end
                    end
                end
            else
                hMarks = {};
            end
        end
        
        function [hMarks, hIntervals] = drawOnRasterByTrial(ad, varargin)
            % annotate per-trial raster plots with markers according to the labels indicated
            % by this align descriptor.
            %
            % timeCell is nTrials x 1 cell of time vectors
            % it is assumed that each row is drawn below the last, such that timeCell{1} is drawn from
            % yOffsetTop to yOffsetTop + 1, with timeCell{2} drawn 1 below that
            % yOffsetTop refers to the top of the first trial's spikes
            p = inputParser();
            p.addParameter('startByTrial', [], @(x) isnumeric(x) && isvector(x));
            p.addParameter('stopByTrial', [], @(x) isnumeric(x) && isvector(x));
            p.addParameter('axh', gca, @ishandle);
            p.addParameter('tOffsetZero', 0, @isscalar);
            p.addParameter('yOffsetTop', 0, @isscalar);
            p.addParameter('tickHeight', 1, @isscalar);
            p.addParameter('intervalMinWidth', NaN, @isscalar); % if specified, draws intervals wider than they really are to be more readily visible
            
            p.addParameter('markAsTicks', true, @islogical);
            p.addParameter('markAlpha', 0.5, @isscalar);
            p.addParameter('markTickWidth', 2, @isscalar);
            p.addParameter('intervalAlpha', 0.5, @isscalar);
            p.addParameter('trialIdx', 1:ad.nTrials, @isnumeric);
            p.addParameter('showInLegend', true, @islogical);
            
            p.addParameter('shadeStartStopInterval', false, @islogical);
            p.addParameter('shadeOutsideStartStopInterval', false, @islogical);
            
            p.addParameter('showMarks', true, @islogical);
            p.addParameter('showIntervals', true, @islogical);
            
            p.addParameter('fullTimeLimits', [], @isnumeric);
            p.parse(varargin{:});
            
            axh = p.Results.axh;
            tOffsetZero = p.Results.tOffsetZero;
            yOffsetTop = p.Results.yOffsetTop;
            trialIdx = p.Results.trialIdx;
            startByTrial = p.Results.startByTrial;
            stopByTrial = p.Results.stopByTrial;
            
            hold(axh, 'on');
            
            N = numel(trialIdx);
            assert(N == size(startByTrial, 1), 'numel(startByTrial) must match nTrials');
            assert(N == size(stopByTrial, 1), 'numel(startByTrial) must match nTrials');
            
            % optionally shade start:stop interval to visualize valid data
            % region
            if p.Results.shadeStartStopInterval
                % must be row vectors because expects nOccurrences x
                % nTrials
                h = TrialDataUtilities.Plotting.DrawOnData.plotIntervalOnRaster(axh, startByTrial', stopByTrial', ...
                    AppearanceSpec('Color', [0.7 0.7 0.7]), p.Results.intervalAlpha, 'xOffset', tOffsetZero, 'yOffset', yOffsetTop, ...
                    'intervalHeight', p.Results.tickHeight, 'intervalMinWidth', 0);
                if p.Results.showInLegend
                    TrialDataUtilities.Plotting.showFirstInLegend(h, 'Valid time region');
                else
                    TrialDataUtilities.Plotting.hideInLegend(h);
                end
            end
            
            if p.Results.shadeOutsideStartStopInterval
                fullTimeLimits = p.Results.fullTimeLimits;
                if isempty(fullTimeLimits)
                    error('Must specify fullTimeLimits when shadeOutsideStartStopInterval is true');
                end
                % must be row vectors because expects nOccurrences x
                % nTrials
                h = TrialDataUtilities.Plotting.DrawOnData.plotIntervalOnRaster(axh, repmat(fullTimeLimits(1), 1, N), startByTrial', ...
                    AppearanceSpec('Color', [0.7 0.7 0.7]), p.Results.intervalAlpha, 'xOffset', tOffsetZero, 'yOffset', yOffsetTop, ...
                    'intervalHeight', p.Results.tickHeight, 'intervalMinWidth', 0);
                h2 = TrialDataUtilities.Plotting.DrawOnData.plotIntervalOnRaster(axh, stopByTrial', repmat(fullTimeLimits(2), 1, N), ...
                    AppearanceSpec('Color', [0.7 0.7 0.7]), p.Results.intervalAlpha, 'xOffset', tOffsetZero, 'yOffset', yOffsetTop, ...
                    'intervalHeight', p.Results.tickHeight, 'intervalMinWidth', 0);
                if p.Results.showInLegend
                    TrialDataUtilities.Plotting.showFirstInLegend(h, 'Invalid time region');
                else
                    TrialDataUtilities.Plotting.hideInLegend(h);
                end
                TrialDataUtilities.Plotting.hideInLegend(h2);
            end
            
            
            % plot intervals
            if p.Results.showIntervals
                hIntervals = cell(ad.nIntervals, 1);
                nOccurByInterval = ad.intervalMaxCounts;
                [intStartData, intStopData] = ad.getAlignedIntervalData();
                for iInterval = 1:ad.nIntervals
                    if ~ad.intervalShowOnData(iInterval), continue; end
                    % gather interval locations
                    % nOccur x nTrials set of start, stop by in interval
                    [intStart, intStop] = deal(nan(nOccurByInterval(iInterval), N));
                    
                    for iTrial = 1:N
                        % filter by the time window specified (for this trial)
                        tStart = intStartData{iInterval}(trialIdx(iTrial), :)';
                        tStop = intStopData{iInterval}(trialIdx(iTrial), :)';
                        
                        % constrain the time window to the interval being
                        % plotted as defined by tvec. not valid if it lies entirely outside
                        % the interval for this trial defined by start/stopByTrial
                        valid = ~isnan(tStart) & ~isnan(tStop);
                        valid(tStart > stopByTrial(iTrial)) = false;
                        valid(tStop < startByTrial(iTrial)) = false;
                        
                        if ~any(valid), continue; end
                        
                        tStart(tStart < startByTrial(iTrial)) = startByTrial(iTrial);
                        tStop(tStop > stopByTrial(iTrial)) = stopByTrial(iTrial);
                        
                        intStart(:, iTrial) = tStart;
                        intStop(:, iTrial) = tStop;
                    end
                    
                    app = ad.intervalAppear{iInterval};
                    
                    hIntervals{iInterval} = TrialDataUtilities.Plotting.DrawOnData.plotIntervalOnRaster(axh, intStart, intStop, ...
                        app, p.Results.intervalAlpha, 'xOffset', tOffsetZero, 'yOffset', yOffsetTop, ...
                        'intervalHeight', p.Results.tickHeight, 'intervalMinWidth', p.Results.intervalMinWidth);
                    
                    if p.Results.showInLegend
                        TrialDataUtilities.Plotting.showFirstInLegend(hIntervals{iInterval}, ad.intervalLabels{iInterval});
                    else
                        TrialDataUtilities.Plotting.hideInLegend(hIntervals{iInterval});
                    end
                end
            end
            
            if p.Results.showMarks
                % plot marks
                % grab information about the mark times relative to the trial
                nOccurByMark = ad.markMaxCounts;
                markData = ad.getAlignedMarkData();
                
                hMarks = cell(ad.nMarks, 1);
                for iMark = 1:ad.nMarks
                    if ~ad.markShowOnData(iMark), continue; end
                    % gather mark locations
                    % nOccur x N
                    markLoc = nan(nOccurByMark(iMark), N);
                    
                    for t = 1:N
                        % get the mark times on this trial
                        tMark = markData{iMark}(trialIdx(t), :);
                        
                        % filter by the time window specified (for this trial)
                        maskInvalid = tMark < startByTrial(t) | tMark > stopByTrial(t);
                        tMark(maskInvalid) = NaN;
                        
                        if all(isnan(tMark))
                            % none found in this time window for this condition
                            continue;
                        end
                        
                        markLoc(:, t) = tMark;
                    end
                    
                    app = ad.markAppear{iMark};
                    
                    % plot mark and provide legend hint
                    hMarks{iMark} = TrialDataUtilities.Plotting.DrawOnData.plotMarkOnRaster(axh, markLoc, app, ...
                        p.Results.markAlpha, 'tickWidth', p.Results.markTickWidth, ...
                        'xOffset', tOffsetZero, 'yOffset', yOffsetTop, ...
                        'tickHeight', p.Results.tickHeight);
                    
                    if p.Results.showInLegend
                        TrialDataUtilities.Plotting.showInLegend(hMarks{iMark}(1), ad.markLabels{iMark});
                        TrialDataUtilities.Plotting.hideInLegend(hMarks{iMark}(2:end));
                    else
                        TrialDataUtilities.Plotting.hideInLegend(hMarks{iMark});
                    end
                end
            end
        end
    end
    
    methods(Static) % Default accessor methods
        function [eventData, eventCounts] = defaultGetEventTimesFn(td, ~)
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
