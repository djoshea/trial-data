classdef(ConstructOnLoad) AlignInfo < AlignDescriptor
% AlignInfo is a subclass of AlignDescriptor that may be bound to a set of trials
% where the events have actual timestamps by trial.

    properties(SetAccess=protected)
        % AlignInfoOnDemandCache instance
        odc
    end

    properties(SetAccess=protected)
        % this function maps (R, eventList) --> eventTimes array nTrials x nEvents
        getEventTimesFn = @AlignInfo.defaultGetEventTimesFn;

        % struct array of nTrials x 1 containing the absolute times of each event
        % for the trials, as returned by getEventTimesFn
        eventInfo 

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
        
        % nTrials x nIntervals cell array of interval start times
        intervalStartData
        
        % nTrials x nIntervals cell array of interval stop times
        intervalStopData
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
            ad.eventInfo = ad.requestEventInfo(R);
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
        function [eventInfo] = requestEventInfo(ad, R)
            eventList = ad.getEventList();
            
            % grab the event times for all needed events
            times = ad.getEventTimesFn(R, eventList);

            % times may be:
            % nTrials x nEvents cell array
            % nTrials x nEvents matrix
            % nTrials x 1 struct array with fields == eventList
            %
            % ultimately convert to the struct format
            if iscell(times)
                eventInfo = cell2struct(times, eventList, 2);
            elseif isnumeric(times)
                eventInfo = cell2struct(num2cell(times), eventList, 2);
            elseif isstruct(times)
                eventInfo = times;
            else
                error('AlignInfo.getEventTimesFn returned an invalid value');
            end

            % replace empty cells with NaN --> no longer necessary as this
            % is done with the accessors below
            %eventInfo = structReplaceEmptyValues(eventInfo, NaN);
        end
        
        function assertHasEvent(ad, eventName)
            assert(isfield(ad.eventInfo, eventName), ...
                'AlignInfo has no event named %s', eventName);
        end

        % internal utility functions for accessing specific event times
        
        % n must be a scalar, times is a numeric array
        function times = getEventNthTimeVector(ad, event, n, offset, roundRelativeTo)
            ad.assertHasEvent(event);
            
            if strcmp(n, 'end')
                fn = @(info) info.(event)(end);
            else
                fn = @(info) info.(event)(n);
            end
            times = makecol(arrayfun(fn, ad.eventInfo, ...
                'ErrorHandler', @(varargin) NaN));
            times = times + offset;
            
            if ~isempty(roundRelativeTo) && ad.roundTimes
                % shift event times to maintain an integer multiple of
                % timeDelta offset from .zero
                times = ad.roundTimesRelativeTo(times, roundRelativeTo);
            end
        end
        
        % n may be an index or a selector (e.g. '2:end')
        function timeCell = getEventIndexedTimeCell(ad, event, n, offset, roundRelativeTo)
            ad.assertHasEvent(event);
            
            % similar to above but returns cell array, and n may be be a
            % string of the form '1:2', '1:end', ':', etc
            if ~ischar(n) || strcmp(n, 'end')
                timeCell = num2cell(ad.getEventNthTimeVector(event, n, offset, roundRelativeTo));
            else
                % must have a colon, parse into tokens
                pat = '(?<end1>end)?(?<ind1>-?\d*)?:(?<end2>end)?(?<ind2>-?\d*)?';
                info = regexp(n, pat, 'names', 'once');

                if isempty(info)
                    error('Unable to parse event index %s(%s)', event, n);
                end

                % convert ind1, ind2 to doubles
                if isempty(info.ind1)
                    ind1 = 1;
                else
                    ind1 = str2double(info.ind1);
                end
                if isempty(info.ind2)
                    ind2 = 0;
                else
                    ind2 = str2double(info.ind2);
                end
                if isempty(info.end1)
                    if isempty(info.end2)
                        % no end found before or after the colon
                        if isempty(info.ind2)
                            fn = @(info) info.(event)(ind1:end);
                        else
                            fn = @(info) info.(event)(ind1:ind2);
                        end
                    else
                        % after the colon has an end, use ind2 as an offset
                        fn = @(info) info.(event)(ind1:end+ind2);
                    end
                else
                    % end found before the colon, use ind1 as offset
                    if isempty(info.end2)
                        % no end after the colon
                        if isempty(info.ind2)
                            fn = @(info) info.(event)(end+ind1:end);
                        else
                            fn = @(info) info.(event)(end+ind1:ind2);
                        end
                    else
                        % end after the colon, use ind2 as offset
                        fn = @(info) info.(event)(end+ind1:end+ind2);
                    end
                end
                        
                timeCell = arrayfun(fn, ad.eventInfo, ...
                    'ErrorHandler', @(varargin) [], 'UniformOutput', false);
                %timeCell = cellfun(@ad.eventTimeRoundFn, timeCell, 'UniformOutput', false);
            end
            
            timeCell = cellfun(@(x) x + offset, timeCell, 'UniformOutput', false);
            if ~isempty(roundRelativeTo)
                % shift event times to maintain an integer multiple of
                % timeDelta offset from .zero
                timeCell = ad.roundTimesRelativeTo(timeCell, roundRelativeTo);
            end
        end
        
        % same as above but empty cells are filled with NaN
        function timeCell = getEventIndexedTimeCellFillEmptyWithNaN(ad, event, n, offset, roundRelativeTo)
            timeCell  = ad.getEventIndexedTimeCell(event, n, offset, roundRelativeTo);
            emptyMask = cellfun(@isempty, timeCell);
            [timeCell{emptyMask}] = deal(NaN);
        end
        
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

            nTrials = numel(ad.eventInfo);

            t.valid = truevec(nTrials);
            t.invalidCause = cellvec(nTrials);

            % get zero alignment event without rounding
            t.zero= ad.getEventNthTimeVector(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset, []);
            noZero = isnan(t.zero);
            [t.invalidCause{noZero & t.valid}] = deal(sprintf('Missing zero event %s', ad.zeroUnabbreviatedLabel));
            t.valid(noZero) = false;
            
            t.trialStart = ad.getEventNthTimeVector('TrialStart', 1, 0, t.zero); 
            t.trialStop = ad.getEventNthTimeVector('TrialEnd', 1, 0, t.zero);

            % get start event
            t.start = ad.getEventNthTimeVector(ad.startEvent, ad.startEventIndex, ad.startOffset, t.zero);
            noStart = isnan(t.start);
            [t.invalidCause{noStart}] = deal(sprintf('Missing start event %s', ad.startUnabbreviatedLabel));
            t.valid(noStart & t.valid) = false;

            % get stop event
            t.stop = ad.getEventNthTimeVector(ad.stopEvent, ad.stopEventIndex, ad.stopOffset, t.zero);
            noStop = isnan(t.stop);
            [t.invalidCause{noStop & t.valid}] = deal(sprintf('Missing stop event %s', ad.stopUnabbreviatedLabel));
            t.valid(noStop) = false;
        
            % get pad window
            t.startPad = t.start - padPre;
            t.stopPad = t.stop + padPost;
            
            % truncate trial end (including padding) based on truncateAfter events
            t.isTruncatedStop = falsevec(nTrials);
            for i = 1:length(ad.truncateAfterEvents)
                times = ad.getEventNthTimeVector(ad.truncateAfterEvents{i}, ad.truncateAfterEventsIndex{i}, ad.truncateAfterOffsets(i), t.zero);
                t.isTruncatedStop = t.isTruncatedStop | times < t.stopPad;
                t.stopPad = nanmin(t.stop, times);
                t.stop = t.stopPad - padPost;
            end
            
            % truncate trial start (including padding) based on truncateBefore events
            t.isTruncatedStart = falsevec(nTrials);
            for i = 1:length(ad.truncateBeforeEvents)
                times = ad.getEventNthTimeVector(ad.truncateBeforeEvents{i}, ad.truncateBeforeEventsIndex{i}, ad.truncateBeforeOffsets(i), t.zero);
                t.isTruncatedStart = t.isTruncatedStart | times > t.startPad;
                t.startPad = nanmax(t.startPad, times);
                t.start = t.startPad + padPre;
            end

            % mark trials as invalid if startPad:stopPad includes any invalidateEvents
            for i = 1:length(ad.invalidateEvents)
                timesCell = ad.getEventIndexedTimeCellFillEmptyWithNaN(ad.invalidateEvents{i}, ad.invalidateEventsIndex{i}, ad.invalidateOffsets(i), t.zero);
                maskInvalid = falsevec(length(t.startPad));
                for iT = 1:length(t.startPad)
                    maskInvalid(iT) = any(timesCell{iT} > t.startPad(iT) & timesCell{iT} < t.stopPad(iT));
                end
                t.valid(maskInvalid) = false;
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
                    t.valid(mask) = false;
                    [t.invalidCause{mask}] = deal('Start plus padding occurs before trial start');
                    
                    mask = t.stopPad > t.trialStop;
                    t.valid(mask) = false;
                    [t.invalidCause{mask}] = deal('Stop plus padding occurs after trial end');
                end
            end

            % handle minimum duration window
            mask = t.stop - t.start < ad.minDuration;
            t.valid(mask) = false;
            [t.invalidCause{mask}] = deal(sprintf('Trial duration is less than minDuration %g', ad.minDuration));
        
            % clear out values for invalid trials to avoid hard to catch bugs
            t.startPad(~t.valid) = NaN;
            t.stopPad(~t.valid) = NaN;
            t.start(~t.valid) = NaN;
            t.stop(~t.valid) = NaN;
            t.zero(~t.valid) = NaN;
            
            % now build the final timeInfo struct
            valid = t.valid;
            %t = rmfield(t, 'valid');
            timeInfo = makecol(structOfArraysToStructArray(t)); %#ok<*PROP>
            
            c = ad.odc;
            c.timeInfo = timeInfo;
            c.computedValid = valid;
        end
        
        function buildMarkData(ad)
            % compute the mark times and store in odc.markInfo
            
            markData = cell(ad.nTrials, ad.nMarks);
            for i = 1:length(ad.markEvents)
                markData(:, i) = ...
                    ad.getEventIndexedTimeCell(ad.markEvents{i}, ...
                    ad.markEventsIndex{i}, ad.markOffsets(i), ad.getZeroByTrial());
            end
            
            c = ad.odc;
            c.markData = markData;
        end
        
        function buildIntervalData(ad)
            % compute the interval times by trial and store in
            % odc.intervalStartData and odc.intervalStopData
            [intervalStartData, intervalStopData] = deal(cell(ad.nTrials, ad.nIntervals));
            
            % include the interval times
            for iInt = 1:ad.nIntervals
                intervalStartData(:, iInt) = ad.getEventIndexedTimeCell(...
                    ad.intervalEventsStart{iInt}, ad.intervalEventsIndexStart{iInt}, ...
                    ad.intervalOffsetsStart(iInt), ad.getZeroByTrial());
                intervalStopData(:, iInt) = ad.getEventIndexedTimeCell(...
                    ad.intervalEventsStop{iInt}, ad.intervalEventsIndexStop{iInt},  ...
                    ad.intervalOffsetsStop(iInt), ad.getZeroByTrial());
            end
            
            c = ad.odc;
            c.intervalStartData = intervalStartData;
            c.intervalStopData = intervalStopData;
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
                nt = numel(ad.timeInfo);
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
            ad = ad.updateSummary();
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
            ad.eventInfo = ad.eventInfo(mask);
            ad.timeInfo = ad.timeInfo(mask);
            ad.manualInvalid = ad.manualInvalid(mask);
            ad.computedValid = ad.computedValid(mask);
        end
    end
    
    methods % Time-Aligning data transformations
        
        % use the alignment to shift the times in rawTimesCell to be zero relative
        % and filter by time window determined by getTimeInfo for each trial
        % if includePadding is true, will additional times found in the padWindow, see .setPadWindow
        function [alignedTimes, rawTimesMask] = getAlignedTimes(ad, rawTimesCell, includePadding)
            if nargin < 3
                includePadding = false;
            end
            
            if isempty(rawTimesCell)
                alignedTimes = rawTimesCell;
                rawTimesMask = rawTimesCell;
                return;
            end
            
            % filter the spikes within the window and recenter on zero
            if includePadding
                start = num2cell([ad.timeInfo.startPad]);
                stop = num2cell([ad.timeInfo.stopPad]);
            else
                start = num2cell([ad.timeInfo.start]);
                stop = num2cell([ad.timeInfo.stop]);
            end
            
            repMatMultiplier = size(rawTimesCell);
            repMatMultiplier(1) = 1;
            repFn = @(z) repmat(makecol(z), repMatMultiplier);
            
            [alignedTimes, rawTimesMask] = cellfun(@fn, ...
                    makecol(rawTimesCell), ...
                    repFn(start), ...
                    repFn(stop), ...
                    repFn(num2cell([ad.timeInfo.zero])), ...
                    'UniformOutput', false);
                
            alignedTimes(~ad.valid) = {[]};
            
            function [alignedTimes, mask] = fn(rawTimes, tStart, tEnd, tZero)
                mask = rawTimes >= tStart & rawTimes <= tEnd;
                alignedTimes = rawTimes(mask) - tZero;
            end
        end
        
        function [alignedData, alignedTime] = getAlignedTimeseries(ad, dataCell, timeCell, includePadding, varargin)
            [alignedTime, rawTimesMask] = ad.getAlignedTimes(timeCell, includePadding);
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
            [markData, markDataMask] = ad.getAlignedTimes(ad.markData);
        end
        
        function [intervalStartData, intervalStopData, intervalStartMask, intervalStopMask] = getAlignedIntervalData(ad)
            [intervalStartData, intervalStartMask] = ad.getAlignedTimes(ad.intervalStartData);
            [intervalStopData, intervalStopMask] = ad.getAlignedTimes(ad.intervalStopData);
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
        function [timeData] = defaultGetEventTimesFn(R, eventNameList)
            % default function for accessing event times
            % Returns a cell array of size nTrials x nEvents
            % with each containing all events for event i, trial j
            % unless either 'takeFirst', or 'takeLast' is true, in which case not
            % regardless, if no events are found, the cell will contain a NaN rather than be empty

            if ~iscell(eventNameList)
                eventNameList = {eventNameList};
            end
            
            if isstruct(R)
                % keep only the appropriate fields
                timeData = keepfields(R, eventNameList);

            elseif isa(R, 'TrialData')
                timeData = R.getRawEventStruct();
               
            else
                error('Unsupported trial data type, please specify .getEventTimesFn');
            end
        end
    end
end

