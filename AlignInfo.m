classdef AlignInfo < AlignDescriptor
% AlignInfo is a subclass of AlignDescriptor that may be bound to a set of trials
% where the events have actual timestamps by trial.

    properties
        % marks not using the same event as zero will be surrounded by < > 
        % to indicate the mean time is being plotted. if false, they will
        % not be plotted at all.
        markPlotMedians = true;
        
        % for such marks, if the range of relative event times (in either 
        % direction is less than this threshold, the < > marks will be omitted.
        markRelativeDeltaIgnore = 7.5;

        % this function maps (R, eventList) --> eventTimes array nTrials x nEvents
        getEventTimesFn = @AlignInfo.defaultGetEventTimesFn;

        % struct array of nTrials x 1 containing the absolute times of each event
        % for the trials, as returned by getEventTimesFn
        eventInfo 

        % struct array of nTrials x 1 containing the absolute times of the
        % start, stop, zero events, as well as intervals and marks 
        timeInfo 

        % valid is a dependent property formed by merging these two
        computedValid
        manualInvalid

        applied = false; 
    end
    
    properties(Dependent)
        valid
        nTrials
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
            validStr = sprintf('(%d valid)', nnz(ad.computedValid));
            if ~ad.nameDefault
                tcprintf('inline', '{yellow}%s: {bright blue}%s : {none}%s %s\n', class(ad), ad.name, desc, validStr);
            else
                % name will just match desc, so don't print it twice
                tcprintf('inline', '{yellow}%s: {none}%s %s\n', class(ad), desc, validStr);
            end
        end

        function ad = update(ad)
            if ~ad.applied
                % nothing to udpate if we haven't applied to trial data yet
                return;
            end
            
            ad.warnIfNoArgOut(nargout);
            [ad.timeInfo, ad.computedValid] = ad.buildTimeInfo();
            %ad.summary = ad.summarizeTimeInfo(ad.timeInfo);
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

            % replace empty cells with NaN
            %eventInfo = structReplaceEmptyValues(eventInfo, NaN);
        end

        % internal utility functions for accessing specific event times
        
        % n must be a scalar, times is a numeric array
        function times = getEventNthTimeVector(ad, event, n, offset, roundRelativeTo)
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
            % similar to above but returns cell array, and n may be be a
            % string of the form '1:2', '1:end', ':', etc
            if ~ischar(n) || strcmp(n, 'end')
                timeCell = num2cell(ad.getEventNthTimeVector(event, n, offset));
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
                        fn = @(info) info.(event)(ind1:ind2);
                    else
                        fn = @(info) info.(event)(ind1:end+ind2);
                    end
                else
                    if isempty(info.end2)
                        fn = @(info) info.(event)(end+ind1:ind2);
                    else
                        fn = @(info) info.(event)(end+ind1:end+ind2);
                    end
                end
                        
                timeCell = arrayfun(fn, ad.eventInfo, ...
                    'ErrorHandler', @(varargin) [], 'UniformOutput', false);
                %timeCell = cellfun(@ad.eventTimeRoundFn, timeCell, 'UniformOutput', false);
            end
            
            timeCell = cellfun(@(x) x + offset, timeCell, 'UniformOutput', false);
            if ~isempty(roundRelativeTo) && ~isempty(roundRelativeTo)
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
            delta = ad.minTimeDelta;
            roundFn = @(times, ref) round((times - ref) / delta) * delta + ref;
            if isnumeric(times)
                times = roundFn(times, ref);
            else
                times = cellfun(roundFn, times, num2cell(ref), 'UniformOutput', false);
            end
        end
        
        % get the aligned start/stop/zero/mark time windows for each trial, 
        % respecting all truncation and invalidation instructions
        function [timeInfo, valid] = buildTimeInfo(ad)
            % returns a struct array with the actual time window and time of zero for trial i as
            %   timeInfo(i).start, .stop, .zero
            %
            % timeInfo(i).valid and valid(i) indicate whether trial i satisfied inclusion criteria
            % as specified by the various means of trial invalidation 
           
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
            timeInfo = structOfArraysToStructArray(t);
            nTrials = numel(timeInfo);
            
            % include the mark times
            for i = 1:length(ad.markEvents)
                markTimesCell = ad.getEventIndexedTimeVector(ad.markEvents{i}, ad.markEventsIndex{i}, ad.markOffsets(i));
                for iT = 1:nTrials
                    timeInfo(iT).mark{i} = markTimesCell{i};
                end
            end
            
            % include the interval times
            for iInt = 1:size(ad.intervalEventsStart, 1)
                startTimes = ad.getEventIndexedTimeCellFillEmptyWithNaN(ad.intervalEventsStart{iInt}, ad.intervalEventsIndexStart, ad.intervalOffsetsStart(iInt), t.zero);
                stopTimes = ad.getEventIndexedTimeCellFillEmptyWithNaN(ad.intervalEventsStop{iInt}, ad.intervalEventsIndexStop,  ad.intervalOffsetsStop(iInt), t.zero);
                for iTrial = 1:nTrials
                    timeInfo(iTrial).intervalStart{iInt} = startTimes{iTrial};
                    timeInfo(iTrial).intervalStop{iInt} = stopTimes{iTrial};
                end 
            end
        end
        
        function ad = updateSummary(ad)
            % look over the timeInfo struct and compute aggregate statistics about
            % the timing of each event relative to .zero
            %
            % .summaryInfo(i) looks like
            %     .name 
            %     .median
            %     .mean
            %     .list
            return;
            
%             events = ad.getEventList(); 
% 
%             zeroTimes = [ad.timeInfo.(ad.zeroEvent)];
% 
%             for iEv = 1:length(events)
%                 event = events{iEv};
% 
%                 times = ad.timeInfo;
% 
%                 evi.fixed = true;
%                 evi.relativeMedian = ad.startOffset - ad.zeroOffset;
%                 evi.relativeList = repmat(ad.nTrials, 1, evi.relativeMedian); 
%                 evi.relativeMin = evi.relativeMedian;
%                 evi.relativeMax = evi.relativeMedian;
% 
%                 ad.startOffset - ad.zeroOffset;
%                 labelInfo(counter).name = ad.startLabel;
%                 labelInfo(counter).time = ad.startOffset - ad.zeroOffset;
%                 labelInfo(counter).align = 'left';
%                 labelInfo(counter).info = ad.startInfo;
%                 labelInfo(counter).markData = ad.startMarkData;
%                 labelInfo(counter).fixed = true;
%                 counter = counter + 1;
%                 drewStartLabel = true;
%             end
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
       
       function [startRel, stopRel] = getStartStopRelativeToZeroByTrial(ad)
           ad.assertApplied();
           [start, stop, zero] = ad.getStartStopZeroByTrial();
           startRel = start - zero;
           stopRel = stop - zero;
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
            ad = ad.updateSummary();
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
            ad = ad.updateSummary();
        end
    end
    
    methods % Time-Aligning data transformations
        
        % use the alignment to shift the times in rawTimesCell to be zero relative
        % and filter by time window determined by getTimeInfo for each trial
        % if includePadding is true, will additional times found in the padWindow, see .setPadWindow
        function [alignedTimes, rawTimesMask] = getAlignedTimes(ad, rawTimesCell, includePadding)
            % filter the spikes within the window and recenter on zero
            if includePadding
                start = num2cell([ad.timeInfo.startPad]);
                stop = num2cell([ad.timeInfo.stopPad]);
            else
                start = num2cell([ad.timeInfo.start]);
                stop = num2cell([ad.timeInfo.stop]);
            end
            
            [alignedTimes, rawTimesMask] = cellfun(@fn, ...
                    makecol(rawTimesCell), ...
                    makecol(start), ...
                    makecol(stop), ...
                    makecol(num2cell([ad.timeInfo.zero])), ...
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

