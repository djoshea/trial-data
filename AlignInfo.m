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
        getEventTimesFn = @AlignInfo.getEventTimes;

        eventTimeRoundFn = @ceil;
        
        % struct array of nTrials x 1 containing the absolute times of each event
        % for the trials, as returned by getEventTimesFn
        eventInfo 

        % struct array of nTrials x 1 containing the absolute times of the
        % start, stop, zero events, as well as intervals and marks 
        timeInfo 

        % struct where each .eventName carries summary info about that event's relative
        % timing with respect to 0. Derives from timeInfo but speeds up plotting / labeling operations 
        summaryInfo

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

    methods % Binding to trial data
        function ad = AlignInfo(varargin)
            ad = ad@AlignDescriptor(varargin{:});
        end

        % bind this AlignInfo to a set of trials
        function ad = applyToTrialData(ad, R)
            ad.warnIfNoArgOut(nargout);
            ad.eventInfo = ad.requestEventInfo(R);
            [ad.timeInfo, ad.computedValid] = ad.buildTimeInfo();
            %ad.summary = ad.summarizeTimeInfo(ad.timeInfo);
            ad.applied = true;
            ad.manualInvalid = falsevec(ad.nTrials);
        end
        
        function nt = get.nTrials(ad)
            if ~ad.applied
                nt = 0;
            else
                nt = numel(ad.timeInfo);
            end
        end
        
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

        function times = getEventNthTimeVector(ad, event, n)
            if strcmp(n, 'end')
                fn = @(info) info.(event)(end);
            else
                fn = @(info) info.(event)(n);
            end
            times = arrayfun(fn, ad.eventInfo, ...
                'ErrorHandler', @(varargin) NaN);
            times = ad.eventTimeRoundFn(times);
        end
        
        function timeCell = getEventIndexedTimeCell(ad, event, n)
            % similar to above but returns cell array, and n may be be a
            % string of the form '1:2', '1:end', ':', etc
            if ~ischar(n) || strcmp(n, 'end')
                times = num2cell(ad.getEventNthTimeVector(event, n));
            else
                % must have a colon, parse into tokens
                pat = '(?<end1>end)?(?<ind1>-?\d*)?:(?<end2>end)?(?<ind2>-?\d*)?';
                str = 'end-1:end';
                info = regexp(n, pat, 'names', 'once');
                
                if isempty(info)
                    error('Unable to parse event index %s(%s)', event, n);
                end
                
                % convert ind1, ind2 to doubles
                if isempty(info.ind1)
                    ind1 = 0;
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
                timeCell = cellfun(@ad.eventTimeRoundFn, timeCell, 'UniformOutput', false);
            end
        end
        
        function timeCell = getEventIndexedTimeCellFillEmptyWithNaN(ad, event, n)
            timeCell  = ad.getEventIndexedTimeCell(event, n);
            emptyMask = cellfun(@isempty, timeCell);
            [timeCell{emptyMask}] = deal(NaN);
        end
        
        % get the aligned start/stop/zero/mark time windows for each trial, 
        % respecting all truncation and invalidation instructions
        function [timeInfo valid] = buildTimeInfo(ad, R)
            % returns a struct array with the actual time window and time of zero for trial i as
            %   timeInfo(i).start, .stop, .zero
            %
            % timeInfo(i).valid and valid(i) indicate whether trial i satisfied inclusion criteria
            % as specified by the various means of trial invalidation 
           
            padPre = ad.padPre;
            padPost = ad.padPost;

            eventNameFn = @(event) event;
            nTrials = numel(ad.eventInfo);

            t.valid = truevec(nTrials);
            t.invalidCause = cellvec(nTrials);

            t.trialStart = ad.getEventNthTimeVector('TrialStart', 1); 
            t.trialStop = ad.getEventNthTimeVector('TrialEnd', 1);

            % get start event
            t.start = ad.getEventNthTimeVector(ad.startEvent, ad.startEventIndex) + ad.startOffset;
            noStart = isnan(t.start);
            [t.invalidCause{noStart}] = deal(sprintf('Missing start event %s', ad.startUnabbreviatedLabel));
            t.valid(noStart & t.valid) = false;

            % get stop event
            t.stop = ad.getEventNthTimeVector(ad.stopEvent, ad.stopEventIndex) + ad.stopOffset;
            noStop = isnan(t.stop);
            [t.invalidCause{noStop & t.valid}] = deal(sprintf('Missing stop event %s', ad.stopUnabbreviatedLabel));
            t.valid(noStop) = false;

            % get zero alignment event
            t.zero= ad.getEventNthTimeVector(ad.zeroEvent, ad.zeroEventIndex) + ad.zeroOffset;
            noZero = isnan(t.zero);
            [t.invalidCause{noZero & t.valid}] = deal(sprintf('Missing zero event %s', ad.zeroUnabbreviatedLabel));
            t.valid(noZero) = false;
                        
            % get pad window
            t.startPad = t.start - padPre;
            t.stopPad = t.stop + padPost;
            
            % truncate trial end (including padding) based on truncateAfter events
            t.isTruncatedStop = falsevec(nTrials);
            for i = 1:length(ad.truncateAfterEvents)
                times = ad.getEventNthTimeVector(ad.truncateAfterEvents{i}, ad.truncateAfterEventsIndex{i}) + ad.truncateAfterOffsets(i);
                t.isTruncatedStop = t.isTruncatedStop | times < t.stopPad;
                t.stopPad = nanmin(t.stop, times);
                t.stop = t.stopPad - padPost;
            end
            
            % truncate trial start (including padding) based on truncateBefore events
            t.isTruncatedStart = falsevec(nTrials);
            for i = 1:length(ad.truncateBeforeEvents)
                times = ad.getEventNthTimeVector(ad.truncateBeforeEvents{i}, ad.truncateBeforeEventsIndex{i}) + ad.truncateBeforeOffsets(i);
                t.isTruncatedStart = t.isTruncatedStart | times > t.startPad;
                t.startPad = nanmax(t.startPad, times);
                t.start = t.startPad + padPre;
            end

            % mark trials as invalid if startPad:stopPad includes any invalidateEvents
            for i = 1:length(ad.invalidateEvents)
                timesCell = ad.getEventIndexedTimeVectorFillEmptyWithNan(ad.invalidateEvents{i}, ad.invalidOffset(i))
                maskInvalid = falsevec(length(t.startPad))
                for iT = 1:length(t.startPad)
                    maskInvalid(iT) = any(timesCell{iT} > t.startPad(iT) & timesCell{iT} < t.stopPad(iT));
                end
                t.valid(maskInvalid) = false;
                [t.invalidCause{maskInvalid}] = deal(sprintf('Overlapped invalidating event %s', ad.invalidUnabbreviatedLabels{i}));
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
            for iInt = 1:size(ad.intervalEvents, 1)
                startTimes = ad.getEventIndexedTimeVector(ad.intervalEventsStart{iInt}, ad.intervalEventsIndexStart, ad.intervalOffsetsStart(iInt));
                stopTimes = ad.getEventIndexedTimeVector(ad.intervalEventsStop{iInt}, ad.intervalEventsIndexStop, ad.intervalOffsetsStop(iInt));
                for iTrial = 1:nTrials
                    timeInfo(iTrial).intervalStart{iInt} = startTimes{iTrial};
                    timeInfo(iTrial).intervalStop{iInt} = stopTimes{iTrial};
                end 
            end
            
        end
    end

    methods % Labeling and axis drawing
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
            
            ti = ad.timeInfo;
            events = ad.getEventList(); 

            zeroTimes = [ad.timeInfo.(ad.zeroEvent)];

            for iEv = 1:length(events)
                event = events{iEv};

                times = ai.timeInfo;

                evi.fixed = true;
                evi.relativeMedian = ad.startOffset - ad.zeroOffset;
                evi.relativeList = repmat(ad.nTrials, 1, evi.relativeMedian); 
                evi.relativeMin = evi.relativeMedian;
                evi.relativeMax = evi.relativeMedian;

                ad.startOffset - ad.zeroOffset;
                labelInfo(counter).name = ad.startLabel;
                labelInfo(counter).time = ad.startOffset - ad.zeroOffset;
                labelInfo(counter).align = 'left';
                labelInfo(counter).info = ad.startInfo;
                labelInfo(counter).markData = ad.startMarkData;
                labelInfo(counter).fixed = true;
                counter = counter + 1;
                drewStartLabel = true;
            end
        end

        % struct with fields .time and and .name with where to label the time axis appropriately
        % pass along timeInfo so that the medians of non-fixed events can be labeled as well 
        function [labelInfo] = getLabelInfo(ad, varargin)
            % build a list of label names / times to mark on the time axis 
            % when using this alignment. Essentially, all events which are fixed relative
            % to zero (i.e. reference off the same event) will be included as labels
            % also, if 'tMin' and 'tMax' are specified, they will be drawn as well provided that 
            % the start and stop events aren't fixed relative to zero 
            
            timeInfo = ad.timeInfo; 
            
            % optionally provide time window which filters which labels will be included
            tMin = [];
            tMax = [];
            assignargs(varargin);

            labelInfo = struct('name', {}, 'time', {}, 'align', {}, 'info', {}, ...
                'markData', {}, 'fixed', {});
            counter = 1;

            drewStartLabel = false;

            % label the start event / min limit
            if ad.isStartFixedTime 
                if ~ad.isStartZero
                    % fixed but not redundant with zero
                    labelInfo(counter).name = ad.startLabel;
                    labelInfo(counter).time = ad.startOffset - ad.zeroOffset;
                    labelInfo(counter).align = 'left';
                    labelInfo(counter).info = ad.startInfo;
                    labelInfo(counter).markData = ad.startMarkData;
                    labelInfo(counter).fixed = true;
                    counter = counter + 1;
                    drewStartLabel = true;
                end
            else
                % label at minimum time
                labelInfo(counter).name = ad.startLabel;
                labelInfo(counter).time = min([timeInfo.start] - [timeInfo.zero]);
                labelInfo(counter).align = 'left';
                labelInfo(counter).info = ad.startInfo;
                labelInfo(counter).markData = ad.startMarkData;
                labelInfo(counter).fixed = true;
                counter = counter + 1;
                drewStartLabel = true;
            end
            
            % label the zero event provided that it lies within the start/stop window
            if ~ad.isZeroOutsideStartStop
                labelInfo(counter).name = ad.zeroLabel; 
                labelInfo(counter).time = 0;
                if drewStartLabel
                    labelInfo(counter).align = 'center';
                else
                    labelInfo(counter).align = 'left';
                end
                labelInfo(counter).info = ad.zeroInfo;
                labelInfo(counter).markData = ad.zeroMarkData;
                labelInfo(counter).fixed = true;
                counter = counter + 1;
            end

            % label the stop event / max limit
            if ad.isStopFixedTime
                if ~ad.isStopZero
                    % fixed but not redundant with zero
                    labelInfo(counter).name = ad.stopLabel; 
                    labelInfo(counter).time = ad.stopOffset - ad.zeroOffset;
                    labelInfo(counter).align = 'right';
                    labelInfo(counter).info = ad.stopInfo;
                    labelInfo(counter).markData = ad.stopMarkData;
                    labelInfo(counter).fixed = true;
                    counter = counter + 1;
                end
            else
                labelInfo(counter).name = ad.stopLabel; 
                labelInfo(counter).time = max([timeInfo.stop] - [timeInfo.zero]);
                labelInfo(counter).align = 'right';
                labelInfo(counter).info = ad.stopInfo;
                labelInfo(counter).markData = ad.stopMarkData;
                labelInfo(counter).fixed = true;
                counter = counter + 1;
            end

            % label each of the event marks that are fixed with respect to the zero event
            isMarkFixed = ad.isMarkFixedTime;
            for iMark = 1:length(ad.markEvents)
                if isMarkFixed(iMark)
                    % mark time is identical for each trial
                    labelInfo(counter).name = ad.markLabels{iMark};
                    labelInfo(counter).time = ad.markOffsets(iMark) - ad.zeroOffset;
                    labelInfo(counter).align = 'center';
                    labelInfo(counter).info = ad.markInfo{iMark};
                    labelInfo(counter).fixed = true;
                    labelInfo(counter).markData = ad.markMarkData(iMark);
                    counter = counter + 1;

                elseif ~isempty(timeInfo)
                    % compute the median time for this mark and put it down with <brackets>
                    % for now, mark only the first event occurrence
                    if iscell(timeInfo(1).mark)
                        markTimes = arrayfun(@(ti) ti.mark{iMark}(1) - ti.zero, ...
                            timeInfo([timeInfo.valid]), 'ErrorHandler', @(varargin) NaN);
                    else
                        markTimes = arrayfun(@(ti) ti.mark(iMark) - ti.zero, ...
                            timeInfo([timeInfo.valid]));
                    end
                    
                    medianMarkTime = nanmedian(markTimes);
                    minMarkTime = nanmin(markTimes);
                    maxMarkTime = nanmax(markTimes);

                    if abs(minMarkTime - medianMarkTime) <= ad.markRelativeDeltaIgnore && ...
                       abs(maxMarkTime - medianMarkTime) <= ad.markRelativeDeltaIgnore
                       % range acceptable, don't use < > brackets
                       fmat = '%s';
                       labelInfo(counter).fixed = true;
                    else
                       % mark median with < > brackets to indicate there is
                       % a range
                        if ad.markPlotMedians
                            fmat = '<%s>';
                            labelInfo(counter).fixed = false;
                        else
                            % not plotting medians, move to next mark
                            continue;
                        end
                    end
                    
                    labelInfo(counter).name = sprintf(fmat, ad.markLabels{iMark});
                    labelInfo(counter).time = medianMarkTime;
                    labelInfo(counter).align = 'center';
                    labelInfo(counter).markData = ad.markMarkData(iMark);
                    labelInfo(counter).info = ad.markInfo{iMark};
                    counter = counter + 1;
                end
            end
                 
            if ~isempty(tMin)
                % start not fixed, include one for the lower limit tMin
                % first check whether there is an existing label (from a
                % mark) at this point already...
                if ~any(floor([labelInfo.time]) == floor(tMin))
                    labelInfo(counter).name = ad.buildLabel(ad.zeroEvent, ad.zeroEventIndex, tMin);
                    labelInfo(counter).time = tMin;
                    labelInfo(counter).align = 'left';
                    labelInfo(counter).markData = false; % this is just convenience, don't mark on data
                    labelInfo(counter).fixed = true;
                    counter = counter + 1;
                    drewStartLabel = true;
                end
            end
                     
            if ~isempty(tMax)
                % stop not fixed, include one for the upper limit tMin
                % first check whether there is an existing label (from a
                % mark) at this point already...
                if ~any(floor([labelInfo.time]) == floor(tMax))
                    labelInfo(counter).name = ad.buildLabel(ad.zeroEvent, ad.zeroEventIndex, tMax);
                    labelInfo(counter).time = tMax;
                    labelInfo(counter).align = 'right';
                    labelInfo(counter).markData = false;
                    labelInfo(counter).fixed = true;
                    counter = counter + 1;
                end
            end

            % generate default label info where missing
            cmap = jet(length(labelInfo));
            for i = 1:length(labelInfo)
                default = struct('color', cmap(i, :), 'size', 10, 'marker', 'o');
                if isempty(labelInfo(i).info)
                    labelInfo(i).info = default;
                else
                    labelInfo(i).info = structMerge(default, labelInfo(i).info, 'warnOnOverwrite', false);
                end
            end
            
            times = [labelInfo.time];
            timeMask = true(size(times));
            if ~isempty(tMin)
                timeMask = timeMask & times >= tMin;
            end
            if ~isempty(tMax)
                timeMask = timeMask & times <= tMax;
            end
            labelInfo = labelInfo(timeMask);
            
            labelInfo = makecol(labelInfo);
        end
        
        function [info valid] = getIntervalInfoByCondition(ad, ciOrig, varargin)
            % using the conditions picked out by ConditionInfo ci and the event info found in data
            % R, compute the start and stop times for each interval defined by this AlignDescriptor
            % 
            % info is a nConditions x 1 cell containing a struct with the data on each
            % of the intervals for that condition
            %
            % currently uses first trial from each condition's intervals,
            % but double checks that all trials have the same number of
            % start/stops within each interval
        
            timeInfo = ad.timeInfo;
            
            % TODO Implement condition matching
            valid = [timeInfo.valid];
            ci = ciOrig.copy();
            ci.markInvalid(~valid);
            
            nIntervals = size(ad.intervalEvents,1);
            info = struct();
             
            for iC = 1:ci.nConditions
                rMask = ci.listByCondition{iC};
                info(iC).interval = cell(nIntervals, 1);
                
                for iInt = 1:nIntervals
                    if isempty(rMask)
                        % no trials, fill with nan
                        info(iC).interval{iInt} = [NaN NaN];
                    else
                        % check that all trials within condition have the same number
                        % of periods for this interval
                        nPeriods = arrayfun(@(ti) size(ti.interval{iInt}, 1), timeInfo(rMask));
                        if length(unique(nPeriods)) > 1
                            debug('WARNING: Trials within condition have differing number of periods for interval %d\n', iInt);
                        end

                        % grab the interval info from the first trial
                        info(iC).interval{iInt} = timeInfo(rMask(1)).interval{iInt} - timeInfo(rMask(1)).zero;
                    end
                end
            end
        end

        % use the alignment to shift the times in rawTimesCell to be zero relative
        % and filter by time window determined by getTimeInfo for each trial
        % INCLUDES additional times found in the padWindow, see .setPadWindow
        function [alignedTimes rawTimesMask] = getAlignedTimes(ad, rawTimesCell)
            timeInfo = ad.timeInfo;
            
            % filter the spikes within the window and recenter on zero
            [alignedTimes, rawTimesMask] = cellfun(@fn, ...
                    makecol(rawTimesCell), ...
                    makecol(num2cell([timeInfo.startPad])), ...
                    makecol(num2cell([timeInfo.stopPad])), ...
                    makecol(num2cell([timeInfo.zero])), ...
                    'UniformOutput', false);
                
            alignedTimes(~ad.valid) = {[]};
            
            function [alignedTimes, mask] = fn(rawTimes, tStart, tEnd, tZero)
                mask = rawTimes >= tStart & rawTimes <= tEnd;
                alignedTimes = rawTimes(mask) - tZero;
            end
        end
        
        function [alignedData alignedTime] = getAlignedTimeseries(ad, dataCell, timeCell, varargin)
            [alignedTime rawTimesMask] = ad.getAlignedTimes(timeCell, varargin{:});
            %dataCell = cellfun(@makecol, dataCell, 'UniformOutput', false);
            alignedData = cellfun(@(data, mask) data(mask, :), dataCell, rawTimesMask, ...
                'UniformOutput', false, 'ErrorHandler', @(varargin) []);
        end

        % used to annotate a time axis with the relevant start/stop/zero/marks
        % non-fixed marks as <markLabel> unless the range is less than a specified 
        % noise-threshold, in which case it is marked as though it were fixed
        function drawTimeAxis(ad, varargin)
            timeInfo = ad.timeInfo;

            % uses ad.labelInfo to call drawPrettyAxis
            tLims = [];
            xLabel = ''; 
            axh = [];
            drawY = true; % also draw the y axis while we're here? otherwise they'll be nothing there
            setXLim = false;
            assignargs(varargin);

            if isempty(axh)
                axh = gca;
            end
            if isempty(tLims)
                if setXLim
                    tLims = ad.getTimeAxisLims(timeInfo);
                else
                    tLims = xlim(axh);
                end
            end
            tMin = tLims(1);
            tMax = tLims(2);
              
            labelInfo = ad.getLabelInfo(timeInfo, 'tMin', tMin, 'tMax', tMax);
            tickPos = [labelInfo.time];
            tickLabels = {labelInfo.name};
            tickAlignments = {labelInfo.align};
                     
            if setXLim
                xlim([min(tickPos), max(tickPos)]);
            end
            
            if drawY
                makePrettyAxis('yOnly', true); 
            else
                axis(axh, 'off');
                box(axh, 'off');
            end
            
            if all(~isnan(tLims))
                drawAxis(tickPos, 'tickLabels', tickLabels, 'tickAlignments', tickAlignments, 'axh', axh); 
            end
        end

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
        function drawOnData(ad, timeData, data, varargin)
            p = inputParser();
            p.addParamValue('drawLegend', false, @islogical);
            p.parse(varargin{:});
            
            timeInfo = ad.timeInfo;

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

    methods(Static)
        function [timeData] = getEventTimes(R, eventNameList)
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
       
        function [startMs stopMs] = getTrialLengths(R)
            nTrials = AlignInfo.getNTrialsFromData(R);
            if isstruct(R)
                
                startMs = zeros(nTrials, 1);
                stopMs = [R.length];

                assert(numel(stopMs) == nTrials, 'R.length does not contain a value in all trials');

            elseif isa(R, 'TrialData')
                startMs = zeros(nTrials, 1);
                stopMs = R.getParam('duration');

            else
                error('Unsupported trial data type, please specify .getEventTimesFn');
            end
        end

    end
end

