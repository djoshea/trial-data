classdef AlignDescriptor 
    % Describes a particular way to slice and align trial data in time based off 
    % of event times in each trial. Typically analog data and neural data are 
    % aligned relative to the times of certain events before further analysis takes
    % place. This class provides a common language for expressing these different
    % types of aligning data.

    properties
        omitRedundantLabel = false; % e.g. when zero=event1, start=event1-20, label start as "-20" rather than "event1-20"

        % abbreviate event names to capital/lowercase initials automatically
        % override this by calling addEventAbbrev
        autoAbbreviateLabels = false;
        
        % use capitalized initials instead of lowercased initials
        autoAbbreviateMakeUpper = true;
        
        % for description
        name % set manually, or will be auto-populated
        
        % minimum duration to mark as valid
        minDuration = 0;
         
    end
        
    properties(SetAccess=protected, Hidden)
        nameDefault
        
        % how to handle the window falling outside of the trial
        outsideOfTrialMode = AlignDescriptor.TRUNCATE; % TRUNCATE or INVALIDATE or IGNORE
        
        % allow an additional time window around the alignment interval
        % see .setPadWindow below
        padPre = 0;
        padPost = 0;
    end
    
    properties(Access=protected)
        % eventAbbrevLookup.eventName == abbreviatedName, overrides auto abbreviation
        % use ad = ad.abbrev(eventName, abbrev) to add an entry
        eventAbbrevLookup
    end

    properties(SetAccess=protected) % SetAccess=protected
        % slice times >= t_startEvent + startOffset
        startEvent;
        startEventIndex = 1;
        startOffset = 0;
        startAppear
        startMark % mark this point on a data trace in drawOnData
        startDefault = true;
        startLabelStored = '';

        % slice times <= t_stopEvent + stopOffset
        stopEvent;
        stopOffset = 0;
        stopEventIndex = 1;
        stopAppear
        stopMark 
        stopDefault = true;
        stopLabelStored = '';

        % shift time such that t_zeroEvent + zeroOffset is t=0
        zeroEvent;
        zeroOffset = 0;
        zeroEventIndex = 1;
        zeroAppear
        zeroMark
        zeroDefault = true;
        zeroLabelStored = '';

        % move interval start forward in time to avoid these times
        truncateBeforeEvents = {};
        truncateBeforeEventsIndex = {}
        truncateBeforeOffsets

        % move interval stop backward in time to avoid these times
        truncateAfterEvents = {};
        truncateAfterEventsIndex = {};
        truncateAfterOffsets 

        % invalidate any trial whose window includes these events
        invalidateEvents = {};
        invalidateEventsIndex = {};
        invalidateOffsets

        % purely for annotation, used by getTimeInfo and drawTimeAxis methods 
        markEvents = {};
        markEventsIndex
        markOffsets 
        markLabelsStored = {};
        markAppear
        markShowOnData % whether to mark on data traces
        markShowOnAxis % whether to mark on the axis

        % for marking time intervals on the time axis as colored rectangles
        intervalEventsStart = {}; % n x 1 cell array of start/stop events
        intervalEventsStop = {}; % n x 1 cell array of start/stop events
        intervalEventsIndexStart % n x 1 array of event indices
        intervalEventsIndexStop % n x 1 array of event indices
        intervalOffsetsStart % n x 1 array of offsets
        intervalOffsetsStop % n x 1 array of offsets
        intervalLabelsStored = {};

        % 'first' means use events from first trial in each condition
        %intevalMultiTrialMode = 'first';
        intervalConditionMatch = {}; % n x 1 cell array of structs with .attrName = attrValue(s)
        intervalAppear % miscellaneous
        
        intervalShowOnData % whether to mark on data traces
        intervalShowOnAxis % whether to mark on the axis
        
        % Timestamp rounding. If active,
        % timeseries and timestamps will be shifted such that the
        % difference between each timestamp and zero is an exact integer multiple of
        % minTimeDelta. The alignment will be performed to maintain single
        % timestamp alignments (i.e. start == stop) as a single time point
        roundTimes = true; % boolean indicating whether resampling is active
        minTimeDelta = 1; % the minimum acceptable spacing between timestamps to which alignment will be enforced, relative to the zero event
    end

    properties(Constant, Hidden)
        TRUNCATE = 'truncate';
        INVALIDATE = 'invalidate';
        IGNORE = 'ignore';
        AUTO = '<auto>';
    end

    properties(Dependent, Hidden)
        padWindow
        markLabels
        intervalLabels
        
        startLabel
        stopLabel
        zeroLabel

        startUnabbreviatedLabel
        stopUnabbreviatedLabel
        zeroUnabbreviatedLabel
        markUnabbreviatedLabels
        intervalUnabbreviatedLabels
        truncateBeforeUnabbreviatedLabels
        truncateAfterUnabbreviatedLabels
        invalidateUnabbreviatedLabels
    end

    properties(Dependent)
        isFixedLength % is the time window guaranteed fixed length
        isStartFixedTime % is the start event always at a fixed time relative to zero
        isStopFixedTime % is the start event always at a fixed time relative to zero
        isMarkFixedTime % is each mark event always at a fixed time relative to zero
        isIntervalFixedTime % n x 2 is interval start/end fixed time relative to zero 

        isStartZero % is the start event always == 0
        isStopZero % is the stop event always == 0
        isMarkZero % is each mark event always == 0
        
        isStartStopEqual % is start event same as start event

        isZeroOutsideStartStop % is the zero event guaranteed to lie outside the start/stop window
        
        isFullTrial
        
        nMarks % number of marked events (length .markEvents)
        nIntervals % number of marked intervals (length .markIntervals)

    end

    methods % Dependent properties
        function tf = get.isFixedLength(ad)
            % the window is always the same length if the start / stop events are the same
            if strcmp(ad.outsideOfTrialMode, ad.INVALIDATE) ||...
                strcmp(ad.outsideOfTrialMode, ad.IGNORE) && ...
                isempty(ad.truncateBeforeEvents) && isempty(ad.truncateAfterEvents)

                tf = strcmp(ad.startEvent, ad.stopEvent) && ad.startEventIndex == ad.stopEventIndex;
            else
                tf = false;
            end
        end

        function tf = get.isStartFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.startEvent) && isequal(ad.zeroEventIndex, ad.startEventIndex);
        end

        function tf = get.isStopFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.stopEvent) && isequal(ad.zeroEventIndex, ad.stopEventIndex);
        end

        function tf = get.isMarkFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.markEvents) & isequal(ad.zeroEventIndex, ad.markEventsIndex);
        end

        function tf = get.isIntervalFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.intervalEventsStart) & strcmp(ad.zeroEvent, ad.intervalEventsStop) & ...
                isequal(ad.zeroEventIndex, ad.intervalEventsIndexStart) & ...
                isequal(ad.zeroEventIndex, ad.intervalEventsIndexStop);
        end

        function tf = get.isStartZero(ad)
            tf = strcmp(ad.zeroEvent, ad.startEvent) && isequal(ad.zeroEventIndex, ad.startEventIndex) && ...
                ad.zeroOffset == ad.startOffset;
        end

        function tf = get.isStopZero(ad)
            tf = strcmp(ad.zeroEvent, ad.stopEvent) && isequal(ad.zeroEventIndex, ad.stopEventIndex) && ...
                ad.zeroOffset == ad.stopOffset;
        end

        function tf = get.isMarkZero(ad)
            if isempty(ad.nMarks)
                tf = false;
            else
                tf = strcmp(ad.zeroEvent, ad.markEvents) & cellfun(@(x) strcmp(x, ':') || isequal(x, ad.zeroEventIndex), ad.markEventsIndex) & ...
                    ad.markOffsets == ad.zeroOffset;
            end
        end
        
        function tf = get.isStartStopEqual(ad)
            tf = strcmp(ad.startEvent, ad.stopEvent) && isequal(ad.startEventIndex, ad.stopEventIndex) && ...
                ad.startOffset == ad.stopOffset;
        end

        function tf = get.isZeroOutsideStartStop(ad)
            tf = (ad.isStartFixedTime && ad.startOffset > 0) || (ad.isStopFixedTime && ad.stopOffset < 0);
        end
        
        function tf = get.isFullTrial(ad)
            tf = strcmp(ad.startEvent, 'TrialStart') && ad.startOffset == 0 && ...
                strcmp(ad.stopEvent, 'TrialEnd') && ad.stopOffset == 0;
        end
        
        function n = get.nMarks(ad)
            n = numel(ad.markEvents);
        end
        
        function n = get.nIntervals(ad)
            n = numel(ad.intervalEventsStart);
        end
        
        function name = get.name(ad)
            if isempty(ad.name)
                name = ad.getStartStopZeroDescription();
            else
                name = ad.name;
            end
        end
        
        function ad = setName(ad, name)
            ad.name = name;
            ad.nameDefault = false;
        end

        function padWindow = get.padWindow(ad)
            padWindow = [ad.padPre ad.padPost];
        end
    end
    
    methods % Auto-generated / manually set abbreviated labels
        function str = get.startLabel(ad)
            if isempty(ad.startLabelStored)
                if ad.isStartFixedTime && ~ad.isStartZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', [], ad.startOffset);
                else
                    str = ad.buildLabel(ad.startEvent, ad.startEventIndex, ad.startOffset);
                end
            else
                str = ad.startLabelStored;
            end
        end
        
        function ad = set.startLabel(ad, v)
            ad.startLabelStored = v;
        end
        
        function str = get.stopLabel(ad)
            if isempty(ad.stopLabelStored)
                if ad.isStopFixedTime && ~ad.isStopZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', [], ad.stopOffset);
                else
                    str = ad.buildLabel(ad.stopEvent, ad.stopEventIndex, ad.stopOffset);
                end
            else
                str = ad.stopLabelStored;
            end
        end
        
        function ad = set.stopLabel(ad, v)
            ad.stopLabelStored = v;
        end

        function str = get.zeroLabel(ad)
            if isempty(ad.zeroLabelStored)
                str = ad.buildLabel(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset);
            else
                str = ad.zeroLabelStored;
            end
        end
        
        function ad = set.zeroLabel(ad, v)
            ad.zeroLabelStored = v;
        end
            
        function markLabels = get.markLabels(ad)
            if isempty(ad.markLabelsStored)
                ad.markLabelsStored = cell(length(ad.markEvents), 1);
            end
            
            markLabels = cell(length(ad.markLabelsStored), 1);
            for iMark = 1:length(ad.markLabelsStored)
                if ~strcmp(ad.markLabelsStored{iMark}, AlignDescriptor.AUTO)
                    % manual label
                    markLabels{iMark} = ad.markLabelsStored{iMark};
                elseif ad.isMarkFixedTime(iMark) && ad.omitRedundantLabel
                    % same as zero, just include the offset (e.g. '+100')
                    markLabels{iMark} = ad.buildLabel('', [], ad.markOffsets(iMark));
                else
                    % full label string, omit the index since everything
                    % will be labeled the same
                    markLabels{iMark} = ad.buildLabel(ad.markEvents{iMark}, ...
                        [], ad.markOffsets(iMark));
                end
            end
        end
        
        function ad = set.markLabels(ad, v)
            ad.markLabelsStored = v;
        end

        function intervalLabels = get.intervalLabels(ad)
            nIntervals = size(ad.intervalEventsStart, 1); %#ok<*PROP>
            if isempty(ad.intervalLabelsStored)
                ad.intervalLabelsStored = cell(nIntervals, 1);
            end
            
            intervalLabels = cell(length(ad.intervalLabelsStored), 1);
            for i = 1:nIntervals
                if ~strcmp(ad.intervalLabelsStored{i}, AlignDescriptor.AUTO)
                    % manual label
                    intervalLabels{i} = ad.intervalLabelsStored{i};
                else
                    if ad.isIntervalFixedTime(i) && ad.omitRedundantLabel
                        % same as zero, just include the offset (e.g. '+100')
                        label1 = ad.buildLabel('', [], ad.intervalOffsetsStart(i));
                    else
                        label1 = ad.buildLabel(ad.intervalEventsStart{i}, ...
                            [], ad.intervalOffsetsStart(i));
                    end
                    if ad.isIntervalFixedTime(i) && ad.omitRedundantLabel
                        % same as zero, just include the offset (e.g. '+100')
                        label2 = ad.buildLabel('', [], ad.intervalOffsetsStop(i));
                    else
                        label2 = ad.buildLabel(ad.intervalEventsStop{i}, ...
                            [], ad.intervalOffsetsStop(i));
                    end

                    intervalLabels{i} = sprintf('%s : %s', label1, label2);
                end
            end
        end
        
        function ad = set.intervalLabels(ad, v)
            ad.intervalLabelsStored = v;
        end
    end

    methods % Unabbreviated raw labels
        function str = get.startUnabbreviatedLabel(ad)
            str = ad.buildUnabbreviatedLabel(ad.startEvent, ad.startEventIndex, ad.startOffset);
        end

        function str = get.stopUnabbreviatedLabel(ad)
            str = ad.buildUnabbreviatedLabel(ad.stopEvent, ad.stopEventIndex, ad.stopOffset);
        end

        function str = get.zeroUnabbreviatedLabel(ad)
            str = ad.buildUnabbreviatedLabel(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset);
        end

        function strCell = get.markUnabbreviatedLabels(ad)
            strCell = cellfun(@ad.buildUnabbreviatedLabel, ad.markEvents, ...
                ad.markEventsIndex, num2cell(ad.markOffsets), 'UniformOutput', false);
        end

        function strCell = get.intervalUnabbreviatedLabels(ad)
            % TODO fix this
            strCellStart = cellfun(@ad.buildUnabbreviatedLabel, ad.intervalEventsStart, ...
                ad.intervalEventsIndexStart, num2cell(ad.intervalOffsetsStart), 'UniformOutput', false);
            strCellStop = cellfun(@ad.buildUnabbreviatedLabel, ad.intervalEventsStop, ...
                ad.intervalEventsIndexStop, num2cell(ad.intervalOffsetsStop), 'UniformOutput', false);
            
            strCell = cellfun(@(s1, s2) sprintf('%s : %s', s1, s2), strCellStart, strCellStop, ...
                'UniformOutput', false);
        end

        function strCell = get.invalidateUnabbreviatedLabels(ad)
            strCell = cellfun(@ad.buildUnabbreviatedLabel, ad.invalidateEvents, ...
                 ad.invalidateEventsIndex, num2cell(ad.invalidateOffsets), 'UniformOutput', false);
        end

        function strCell = get.truncateBeforeUnabbreviatedLabels(ad)
            strCell = cellfun(@ad.buildUnabbreviatedLabel, ad.truncateBeforeEvents, ...
                ad.truncateBeforeEventsIndex, num2cell(ad.truncateBeforeOffsets), 'UniformOutput', false);
        end

        function strCell = get.truncateAfterUnabbreviatedLabels(ad)
            strCell = cellfun(@ad.buildUnabbreviatedLabel, ad.truncateAfterEvents, ...
                ad.truncateAfterEventsIndex, num2cell(ad.truncateAfterOffsets), 'UniformOutput', false);
        end

    end

    methods % constructor, manual event specification
        function ad = AlignDescriptor(varargin)
            ad.eventAbbrevLookup = struct();

            ad.startEvent = 'TrialStart';
            ad.stopEvent = 'TrialEnd';
            ad.zeroEvent = 'TimeZero';

            % either use named property value pairs or string syntax
            if nargin == 1 && ischar(varargin{1});
                ad = ad.parseDescriptionString(varargin{1});
            else
                ad = structargs(ad, varargin);
            end
        end

        function ad = update(ad)
            ad.warnIfNoArgOut(nargout);
            % does nothing here, used primarily in AlignInfo
        end
        
        function ad = postUpdateMark(ad)
            ad.warnIfNoArgOut(nargout);
            % does nothing here, used primarily in AlignInfo
        end
        
        function ad = postUpdateInterval(ad)
            ad.warnIfNoArgOut(nargout);
            % does nothing here, used primarily in AlignInfo
        end

        % return a list of event names this AlignDescriptor references, plus 'start' and 'end'
        function eventList = getEventList(ad)
            eventList = unique([{'TrialStart'; 'TrialEnd'}; ...
                ad.startEvent; ad.stopEvent; ad.zeroEvent; ...
                ad.truncateBeforeEvents; ad.truncateAfterEvents; ad.invalidateEvents; ...
                ad.markEvents; ad.intervalEventsStart; ad.intervalEventsStop]);
            eventList = setdiff(eventList, 'TimeZero');
        end
        
        function assertHasEvent(ad, eventName) %#ok<INUSD>
            % this does nothing here, but is overriden in AlignInfo to
            % prevent post-hoc modifications that refer to non-existent
            % events
        end

        function ad = start(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParameter('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParameter('as', '', @ischar);
            p.addParameter('mark', false, @islogical);
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.assertHasEvent(eventName);
            ad.startEvent = eventName;
            ad.startEventIndex = p.Results.index;
            ad.startOffset = offset;
            ad.startMark = p.Results.mark;
            
            if ~isempty(p.Results.as)
                ad.startLabel = p.Results.as;
            else
                % use default again
                ad.startLabel = '';
            end
            appear = p.Results.appear;
            if ~isempty(p.Results.color)
                appear.Color = p.Results.color;
            end
            ad.startAppear = appear;
            ad.startDefault = false;

            % if no zero is specified, determine it from the start event
            ad = ad.setDefaultZero();

            ad = ad.update();
        end

        function ad = stop(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParameter('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParameter('as', '', @ischar);
            p.addParameter('mark', false, @islogical);
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.assertHasEvent(eventName);
            ad.stopEvent = eventName;
            ad.stopEventIndex = p.Results.index;
            ad.stopOffset = offset;
            ad.stopMark = p.Results.mark;

            if ~isempty(p.Results.as)
                ad.stopLabel = p.Results.as;
            else
                % use default again
                ad.stopLabel = '';
            end
            appear = p.Results.appear;
            if ~isempty(p.Results.color)
                appear.Color = p.Results.color;
            end
            ad.stopAppear = appear;
            ad.stopDefault = true;

            ad = ad.update();
        end

        function ad = pad(ad, pre, post, varargin)
            ad.warnIfNoArgOut(nargout);
            if nargin < 3
                if isscalar(pre)
                    post = pre;
                else
                    post = pre(2);
                    pre = pre(1);
                end
            end
               
            ad.padPre = pre;
            ad.padPost = post;

            ad = ad.update();
        end

        function ad = zero(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParameter('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParameter('as', '', @ischar);
            %p.addParameter('mark', true, @islogical);
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.assertHasEvent(eventName);
            ad.zeroEvent = eventName;
            ad.zeroEventIndex = p.Results.index;
            ad.zeroOffset = offset;
            %ad.zeroMark = p.Results.mark;

            if ~isempty(p.Results.as)
                ad.zeroLabel = p.Results.as;
            else
                % use default again
                ad.zeroLabel = '';
            end
            appear = p.Results.appear;
            if ~isempty(p.Results.color)
                appear.Color = p.Results.color;
            end
            ad.zeroAppear = p.Results.appear;
            ad.zeroDefault = true;

            ad = ad.update();
        end
        
        function idx = findInterval(ad, eventStart, indexStart, offsetStart, ...
                eventStop, indexStop, offsetStop)
            if ad.nIntervals == 0
                idx = [];
                return;
            end
            
            match = truevec(ad.nIntervals);
            
            match = match & makecol(strcmp(ad.intervalEventsStart, eventStart));
            match = match & makecol(cellfun(@(x) isequal(x, indexStart), ad.intervalEventsIndexStart));
            match = match & makecol(ad.intervalOffsetsStart == offsetStart);
            match = match & makecol(strcmp(ad.intervalEventsStop, eventStop));
            match = match & makecol(cellfun(@(x) isequal(x, indexStop), ad.intervalEventsIndexStop));
            match = match & makecol(ad.intervalOffsetsStop == offsetStop);
            idx = find(match);
        end

        function ad = interval(ad, eventStart, eventStop, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addParameter('offsetStart', 0, @isscalar);
            p.addParameter('offsetStop', 0, @isscalar);
            p.addParameter('indexStart', ':', @(x) ischar(x) || isscalar(x));
            p.addParameter('indexStop', ':', @(x) ischar(x) || isscalar(x));
            p.addParameter('as', AlignDescriptor.AUTO, @ischar);
            p.addParameter('appear', AppearanceSpec(), @(x) isa(x, 'AppearanceSpec'));
            p.addParameter('color', [], @(x) true);
            p.addParameter('showOnData', true, @islogical);
            p.addParameter('showOnAxis', true, @islogical);
            %p.addParameter('conditionMatch', struct(), @(x) isstruct(x) && isscalar(x));
            p.parse(varargin{:});
            %conditionMatch = p.Results.conditionMatch;
            as = p.Results.as;

            % try inferring offsets, indexes from provided strings, then
            % overwrite with provided params
            [eventNameStart, indexStart, offsetStart] = ad.parseEventOffsetString(eventStart, 'interval start', 'defaultIndex', ':');
            [eventNameStop, indexStop, offsetStop] = ad.parseEventOffsetString(eventStop, 'interval stop', 'defaultIndex', ':');
            if ~ismember('offsetStart', p.UsingDefaults), offsetStart = p.Results.offsetStart; end
            if ~ismember('offsetStop', p.UsingDefaults), offsetStop = p.Results.offsetStop; end
            if ~ismember('indexStart', p.UsingDefaults), indexStart = p.Results.indexStart; end
            if ~ismember('indexStop', p.UsingDefaults), indexStop = p.Results.indexStop; end
            ad.assertHasEvent(eventNameStart);
            ad.assertHasEvent(eventNameStop);
            
            idx = ad.findInterval(eventNameStart, indexStart, offsetStart, ...
                eventNameStop, indexStop, offsetStop);
            if ~isempty(idx)
                warning('Replacing existing matching interval');
                iInterval = idx(1);
            else
                iInterval = size(ad.intervalEventsStart,1)+1;
            end
            
            ad.intervalEventsStart{iInterval} = eventNameStart; 
            ad.intervalEventsStop{iInterval} = eventNameStop; 
            ad.intervalEventsIndexStart{iInterval} = indexStart;
            ad.intervalEventsIndexStop{iInterval} = indexStop;
            ad.intervalOffsetsStart(iInterval) = offsetStart;
            ad.intervalOffsetsStop(iInterval) = offsetStop;
            ad.intervalLabelsStored{iInterval} = as;
            
            appear = p.Results.appear;
            if ~isempty(p.Results.color)
                appear.Color = p.Results.color;
            end
            ad.intervalAppear{iInterval} = appear;
            
            ad.intervalEventsStart = makecol(ad.intervalEventsStart);
            ad.intervalEventsStop = makecol(ad.intervalEventsStop);
            ad.intervalEventsIndexStart = makecol(ad.intervalEventsIndexStart);
            ad.intervalEventsIndexStop = makecol(ad.intervalEventsIndexStop);
            ad.intervalOffsetsStart = makecol(ad.intervalOffsetsStart);
            ad.intervalOffsetsStop = makecol(ad.intervalOffsetsStop);
            
            ad.intervalAppear = makecol(ad.intervalAppear);
            ad.intervalLabelsStored = makecol(ad.intervalLabelsStored);
            ad.intervalShowOnData(iInterval) = p.Results.showOnData;
            ad.intervalShowOnData = makecol(ad.intervalShowOnData);
            ad.intervalShowOnAxis(iInterval) = p.Results.showOnAxis;
            ad.intervalShowOnAxis = makecol(ad.intervalShowOnAxis);

            ad = ad.postUpdateInterval();
        end
        
        function ad = removeIntervalsByIdx(ad, mask)
            ad.warnIfNoArgOut(nargout);

            mask = makecol(TensorUtils.vectorIndicesToMask(mask, ad.nIntervals));
            
            ad.intervalEventsStart = ad.intervalEventsStart(~mask);
            ad.intervalEventsStop = ad.intervalEventsStop(~mask);
            ad.intervalEventsIndexStart = ad.intervalEventsIndexStart(~mask);
            ad.intervalEventsIndexStop = ad.intervalEventsIndexStop(~mask);
            ad.intervalOffsetsStart = ad.intervalOffsetsStart(~mask);
            ad.intervalOffsetsStop = ad.intervalOffsetsStop(~mask);
            ad.intervalAppear = ad.intervalAppear(~mask);
            ad.intervalLabelsStored = ad.intervalLabelsStored(~mask);
            ad.intervalShowOnData = ad.intervalShowOnData(~mask);
            ad.intervalShowOnAxis = ad.intervalShowOnAxis(~mask);
            
            ad = ad.postUpdateInterval();
        end
        
        function ad = removeInterval(ad, eventStart, indexStart, offsetStart, ...
                eventStop, indexStop, offsetStop)
            ad.warnIfNoArgOut(nargout);

            idx = ad.findInterval(eventStart, indexStart, offsetStart, ...
                eventStop, indexStop, offsetStop);
            ad = ad.removeIntervalsByIdx(idx);
        end
        
        
        function idx = findMarkByString(ad, str)
            [eventName, index, offset] = ad.parseEventOffsetString(str, ...
                'mark', 'defaultIndex', ':');
            idx = ad.findMark(eventName, index, offset);
        end
        
        function idx = findMark(ad, eventName, varargin)
            p = inputParser();
            p.addOptional('index', ':', @(x) isempty(x) || ischar(x) || isscalar(x));
            p.addOptional('offset', 0, @isscalar);
            p.parse(varargin{:});
            
            if ad.nMarks == 0
                idx = [];
                return;
            end
            
            match = truevec(ad.nMarks);
            
            match = match & makecol(strcmp(ad.markEvents, eventName));
            match = match & makecol(cellfun(@(x) isequal(x, p.Results.index), ad.markEventsIndex));
            match = match & makecol(ad.markOffsets == p.Results.offset);
            idx = find(match);
        end
        
        function ad = mark(ad, eventStr, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParameter('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.addParameter('as', AlignDescriptor.AUTO, @ischar);
            p.addParameter('color', [], @(x) isempty(x) || ischar(x) || isvector(x));
            p.addParameter('appear', [], @(x) isempty(x) || isa(x, 'AppearanceSpec'));
            p.addParameter('showOnData', true, @islogical);
            p.addParameter('showOnAxis', true, @islogical);
            p.parse(varargin{:});
            
            [eventName, index, offset] = ad.parseEventOffsetString(eventStr, ...
                'mark', 'defaultIndex', ':');
            if ~ismember('offset', p.UsingDefaults), offset = p.Results.offset; end
            if ~ismember('index', p.UsingDefaults), index = p.Results.index; end

            ad.assertHasEvent(eventName);

            % check for existing mark which matches
            idx = ad.findMark(eventName, index, offset);
            
            if ~isempty(idx)
                warning('Replacing existing mark');
                iMark = idx(1);
            else
                iMark = length(ad.markEvents)+1;
            end
            
            ad.markEvents{iMark} = eventName;
            ad.markEvents = makecol(ad.markEvents);
            ad.markEventsIndex{iMark} = index;
            ad.markEventsIndex = makecol(ad.markEventsIndex);
            ad.markOffsets(iMark) = offset;
            ad.markOffsets = makecol(ad.markOffsets);
            
            appear = p.Results.appear;
            if isempty(appear)
                appear = AppearanceSpec();
                appear.Color = ad.getNextMarkColor();
            end
            if ~isempty(p.Results.color)
                appear.Color = p.Results.color;
            end
            ad.markAppear{iMark} = appear;
            ad.markAppear = makecol(ad.markAppear);
            ad.markShowOnData(iMark) = p.Results.showOnData;
            ad.markShowOnData = makecol(ad.markShowOnData);
            ad.markShowOnAxis(iMark) = p.Results.showOnAxis;
            ad.markShowOnAxis = makecol(ad.markShowOnAxis);
            ad.markLabelsStored{iMark,1} = p.Results.as;

            ad = ad.postUpdateMark();
        end
        
        function color = getNextMarkColor(ad)
            nMarks = ad.nMarks;
            map = TrialDataUtilities.Color.cbrewer('qual', 'Set1');
            color = map(mod(nMarks + 1, size(map, 1)), :);
        end
        
        function ad = removeMarksByIdx(ad, mask)
            ad.warnIfNoArgOut(nargout);

            mask = TensorUtils.vectorIndicesToMask(mask, ad.nMarks);
            % remove all marks and intervals
            ad.warnIfNoArgOut(nargout);
            ad.markEvents = ad.markEvents(~mask);
            ad.markEventsIndex = ad.markEventsIndex(~mask);
            ad.markOffsets = ad.markOffsets(~mask);
            ad.markLabelsStored = ad.markLabelsStored(~mask);
            ad.markAppear= ad.markAppear(~mask);
            ad.markShowOnData = ad.markShowOnData(~mask); 
            ad.markShowOnAxis = ad.markShowOnAxis(~mask);
            
            ad = ad.postUpdateMark();
        end
        
        function ad = removeMarkByString(ad, str)
            ad.warnIfNoArgOut(nargout);

            idx = ad.findMarkByString(str);
            ad = ad.removeMarksByIdx(idx);
        end
        
        function ad = removeMark(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            idx = ad.findMark(eventName, varargin{:});
            ad = ad.removeMarksByIdx(idx);
        end

        function ad = truncateBefore(ad, eventStr, varargin)
            p = inputParser;
            p.addParameter('offset', 0, @isnumeric);
            p.addParameter('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            [eventName, index, offset] = ad.parseEventOffsetString(eventStr, ...
                'truncateBefore', 'defaultIndex', 'end');
            if ~ismember('offset', p.UsingDefaults), offset = p.Results.offset; end
            if ~ismember('index', p.UsingDefaults), index = p.Results.index; end

            ad.assertHasEvent(eventName);
            
            ad.warnIfNoArgOut(nargout);
            ad.truncateBeforeEvents{end+1} = eventName;
            ad.truncateBeforeEventsIndex{end+1} = index;
            ad.truncateBeforeOffsets(end+1) = offset;

            ad.truncateBeforeEvents = makecol(ad.truncateBeforeEvents);
            ad.truncateBeforeEventsIndex = makecol(ad.truncateBeforeEventsIndex);
            ad.truncateBeforeOffsets = makecol(ad.truncateBeforeOffsets);
            
            ad = ad.update();
        end

        function ad = truncateAfter(ad, eventStr, varargin)
            p = inputParser;
            p.addParameter('offset', 0, @isscalar);
            p.addParameter('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            [eventName, index, offset] = ad.parseEventOffsetString(eventStr, ...
                'truncateAfter', 'defaultIndex', '1');
            if ~ismember('offset', p.UsingDefaults), offset = p.Results.offset; end
            if ~ismember('index', p.UsingDefaults), index = p.Results.index; end

            ad.assertHasEvent(eventName);
            
            ad.warnIfNoArgOut(nargout);
            ad.truncateAfterEvents{end+1} = eventName;
            ad.truncateAfterEventsIndex{end+1} = index;
            ad.truncateAfterOffsets(end+1) = offset;
            
            ad.truncateAfterEvents = makecol(ad.truncateAfterEvents);
            ad.truncateAfterEventsIndex = makecol(ad.truncateAfterEventsIndex);
            ad.truncateAfterOffsets = makecol(ad.truncateAfterOffsets);

            ad = ad.update();
        end

        function ad = invalidateOverlap(ad, eventStr, varargin)
            p = inputParser;
            p.addParameter('offset', 0, @isscalar);
            p.addParameter('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            [eventName, index, offset] = ad.parseEventOffsetString(eventStr, ...
                'truncateBefore', 'defaultIndex', ':');
            if ~ismember('offset', p.UsingDefaults), offset = p.Results.offset; end
            if ~ismember('index', p.UsingDefaults), index = p.Results.index; end
            
            ad.assertHasEvent(eventName);

            ad.warnIfNoArgOut(nargout);
            ad.invalidateEvents{end+1} = eventName;
            ad.invalidateEventsIndex{end+1} = index;
            ad.invalidateOffsets(end+1) = offset;
            
            ad.invalidateEvents = makecol(ad.invalidateEvents);
            ad.invalidateEventsIndex = makecol(ad.invalidateEventsIndex);
            ad.invalidateOffsets = makecol(ad.invalidateOffsets);

            ad = ad.update();
        end

%         function ad = setOutsideOfTrialIgnore(ad)
%             ad.warnIfNoArgOut(nargout);
%             ad.outsideOfTrialMode = AlignDescriptor.IGNORE;
%             ad = ad.update();
%         end

        function ad = setOutsideOfTrialInvalidate(ad)
            ad.warnIfNoArgOut(nargout);
            ad.outsideOfTrialMode = AlignDescriptor.INVALIDATE;
            ad = ad.update();
        end

        function ad = setOutsideOfTrialTruncate(ad)
            ad.warnIfNoArgOut(nargout);
            ad.outsideOfTrialMode = AlignDescriptor.TRUNCATE;
            ad = ad.update();
        end

        % store an abbreviation for an event name
        function ad = abbrev(ad, eventName, abbrev)
            ad.warnIfNoArgOut(nargout);
            ad.assertHasEvent(eventName);
            ad.eventAbbrevLookup.(eventName) = abbrev;
        end

        % replace the start and end event with a particular zero-relative
        % time window
        function ad = windowAroundZero(ad, tStart, tStop, varargin)
            ad.warnIfNoArgOut(nargout);
            ad = ad.start(ad.zeroEvent, tStart);
            ad = ad.stop(ad.zeroEvent, tStop);
        end
        
        function ad = round(ad, timeDelta)
            ad.warnIfNoArgOut(nargout);
            ad.roundTimes = true;
            ad.minTimeDelta = timeDelta;
            ad = ad.update();
        end
        
        function ad = noRound(ad)
            ad.warnIfNoArgOut(nargout);
            ad.roundTimes = false;
            ad.minTimeDelta = [];
            ad = ad.update();
        end

        function ad = clearIntervals(ad)
            % remove all marks and intervals
            ad.warnIfNoArgOut(nargout);
            
            ad.intervalEventsStart = {}; 
            ad.intervalEventsStop = {}; 
            ad.intervalEventsIndexStart = [];
            ad.intervalEventsIndexStop  = [];
            ad.intervalOffsetsStart  = [];
            ad.intervalOffsetsStop  = [];
            ad.intervalLabelsStored = {};
            ad.intervalConditionMatch = {}; 
            ad.intervalAppear  = [];
            ad.intervalShowOnData = [];
            ad.intervalShowOnAxis = [];

            ad = ad.postUpdateInterval();
        end
        
        function ad = clearMarks(ad)
            % remove all marks and intervals
            ad.warnIfNoArgOut(nargout);
            ad.markEvents = {};
            ad.markEventsIndex = [];
            ad.markOffsets = [];
            ad.markLabelsStored = {};
            ad.markAppear= [];
            ad.markShowOnData = []; 
            ad.markShowOnAxis = [];
            
            ad = ad.postUpdateMark();
        end
        
        function ad = clearMarksIntervals(ad)
            ad.warnIfNoArgOut(nargout);
            ad = ad.clearMarks();
            ad = ad.clearIntervals();
        end
       
        function ad = lag(ad, tDelay)
            % shift the alignment tDelay units forward in time, instead of
            % start/stop at Event, start/stop at Event-tDelay. We must move
            % the zero so that the time reference for interpolation /
            % re-sampling remains fixed when we lag, otherwise we end up
            % with different numbers of samples
            ad.warnIfNoArgOut(nargout);
            ad.startOffset = ad.startOffset - tDelay;
            ad.stopOffset = ad.stopOffset - tDelay;
            ad.zeroOffset = ad.zeroOffset - tDelay;
            ad = ad.update();
        end
    end

    methods % post-hoc appearance specification
        function ad = setStartAppearance(ad, spec)
            % updates the AppearanceSpec for start
            ad.warnIfNoArgOut(nargout);
            assert(isempty(spec) || isa(spec, 'AppearanceSpec'), 'Must provide AppearanceSpec object or []');
            ad.startAppear = spec;
        end
        
        function ad = setStopAppearance(ad, spec)
            % updates the AppearanceSpec for stop
            ad.warnIfNoArgOut(nargout);
            assert(isempty(spec) || isa(spec, 'AppearanceSpec'), 'Must provide AppearanceSpec object or []');
            ad.stopAppear = spec;
        end
        
        function ad = setZeroAppearance(ad, spec)
            % updates the AppearanceSpec for zero
            ad.warnIfNoArgOut(nargout);
            assert(isempty(spec) || isa(spec, 'AppearanceSpec'), 'Must provide AppearanceSpec object or []');
            ad.zeroAppear = spec;
        end
        
        function ad = setMarkAppearance(ad, ind, spec)
            % updates the AppearanceSpec for mark at index ind
            % use .findMark to find the index for a given mark's event,
            ad.warnIfNoArgOut(nargout);
            assert(isempty(spec) || isa(spec, 'AppearanceSpec'), 'Must provide AppearanceSpec object or []');
            assert(ind > 0 && ind < ad.nMarks, 'Index out of range');
            ad.markAppear{ind} = spec;
        end
        
        function ad = setIntervalAppearance(ad, ind, spec)
            % updates the AppearanceSpec for interval at index ind
            % use .findInterval to find the index for a given interval's
            % events
            ad.warnIfNoArgOut(nargout);
            assert(isempty(spec) || isa(spec, 'AppearanceSpec'), 'Must provide AppearanceSpec object or []');
            assert(ind > 0 && ind < ad.nIntervals, 'Index out of range');
            ad.intervalAppear{ind} = spec;
        end
    end
    
    methods % Equivalence to other descriptors
        % determine if this align descriptor is functionally equivalent to another AlignDescriptor,
        % ignoring appearance-level details (e.g. marks, intervals, label info, abbreviations) 
        % and focusing only on time alignment and trial validity details
        function tf = isCompatibleWith(ad1, ad2)
            % check for equivalence of all of the following properties:
            props = { ...
                'startEvent', ...
                'startEventIndex', ...
                'startOffset', ...
                'startEventIndex', ...
                'stopEvent', ...
                'stopEventIndex', ...
                'stopOffset', ...
                'zeroEvent', ...
                'zeroEventIndex', ...
                'zeroOffset', ...
                'truncateBeforeEvents', ...
                'truncateBeforeEventsIndex', ...
                'truncateBeforeOffsets', ...
                'truncateAfterEvents', ...
                'truncateAfterEventsIndex', ...
                'truncateAfterOffsets', ...
                'invalidateEvents', ...
                'invalidateEventsIndex', ...
                'invalidateOffsets', ...
                'outsideOfTrialMode', ...
                'minDuration', ...
                };

            tf = true;
            for i = 1:length(props)
                tf = tf && isequal(ad1.(props{i}), ad2.(props{i}));
            end
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~ishandle(obj)
                warning('%s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.', ...
                    class(obj));
%                 expr = sprintf('debug(''%s'')', message);
%                 evalin('caller', expr); 
            end
        end
    end

    methods(Access=protected) % Defaults, start/stop/zero set dependent on each other if missing
        function ad = setDefaultZero(ad)
            ad.warnIfNoArgOut(nargout);
            % if no zero specified, make zero the start event with no offset
            if isempty(ad.zeroEvent) || ad.zeroDefault
                ad.zeroEvent = ad.startEvent;
                ad.zeroEventIndex = 1;
                ad.zeroOffset = 0;
                ad.zeroLabel = '';
            end
        end
    end

    methods % String parsing and building, disp() 
        function ad = parseDescriptionString(ad, str)
            % parse the string in the format:
            % startEventName+offset : stopEventName+offset @ zeroEventName+offset, truncateBefore event+offset, truncateAfter event+offset, mark event+offset

            ad.warnIfNoArgOut(nargout);

            str = strtrim(str);
            [startStopZero, remain] = strtok(str, ',');

            % split the string into start : stop @ zero
            pat = '(?<start>[^:]+):(?<stop>[^@]+)(?<zero>@.*)?';
            info = regexp(startStopZero, pat, 'names', 'once');

            if isempty(info) || isempty(info.start) || isempty(info.stop)
                % Try parsing it as just a single event, which would be
                % come the zero
                %try 
                    [ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset, ad.zeroLabel] = ...
                        ad.parseEventOffsetString(startStopZero, 'zero', 'defaultIndex', 1);
                    ad.zeroDefault = false;
                    %ad.zeroMark = true;
                    
                    % '@Event' means TrialStart:TrialEnd @ Event (align whole trial to Event), whereas
                    % 'Event' means Event:Event @ Event (single sample at Event)
                    if(str(1) ~= '@')
                        ad = ad.start(ad.zeroEvent, 'index', ad.zeroEventIndex, 'offset', ad.zeroOffset);
                        ad = ad.stop(ad.zeroEvent, 'index', ad.zeroEventIndex, 'offset', ad.zeroOffset);
                    end
                    
                %catch exc
                    % otherwise just fail
                    %error('Error parsing align descriptor %s', string);
                %end
            end

            if ad.zeroDefault
                % parse each of the start, stop, and zero strings
                [ad.startEvent, ad.startEventIndex, ad.startOffset, ad.startLabel] = ...
                    ad.parseEventOffsetString(info.start, 'start', 'defaultIndex', 1);
                ad.startDefault = true;
                [ad.stopEvent, ad.stopEventIndex, ad.stopOffset, ad.stopLabel] = ...
                    ad.parseEventOffsetString(info.stop, 'stop', 'defaultIndex', 1);
                ad.stopDefault = true;
                
                % deal with event:+100, -100:event cases
                if isempty(ad.startEvent) && ~isempty(ad.stopEvent)
                    ad.startEvent = ad.stopEvent;
                elseif isempty(ad.stopEvent) && ~isempty(ad.startEvent)
                    ad.stopEvent = ad.startEvent;
                end

                if isfield(info, 'zero') && ~isempty(info.zero)
                    % zero specified explicitly
                    info.zero = strtrim(info.zero(2:end)); % remove the leading @
                    [ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset, ad.zeroLabel] = ...
                        ad.parseEventOffsetString(info.zero, 'zero', 'defaultIndex', 1);
                    ad.zeroDefault = false;
                    if isempty(ad.startEvent) && isempty(ad.stopEvent)
                        ad.startEvent = ad.zeroEvent;
                        ad.stopEvent = ad.zeroEvent;
                    end 
                else
                    % deal with 500:800 case here
                    if isempty(ad.startEvent) && isempty(ad.stopEvent)
                        ad.startEvent = 'TrialStart';
                        ad.stopEvent = 'TrialStart';
                    end 
                    ad = ad.setDefaultZero();
                end
            end
            
            % now that everything has been manually parsed in a piecemeal
            % fashion, call the functions to ensure everything else gets
            % set up correctly
            ad = ad.start(ad.startEvent, 'index', ad.startEventIndex, 'offset', ad.startOffset);
            ad = ad.stop(ad.stopEvent, 'index', ad.stopEventIndex, 'offset', ad.stopOffset);
            ad = ad.zero(ad.zeroEvent, 'index', ad.zeroEventIndex, 'offset', ad.zeroOffset);

            % parse the remainder strings one by one to get after, before, mark
            while ~isempty(remain)
                % grab the next comma-separated token
                if remain(1) == ','
                    remain = strtrim(remain(2:end));
                end
                [thisToken, remain] = strtok(remain, ','); %#ok<STTOK>

                % parse this as 'typeStr name+offset'
                [typeStr, event, index, offset, label] = ....
                    ad.parseTypedEventOffsetString(thisToken);
                
                if ~isempty(typeStr)
                    % match typeStr against the known types:
                    switch typeStr
                        case {'truncateAfter', 'before'}
                            ad.truncateAfterEvents{end+1} = event;
                            if isempty(index)
                                % choose most conservative default, truncate
                                % everything after the first occurrence
                                ad.truncateAfterEventsIndex{end+1} = 1;
                            else
                                ad.truncateAfterEventsIndex{end+1} = index;
                            end

                            ad.truncateAfterOffsets(end+1) = offset;

                        case {'truncateBefore', 'after'}
                            ad.truncateBeforeEvents{end+1} = event;
                            if isempty(index)
                                % choose most conservative default, truncate
                                % everything before the last occurrence
                                ad.truncateBeforeEventsIndex{end+1} = 'end';
                            else
                                ad.truncateBeforeEventsIndex{end+1} = index;
                            end
                            ad.truncateBeforeOffsets(end+1) = offset;

                        case {'invalidate', 'exclude', 'excluding'}
                            ad.invalidateEvents{end+1} = event;
                            if isempty(index)
                                % choose most conservative default, any
                                % occurrence invalidates the trial
                                ad.invalidateEventsIndex{end+1} = ':';
                            else
                                ad.invalidateEventsIndex{end+1} = index;
                            end

                            ad.invalidateOffsets(end+1) = offset;

                        case {'mark', 'indicate'}
                            ad = ad.mark(event, offset, 'as', label, 'index', index);

                        otherwise
                            error('Unknown descriptor keyword %s', typeStr);
                    end  
                else
                    % try parsing pad command
                    [success, ad.padPre, ad.padPost] = ad.parsePadString(thisToken);
                    if ~success
                        error('Could not parse token %s');
                    end                    
                end
            end

            ad = ad.update();
        end
        
        function [success, padPre, padPost] = parsePadString(ad, str) %#ok<INUSL>
            pat = 'pad\s*\[?(?<padPre>+?-?\d*):(?<padPost>+?-?\d*)';
            str = strtrim(str);
            info = regexp(str, pat, 'names', 'once');
            
            if isempty(info)
                success = false;
            else
                success = true;
                padPre = str2double(info.padPre);
                padPost = str2double(info.padPost);
                
                assert(~isnan(padPre), 'Could not parse padPre %s as number', info.padPre);
                assert(~isnan(padPost), 'Could not parse padPost %s as number', info.padPost);
            end
        end

        function [typeStr, name, index, offset, label] = parseTypedEventOffsetString(ad, str, varargin)
            % looking for something like 'typeStr eventName + offset'
            pat = '(?<typeStr>\w+)\s*(?<offsetStr>.+)';
            str = strtrim(str);
            info = regexp(str, pat, 'names', 'once');

            if isempty(info)
                typeStr = '';
            else
                typeStr = info.typeStr;
                [name, index, offset, label] = ad.parseEventOffsetString(info.offsetStr, typeStr, varargin{:});
            end
        end

        function [name, index, offset, label] = parseEventOffsetString(ad, str, errorName, varargin) %#ok<INUSL>
            % parse a string in a format 'eventName + offset' or 'eventName - offset'
            % or possibly 'eventName + offset as label'
            % or possibly 'eventName(index) ...' where index is 1,2 or
            % end-1, end-2 or 1:end or :
            % errorName provides information for printing error message about what we're trying to parse
            p = inputParser;
            p.addParameter('defaultIndex', [], @isvector);
            p.parse(varargin{:});
            
            % match event name
            eventPat = '(?<event>[a-zA-Z]\w+)?';
            % match anything like V or V:V where V is '#', 'end', 'end-#'
            indexPat = '(?<index>\((end)?-?\d*:?(end)?-?\d*\))?';
            % match anything like +# or -#
            offsetPat = '(?<offset>[+\-\.e\d\s]+)?';
            % match 'as eventLabel'
            asPat = '(as)?\s*(?<label>[\w\s\d+\-\.\<\>\?]+)?';
            
            pat = [eventPat indexPat '\s*' offsetPat '\s*' asPat]; 
            str = strtrim(str);
            info = regexp(str, pat, 'names', 'once');

            if isempty(info)
                error('Error parsing %s descriptor "%s"', errorName, str);
            end

            name = info.event;

            if isempty(info.index)
                index = p.Results.defaultIndex;
            else 
                indexStr = info.index(2:end-1);
                num = str2double(indexStr);
                if ~isnan(num)
                    % convert to number
                    index = num;
                else
                    % keep as string otherwise
                    index = indexStr;
                end
            end

            if isempty(info.offset)
                offset = 0;
            else
                info.offset = info.offset(~isspace(info.offset));
                offset = str2double(info.offset);
                if isnan(offset)
                    error('Error parsing %s offset %s', errorName, info.offset);
                end
            end

            label = info.label;
        end

        function str = buildLabel(ad, name, index, offset)
            if nargin < 3 || isempty(index)
                index = 1;
            end
            if nargin < 4
                offset = 0;
            end

            % determine the appropriate abbreviation
            if isempty(name)
                % just the offset
                abbrev = '';

            elseif isfield(ad.eventAbbrevLookup, name) 
                % abbrev has been specified via addEventAbbrev
                abbrev = ad.eventAbbrevLookup.(name);

            elseif ad.autoAbbreviateLabels 
                % auto abbreviate:
                % here we assume name is WordCased or camelCased, we convert it to the capitalized
                % initials of each word in the event name
                isUpper = upper(name) == name;
                isUpper(1) = true;
                isSpace = name == ' ';
                abbrev = name(isUpper & ~isSpace);
                if ad.autoAbbreviateMakeUpper
                    abbrev = upper(abbrev);
                else
                    abbrev = lower(abbrev);
                end
            else
                % here we assume the name is TitleCased, camelCased, or
                % snake_cased and convert to Spaced Words
                pattern = '([A-Z]*[a-z]+)';
                words = regexp(name, pattern, 'match');
                abbrev = strjoin(upperFirst(words), ' ');
            end

            % build a parenthetical index string if index ~= 1
            if ischar(index)
                indexStr = ['(' index ')'];
            elseif index == 1
                indexStr = '';
            else
                indexStr = ['(' num2str(index) ')'];
            end

            % combine eventAbbrev(index)+offset
            if offset == 0
                str = [abbrev indexStr];
            else
                str = sprintf('%s%s%+.0f', abbrev, indexStr, offset);
            end
        end

        function str = buildUnabbreviatedLabel(ad, name, index, offset) %#ok<INUSL>
            if nargin < 3 || isempty(index)
                index = 1;
            end
            if nargin < 4
                offset = 0;
            end

            % build a parenthetical index string if index ~= 1
            if ischar(index)
                indexStr = ['(' index ')'];
            elseif index == 1
                indexStr = '';
            else
                indexStr = ['(' num2str(index) ')'];
            end

            % combine eventAbbrev(index)+offset
            if offset == 0
                str = [name indexStr];
            else
                str = sprintf('%s%s%+.0f', name, indexStr, offset);
            end
        end

        function str = getDescription(ad)
            desc = ad.getStartStopZeroDescription();
            if strcmp(ad.name, desc)
                str = desc;
            else
                str = sprintf('%s : %s', ad.name, desc);
            end
        end

        function str = getStartStopZeroDescription(ad)
            sStart = sprintf('%s as %s', ad.startUnabbreviatedLabel, ad.startLabel);
            sStop = sprintf('%s as %s', ad.stopUnabbreviatedLabel, ad.stopLabel);
            sZero = sprintf('%s as %s', ad.zeroUnabbreviatedLabel, ad.zeroLabel);

            str = sprintf('%s : %s @ %s', sStart, sStop, sZero);
        end
        
        function str = getStartStopZeroPadDescription(ad)
            if ad.padPre ~= 0 || ad.padPost ~= 0 
                padStr = sprintf(', pad [%d %d]', ad.padPre, ad.padPost);
            else
                padStr = '';
            end
            str = sprintf('%s%s', ad.getStartStopZeroDescription(), padStr);
        end

        function printOneLineDescription(ad)
            desc = ad.getStartStopZeroPadDescription();
            if ~ad.nameDefault
                tcprintf('inline', '{yellow}%s: {bright blue}%s : {none}%s\n', class(ad), ad.name, desc);
            else
                % name will just match desc, so don't print it twice
                tcprintf('inline', '{yellow}%s: {none}%s\n', class(ad), desc);
            end
        end
        
        function printDescription(ad, varargin)
            p = inputParser();
            p.addParameter('active', false, @islogical);
            p.parse(varargin{:});
            
            if ~ad.nameDefault
                if p.Results.active
                    tcprintf('inline', '{yellow}%s: {white}%s {red}(active)\n', class(ad), ad.name);
                else
                    tcprintf('inline', '{yellow}%s: {white}%s\n', class(ad), ad.name);
                end
            else
                % name will just match desc, so don't print it twice
                if p.Results.active
                    tcprintf('inline', '{yellow}%s: {red}(active)\n', class(ad));
                else
                    tcprintf('inline', '{yellow}%s:\n', class(ad));
                end
            end
            
            tcprintf('inline', '  {bright blue}Start {purple}%s {darkGray}as {white}%s\n', ad.startUnabbreviatedLabel, ad.startLabel);
            tcprintf('inline', '  {bright blue}Stop {purple}%s {darkGray}as {white}%s\n', ad.stopUnabbreviatedLabel, ad.stopLabel);
            tcprintf('inline', '  {bright blue}Zero {purple}%s {darkGray}as {white}%s\n', ad.zeroUnabbreviatedLabel, ad.zeroLabel);

            for i = 1:length(ad.markEvents);
                tcprintf('inline', '  {bright blue}Mark {purple}%s {darkGray}as {white}%s\n', ...
                    ad.markUnabbreviatedLabels{i}, ad.markLabels{i}); 
            end
            for i = 1:size(ad.intervalEventsStart, 1)
                tcprintf('inline', '  {bright blue}Interval {purple}%s{darkGray} as {white}%s\n', ...
                    ad.intervalUnabbreviatedLabels{i}, ...
                    ad.intervalLabels{i}); 
                % TODO add condition match description
            end
            for i = 1:length(ad.truncateBeforeEvents);
                tcprintf('inline', '  {bright blue}Truncate before {purple}%s\n', ad.truncateBeforeUnabbreviatedLabels{i});
            end
            for i = 1:length(ad.truncateAfterEvents);
                tcprintf('inline', '  {bright blue}Truncate after {purple}%s\n', ad.truncateAfterUnabbreviatedLabels{i});
            end
            for i = 1:length(ad.invalidateEvents);
                tcprintf('inline', '  {bright blue}Invalidate overlap {purple}%s\n', ad.invalidateUnabbreviatedLabels{i});
            end     
                   
            tcprintf('inline', '  {darkGray}outside trial mode {white}%s\n', ad.outsideOfTrialMode);
            tcprintf('inline', '  {darkGray}minimum duration {white}%d\n', ad.minDuration);
            if ad.roundTimes
                tcprintf('inline', '  {darkGray}round times with min delta {white}%g\n', ad.minTimeDelta);
            end
            if ad.padPre ~= 0 || ad.padPost ~= 0
                tcprintf('inline', '  {darkGray}pad %g pre, %g post\n', ad.padPre, ad.padPost);
            end
        end

        function disp(ad)
            ad.printDescription();
            
            fprintf('\n');
            builtin('disp', ad);
        end
    end

    methods % Utilities for building string descriptors, misc
        function str = buildStringForOffsetFromZero(ad, offset, units)
            if nargin < 3
                units = '';
            else
                units = [' ' units];
            end
            str = sprintf('%s%+g%s', ad.zeroEvent, offset + ad.zeroOffset, units);
        end
    end
    
    methods(Static) % construct from another align descriptor, used primarily by AlignInfo
        function adNew = fromAlignDescriptor(ad, adNew)
            if nargin < 2
                % allow subclasses to provide their own instance to
                % populate
                adNew = AlignDescriptor();
            end

            meta = ?AlignDescriptor;
            props = meta.PropertyList;

            for iProp = 1:length(props)
                prop = props(iProp);
                if prop.Dependent || prop.Constant || prop.Transient
                    continue;
                else
                    adNew.(prop.Name) = ad.(prop.Name);
                end
            end
        end
    end

end
