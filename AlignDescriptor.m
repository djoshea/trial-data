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
        autoAbbreviateLabels = true;
        
        % use capitalized initials instead of lowercased initials
        autoAbbreviateMakeUpper = true;
        
        % for description
        name % set manually, or will be auto-populated
    end
        
    properties(SetAccess=protected, Hidden)
        nameDefault = true; % false if name has been set manually

        % how to handle the window falling outside of the trial
        outsideOfTrialMode = AlignDescriptor.TRUNCATE; % TRUNCATE or INVALIDATE or IGNORE

        % minimum duration to mark as valid
        minDuration = 0;
        
        % allow an additional time window around the alignment interval
        % see .setPadWindow below
        padPre = 0;
        padPost = 0;
    end
    
    properties(Access=protected)
        % ValueMap : eventName -> abbreviatedName, overrides auto abbreviation
        % use ad = ad.abbrev(eventName, abbrev) to add an entry
        eventAbbrevLookup
    end

    properties(Hidden, SetAccess=protected) % SetAccess=protected
        % slice times >= t_startEvent + startOffset
        startEvent;
        startEventIndex = 1;
        startOffset = 0;
        startLabel = ''; 
        startInfo
        startMarkData % mark this point on a data trace in drawOnData
        startDefault = true;

        % slice times <= t_stopEvent + stopOffset
        stopEvent;
        stopOffset = 0;
        stopEventIndex = 1;
        stopLabel = '';
        stopInfo
        stopMarkData 
        stopDefault = true;

        % shift time such that t_zeroEvent + zeroOffset is t=0
        zeroEvent;
        zeroOffset = 0;
        zeroEventIndex = 1;
        zeroLabel = '';
        zeroInfo
        zeroMarkData
        zeroDefault = true;

        % move interval start forward in time to avoid these times
        truncateBeforeEvents = {};
        truncateBeforeEventsIndex
        truncateBeforeOffsets

        % move interval stop backward in time to avoid these times
        truncateAfterEvents = {};
        truncateAfterEventsIndex
        truncateAfterOffsets 

        % invalidate any trial whose window includes these events
        invalidateEvents = {};
        invalidateEventsIndex
        invalidateOffsets

        % purely for annotation, used by getTimeInfo and drawTimeAxis methods 
        markEvents = {};
        markEventsIndex
        markOffsets 
        markLabelsStored = {};
        markInfo
        markMarkData

        % for marking time intervals on the time axis as colored rectangles
        intervalEvents = {}; % n x 2 cell array of start/stop events
        intervalEventsIndex % n x 2 array of event indices
        intervalOffsets % n x 2 array of offsets
        intervalLabelsStored = {};
        intervalColors = {}; % color of rectangle used for the interval

        % 'first' means use events from first trial in each condition
        %intevalMultiTrialMode = 'first';
        intervalConditionMatch = {}; % n x 1 cell array of structs with .attrName = attrValue(s)
        intervalInfo % miscellaneous
    end

    properties(Constant, Hidden)
        TRUNCATE = 'truncate';
        INVALIDATE = 'invalidate';
        IGNORE = 'ignore';
    end

    properties(Dependent, Hidden)
        padWindow
        markLabels
        intervalLabels

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

        isZeroOutsideStartStop % is the zero event guaranteed to lie outside the start/stop window
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
            tf = strcmp(ad.zeroEvent, ad.intervalEvents) & isequal(ad.zeroEventIndex, ad.intervalEventsIndex);
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
            tf = strcmp(ad.zeroEvent, ad.markEvents) & isequal(ad.markEventsIndex, ad.stopEventIndex) & ...
                ad.markOffsets == ad.zeroOffset;
        end

        function tf = get.isZeroOutsideStartStop(ad)
            tf = (ad.isStartFixedTime && ad.startOffset > 0) || (ad.isStopFixedTime && ad.stopOffset < 0);
        end
        
        function name = get.name(ad)
            if isempty(ad.name)
                name = ad.getStartStopZeroDescription();
            else
                name = ad.name;
            end
        end
        
        function ad = set.name(ad, name)
            ad.name = name;
            ad.nameDefault = isempty(name);
        end

        function padWindow = get.padWindow(ad)
            padWindow = [ad.padPre ad.padPost];
        end
    end
    
    methods % Auto-generated / manually set abbreviated labels
        function str = get.startLabel(ad)
            if isempty(ad.startLabel)
                if ad.isStartFixedTime && ~ad.isStartZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', [], ad.startOffset);
                else
                    str = ad.buildLabel(ad.startEvent, ad.startEventIndex, ad.startOffset);
                end
            else
                str = ad.startLabel;
            end
        end
        
        function str = get.stopLabel(ad)
            if isempty(ad.stopLabel)
                if ad.isStopFixedTime && ~ad.isStopZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', [], ad.stopOffset);
                else
                    str = ad.buildLabel(ad.stopEvent, ad.stopEventIndex, ad.stopOffset);
                end
            else
                str = ad.stopLabel;
            end
        end

        function str = get.zeroLabel(ad)
            if isempty(ad.zeroLabel)
                str = ad.buildLabel(ad.zeroEvent, ad.zeroEventIndex, ad.zeroOffset);
            else
                str = ad.zeroLabel;
            end
        end

        function markLabels = get.markLabels(ad)
            if isempty(ad.markLabelsStored)
                ad.markLabelsStored = cell(length(ad.markEvents), 1);
            end
            
            markLabels = cell(length(ad.markLabelsStored), 1);
            for iMark = 1:length(ad.markLabelsStored)
                if ~isempty(ad.markLabelsStored{iMark})
                    % manual label
                    markLabels{iMark} = ad.markLabelsStored{iMark};
                elseif ad.isMarkFixedTime(iMark) && ad.omitRedundantLabel
                    % same as zero, just include the offset (e.g. '+100')
                    markLabels{iMark} = ad.buildLabel('', [], ad.markOffsets(i));
                else
                    % full label string
                    markLabels{iMark} = ad.buildLabel(ad.markEvents{iMark}, ...
                        ad.markEventsIndex{iMark}, ad.markOffsets(iMark));
                end
            end
        end

        function intervalLabels = get.intervalLabels(ad)
            nIntervals = size(ad.intervalEvents, 1);
            if isempty(ad.intervalLabelsStored)
                ad.intervalLabelsStored = cell(nIntervals, 1);
            end
            
            intervalLabels = cell(length(ad.intervalLabelsStored), 1);
            for i = 1:nIntervals
                if ~isempty(ad.intervalLabelsStored{i})
                    % manual label
                    intervalLabels{i} = ad.intervalLabelsStored{i};
                else
                    if ad.isIntervalFixedTime(i, 1) && ad.omitRedundantLabel
                        % same as zero, just include the offset (e.g. '+100')
                        label1 = ad.buildLabel('', [], ad.intervalOffsets(i, 1));
                    else
                        label1 = ad.buildLabel(ad.intervalEvents{i, 1}, ad.intervalEventsIndex{i,1}, ad.intervalOffsets(i, 1));
                    end
                    if ad.isIntervalFixedTime(i, 2) && ad.omitRedundantLabel
                        % same as zero, just include the offset (e.g. '+100')
                        label2 = ad.buildLabel('', [], ad.intervalOffsets(i, 2));
                    else
                        label2 = ad.buildLabel(ad.intervalEvents{i, 2}, ad.intervalEventsIndex{i,2}, ad.intervalOffsets(i, 2));
                    end

                    intervalLabels{i} = sprintf('%s : %s', label1, label2);
                end
            end
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
            strCell = cellfun(@ad.buildUnabbreviatedLabel, ad.intervalEvents, ...
                ad.intervalEventsIndex, num2cell(ad.intervalOffsets), 'UniformOutput', false);
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
            ad.eventAbbrevLookup = ValueMap('KeyType', 'char', 'ValueType', 'char');

            ad.startEvent  = 'TrialStart';
            ad.stopEvent = 'TrialEnd';
            ad.zeroEvent = 'TrialStart';

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

        % return a list of event names this AlignDescriptor references, plus 'start' and 'end'
        function eventList = getEventList(ad)
            eventList = unique({ad.startEvent ad.stopEvent ad.zeroEvent ...
                ad.truncateBeforeEvents{:} ad.truncateAfterEvents{:} ad.invalidateEvents{:} ...
                ad.markEvents{:} ad.intervalEvents{:}});
            eventList = union(eventList, {'TrialStart', 'TrialEnd'});
        end

        function ad = start(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.startEvent = eventName;
            ad.startEventIndex = p.Results.index;
            ad.startOffset = offset;
            if isempty(p.Results.markData)
                ad.startMarkData = offset == 0;
            else
                ad.startMarkData = p.Results.markData;
            end

            if ~isempty(p.Results.as)
                ad.startLabel = p.Results.as;
            end

            ad.startInfo = p.Unmatched;
            ad.startDefault = false;

            % if no zero is specified, determine it from the start event
            ad = ad.setDefaultZero();

            ad = ad.update();
        end

        function ad = stop(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.stopEvent = eventName;
            ad.stopEventIndex = p.Results.index;
            ad.stopOffset = offset;
            if isempty(p.Results.markData)
                ad.stopMarkData = offset == 0;
            else
                ad.stopMarkData = p.Results.markData;
            end

            if ~isempty(p.Results.as)
                ad.stopLabel = p.Results.as;
            end
            ad.stopInfo = p.Unmatched;
            ad.stopDefault = true;

            ad = ad.update();
        end

        function ad = pad(ad, window, varargin)
            ad.warnIfNoArgOut(nargout);
            ad.padPre = window(1);
            ad.padPost = window(2);

            ad = ad.update();
        end

        function ad = zero(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('index', 1, @(x) ischar(x) || isscalar(x));
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.zeroEvent = eventName;
            ad.zeroEventIndex = p.Results.index;
            ad.zeroOffset = offset;
            if isempty(p.Results.markData)
                ad.zeroMarkData = offset == 0;
            else
                ad.zeroMarkData = p.Results.markData;
            end

            if ~isempty(p.Results.as)
                ad.zeroLabel = p.Results.as;
            end
            ad.zeroInfo = p.Unmatched;
            ad.zeroDefault = true;

            ad = ad.update();
        end

        function ad = interval(ad, eventNameStart, offsetStart, eventNameStop, offsetStop, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addParamValue('indexStart', 1, @(x) ischar(x) || isscalar(x));
            p.addParamValue('indexEnd', 1, @(x) ischar(x) || isscalar(x));
            p.addParamValue('as', '', @ischar);
            p.addParamValue('color', 'k', @(x) true);
            p.addParamValue('conditionMatch', struct(), @(x) isstruct(x) && isscalar(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            color = p.Results.color;
            as = p.Results.as;
            conditionMatch = p.Results.conditionMatch;

            iInterval = size(ad.intervalEvents,1)+1;
            ad.intervalEventsStart{iInterval} = eventNameStart; 
            ad.intervalEventsStop{iInterval} = eventNameStop; 
            ad.intervalEventsIndexStart{iInterval} = p.Results.indexStart;
            ad.intervalEventsIndexStop{iInterval} = p.Results.indexEnd;
            ad.intervalOffsetsStart(iInterval) = offsetStart;
            ad.intervalOffsetsStop(iInterval) = offsetStop;
            ad.intervalColors{iInterval} = color;
            ad.intervalInfo{iInterval} = p.Unmatched;
            %ad.intervalConditionMatch{iInterval} = conditionMatch;
            ad.intervalLabelsStored{iInterval} =as;

            ad = ad.update();
        end

        function ad = mark(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;
            
            if isempty(index)
                index = ':';
            else
                index = p.Results.index;
            end

            iMark = length(ad.markEvents)+1;
            ad.markEvents{iMark} = eventName;
            ad.markEventsIndex{iMark} = index;
            ad.markOffsets(iMark) = offset;
            if isempty(p.Results.markData)
                ad.markMarkData(iMark) = offset == 0;
            else
                ad.markMarkData(iMark) = p.Results.markData;
            end
            ad.markInfo{iMark} = p.Unmatched;

            if ~isempty(p.Results.as)
                % store manual label
                ad.markLabelsStored{iMark,1} = p.Results.as;
            else
                ad.markLabelsStored{iMark,1} = '';
            end

            ad = ad.update();
        end

        function ad = truncateBefore(ad, eventName, offset, varargin)
            p = inputParser;
            p.addOptional('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            if isempty(p.Results.index)
                % choose most conservative default
                index = 'end';
            else
                index = p.Results.index;
            end

            ad.warnIfNoArgOut(nargout);
            ad.truncateBeforeEvents{end+1} = eventName;
            ad.truncateBeforeEventsIndex{end+1} = index;
            ad.truncateBeforeOffsets(end+1) = offset;

            ad = ad.update();
        end

        function ad = truncateAfter(ad, eventName, offset, varargin)
            p = inputParser;
            p.addOptional('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            if isempty(p.Results.index)
                index = 1;
            else
                index = p.Results.index;
            end

            ad.warnIfNoArgOut(nargout);
            ad.truncateAfterEvents{end+1} = eventName;
            ad.truncateAfterEventsIndex{end+1} = index;
            ad.truncateAfterOffsets(end+1) = offset;

            ad = ad.update();
        end

        function ad = invalidateOverlap(ad, eventName, offset, varargin)
            p = inputParser;
            p.addOptional('index', [], @(x) isempty(x) || ischar(x) || isscalar(x));
            p.parse(varargin{:});
            
            if isempty(p.Results.index)
                index = ':';
            else
                index = p.Results.index;
            end

            ad.warnIfNoArgOut(nargout);
            ad.invalidateEvents{end+1} = eventName;
            ad.invalidateEventsIndex{end+1} = index;
            ad.invalidateOffsets(end+1) = offset;

            ad = ad.update();
        end

        function ad = setOutsideOfTrialIgnore(ad)
            ad.warnIfNoArgOut(nargout);
            ad.outsideOfTrialMode = AlignDescriptor.IGNORE;
            ad = ad.update();
        end

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
            ad.eventAbbrevLookup(eventName) = abbrev;
        end

        % replace the start and end event with a particular zero-relative
        % time window
        function ad = windowAroundZero(ad, tStart, tStop, varargin)
            ad.warnIfNoArgOut(nargout);
            ad = ad.start(ad.zeroEvent, tStart); 
            ad = ad.stop(ad.zeroEvent, tStop);
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
                message = sprintf('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
                expr = sprintf('debug(''%s'')', message);
                evalin('caller', expr); 
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
            end
        end
    end

    methods % String parsing and building, disp() 
        function ad = parseDescriptionString(ad, str)
            % parse the string in the format:
            % startEventName+offset : stopEventName+offset @ zeroEventName+offset, truncateBefore event+offset, truncateAfter event+offset, mark event+offset

            ad.warnIfNoArgOut(nargout);

            [startStopZero remain] = strtok(str, ',');

            % split the string into start : stop @ zero
            pat = '(?<start>[^:]+):(?<stop>[^@]+)(?<zero>@.*)?';
            info = regexp(startStopZero, pat, 'names', 'once');

            if isempty(info) || isempty(info.start) || isempty(info.stop)
                % Try parsing it as just a single event, which would be
                % come the zero
                %try 
                    [ad.zeroEvent ad.zeroEventIndex ad.zeroOffset ad.zeroLabel] = ...
                        ad.parseEventOffsetString(startStopZero, 'zero', 'defaultIndex', 1);
                    ad.zeroDefault = false;
                    ad = ad.start(ad.zeroEvent, 'index', ad.zeroEventIndex, 'offset', ad.zeroOffset);
                    ad = ad.stop(ad.zeroEvent, 'index', ad.zeroEventIndex, 'offset', ad.zeroOffset);
                %catch exc
                    % otherwise just fail
                    %error('Error parsing align descriptor %s', string);
                %end
            end

            if ad.zeroDefault
                % parse each of the start, stop, and zero strings
                [ad.startEvent ad.startEventIndex ad.startOffset ad.startLabel] = ...
                    ad.parseEventOffsetString(info.start, 'start', 'defaultIndex', 1);
                ad.startDefault = true;
                [ad.stopEvent ad.stopEventIndex ad.stopOffset ad.stopLabel] = ...
                    ad.parseEventOffsetString(info.stop, 'stop', 'defaultIndex', 1);
                ad.stopDefault = true;

                if isfield(info, 'zero') && ~isempty(info.zero)
                    % zero specified explicitly
                    info.zero = strtrim(info.zero(2:end)); % remove the leading @
                    [ad.zeroEvent ad.zeroEventIndex ad.zeroOffset ad.zeroLabel] = ...
                        ad.parseEventOffsetString(info.zero, 'zero', 'defaultIndex', 1);
                    ad.zeroDefault = false;
                else
                    ad = ad.setDefaultZero();
                end
            end

            % parse the remainder strings one by one to get after, before, mark
            while ~isempty(remain)
                % grab the next comma-separated token
                if remain(1) == ','
                    remain = strtrim(remain(2:end));
                end
                [thisToken remain] = strtok(remain, ',');

                % parse this as 'typeStr name+offset'
                [typeStr name index offset label] = ....
                    ad.parseTypedEventOffsetString(thisToken);

                % match typeStr against the known types:
                switch typeStr
                    case {'truncateAfter', 'before'}
                        ad.truncateAfterEvents{end+1} = name;
                        if isempty(index)
                            % choose most conservative default, truncate
                            % everything after the first occurrence
                            ad.truncateAfterEventsIndex{end+1} = '1';
                        else
                            ad.truncateAfterEventsIndex{end+1} = index;
                        end
                            
                        ad.truncateAfterOffsets(end+1) = offset;

                    case {'truncateBefore', 'after'}
                        ad.truncateBeforeEvents{end+1} = name;
                        if isempty(index)
                            % choose most conservative default, truncate
                            % everything before the last occurrence
                            ad.truncateBeforeEventsIndex{end+1} = 'end';
                        else
                            ad.truncateBeforeEventsIndex{end+1} = index;
                        end
                        ad.truncateBeforeOffsets(end+1) = offset;

                    case {'invalidate', 'exclude', 'excluding'}
                        ad.invalidateEvents{end+1} = name;
                        if isempty(index)
                            % choose most conservative default, any
                            % occurrence invalidates the trial
                            ad.invalidateEventsIndex{end+1} = ':';
                        else
                            ad.invalidateEventsIndex{end+1} = index;
                        end
                            
                        ad.invalidateOffsets(end+1) = offset;

                    case {'mark', 'indicate'}
                        ad = ad.mark(name, offset, 'as', label, 'index', index);

                    otherwise
                        error('Unknown descriptor keyword %s', typeStr);
                end    
            end

            ad = ad.update();
        end

        function [typeStr name index offset label] = parseTypedEventOffsetString(ad, str, varargin)
            % looking for something like 'typeStr eventName + offset'
            pat = '(?<typeStr>\w+)\s*(?<offsetStr>.+)';
            str = strtrim(str);
            info = regexp(str, pat, 'names', 'once');

            if isempty(info)
                error('Error parsing offset string %s', str);
            end

            typeStr = info.typeStr;
            [name index offset label] = ad.parseEventOffsetString(info.offsetStr, typeStr, varargin{:});
        end

        function [name index offset label] = parseEventOffsetString(ad, str, errorName, varargin)
            % parse a string in a format 'eventName + offset' or 'eventName - offset'
            % or possibly 'eventName + offset as label'
            % or possibly 'eventName(index) ...' where index is 1,2 or
            % end-1, end-2 or 1:end or :
            % errorName provides information for printing error message about what we're trying to parse
            p = inputParser;
            p.addParamValue('defaultIndex', [], @isvector);
            p.parse(varargin{:});
            
            % match event name
            eventPat = '(?<event>\w+)';
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

            elseif ad.eventAbbrevLookup.isKey(name) 
                % abbrev has been specified via addEventAbbrev
                abbrev = ad.eventAbbrevLookup(name);

            elseif ad.autoAbbreviateLabels 
                % auto abbreviate:
                % here we assume name is camelCased, we convert it to the capitalized
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

        function str = buildUnabbreviatedLabel(ad, name, index, offset)
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

        function printOneLineDescription(ad)
            desc = ad.getStartStopZeroDescription();
            if ~ad.nameDefault
                tcprintf('inline', '{yellow}%s: {white}%s : {bright blue}%s\n', class(ad), ad.name, desc);
            else
                % name will just match desc, so don't print it twice
                tcprintf('inline', '{yellow}%s: {bright blue}%s\n', class(ad), desc);
            end
        end

        function disp(ad)
            ad.printOneLineDescription();
            
            tcprintf('inline', '\tpad {bright blue}[%d %d]\n', ad.padPre, ad.padPost);

            tcprintf('inline', '\toutside trial {bright blue}%s\n', ad.outsideOfTrialMode);
            
            tcprintf('inline', '\tminimum duration {bright blue}%d\n', ad.minDuration);

            for i = 1:length(ad.markEvents);
                tcprintf('inline', '\tmark {white}%s{none} as {bright blue}%s\n', ...
                    ad.markUnabbreviatedLabels{i}, ad.markLabels{i}); 
            end
            for i = 1:size(ad.intervalEvents, 1)
                tcprintf('inline', '\tinterval {white}%s{none} : {white}%s{none} as {bright blue}%s\n', ...
                    ad.intervalUnabbreviatedLabels{i, 1}, ...
                    ad.intervalUnabbreviatedLabels{i, 1}, ...
                    ad.intervalLabels{i}); 
                % TODO add condition match description
            end
            for i = 1:length(ad.truncateBeforeEvents);
                tcprintf('inline', '\ttruncate before {white}%s{none}\n', ad.truncateBeforeUnabbreviatedLabels{i});
            end
            for i = 1:length(ad.truncateAfterEvents);
                tcprintf('inline', '\ttruncate after {white}%s%+.0f{none}\n', ad.truncateAfterUnabbreviatedLabels{i});
            end
            for i = 1:length(ad.invalidateEvents);
                tcprintf('inline', '\tinvalidate overlap {white}%s{none}\n', ad.invalidateUnabbreviatedLabels{i});
            end

            fprintf('\n');
            builtin('disp', ad);
        end
    end

    methods % Label information
        function [labelInfo] = getLabelInfo(ad, varargin)
            % build a list of label names / times to mark on the time axis 
            % when using this alignment. Essentially, all events which are fixed relative
            % to zero (i.e. reference off the same event) will be included as labels
            % also, if 'tMin' and 'tMax' are specified, they will be drawn as well provided that 
            % the start and stop events aren't fixed relative to zero 

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
                end
            end

            if ~isempty(tMin)
                % start not fixed, include one for the lower limit tMin
                % first check whether there is an existing label (from a
                % mark) at this point already...
                if ~any(floor([labelInfo.time]) == floor(tMin))
                    labelInfo(counter).name = ad.buildLabel(ad.zeroEvent, tMin);
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
                    labelInfo(counter).name = ad.buildLabel(ad.zeroEvent, tMax);
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

        function lims = getTimeAxisLims(ad, timeInfo, varargin)
            labelInfo = ad.getLabelInfo(timeInfo);
            tickPos = [labelInfo.time];
            lims = [min(tickPos) max(tickPos)];
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
                    name = prop.Name;
                    adNew.(name) = ad.(name);
                end
            end
        end
    end

end
