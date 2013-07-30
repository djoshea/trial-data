classdef AlignDescriptor 
    % Describes a particular way to slice and align trial data in time based off 
    % of event times in each trial. Typically analog data and neural data are 
    % aligned relative to the times of certain events before further analysis takes
    % place. This class provides a common language for expressing these different
    % types of aligning data.

    properties
        omitRedundantLabel = false; % e.g. when zero=event1, start=event1-20, label start as "-20" rather than "event1-20"

        % marks not using the same event as zero will be surrounded by < > 
        % to indicate the mean time is being plotted. if false, they will
        % not be plotted at all.
        markPlotMedians = true;
        
        % for such marks, if the range of relative event times (in either 
        % direction is less than this threshold, the < > marks will be omitted.
        markRelativeDeltaIgnore = 7.5;

        % abbreviate event names to capital/lowercase initials automatically
        % override this by calling addEventAbbrev
        autoAbbreviateLabels = true;
        
        % use capitalized initials instead of lowercased initials
        autoAbbreviateMakeUpper = true;
        
        % description
        name % set manually, or will be auto-populated
        nameDefault = true; % false if name has been set manually
        
        % how to handle the window falling outside of the trial
        outsideOfTrialMode = AlignDescriptor.INVALIDATE; % TRUNCATE or INVALIDATE or IGNORE

        % minimum duration to mark as valid
        minDuration = 0;
        
        % allow an additional time window around the alignment interval
        % see .setPadWindow below
        padWindow = [];
    end
    
    properties(Access=protected)
        % ValueMap : eventName -> abbreviatedName, overrides auto abbreviation
        % use ad = ad.abbrev(eventName, abbrev) to add an entry
        eventAbbrevLookup
    end

    properties(Hidden) % SetAccess=protected
        % slice times >= t_startEvent + startOffset
        startEvent;
        startOffset = 0;
        startLabel = ''; 
        startInfo
        startMarkData % mark this point on a data trace in drawOnData
        startDefault = true;

        % slice times <= t_stopEvent + stopOffset
        stopEvent;
        stopOffset = 0;
        stopLabel = '';
        stopInfo
        stopMarkData 
        stopDefault = true;

        % shift time such that t_zeroEvent + zeroOffset is t=0
        zeroEvent;
        zeroOffset = 0;
        zeroLabel = '';
        zeroInfo
        zeroMarkData
        zeroDefault = true;

        % move interval start forward in time to avoid these times
        truncateBeforeEvents = {};
        truncateBeforeOffsets

        % move interval stop backward in time to avoid these times
        truncateAfterEvents = {};
        truncateAfterOffsets 

        % invalidate any trial whose window includes these events
        invalidateEvents = {};
        invalidateOffsets

        % purely for annotation, used by getTimeInfo and drawTimeAxis methods 
        markEvents = {};
        markOffsets 
        markLabelsStored = {};
        markInfo
        markMarkData

        % for marking time intervals on the time axis as colored rectangles
        intervalEvents = {}; % n x 2 cell array of start/stop events
        intervalOffsets % n x 2 array of offsets
        intervalLabels = {};
        intervalColors = {}; % color of rectangle used for the interval
        % 'first' means use events from first trial in each condition
        intevalMultiTrialMode = 'first';
        intervalConditionMatch = {}; % n x 1 cell array of structs with .attrName = attrValue(s)
        intervalInfo % miscellaneous
    end


    properties(Constant)
        TRUNCATE = 'truncate';
        INVALIDATE = 'invalidate';
        IGNORE = 'ignore';
    end

    properties(Dependent)
        isFixedLength % is the time window guaranteed fixed length
        isStartFixedTime % is the start event always at a fixed time relative to zero
        isStopFixedTime % is the start event always at a fixed time relative to zero
        isMarkFixedTime % is each mark event always at a fixed time relative to zero

        isStartZero % is the start event always == 0
        isStopZero % is the stop event always == 0
        isMarkZero % is each mark event always == 0

        isZeroOutsideStartStop % is the zero event guaranteed to lie outside the start/stop window
    
        markLabels
    end

    methods % Dependent properties
        function tf = get.isFixedLength(ad)
            % the window is always the same length if the start / stop events are the same
            if strcmp(ad.outsideOfTrialMode, ad.INVALIDATE) ||...
               strcmp(ad.outsideOfTrialMode, ad.IGNORE)
                tf = strcmp(ad.startEvent, ad.stopEvent);
            else
                tf = false;
            end
        end

        function tf = get.isStartFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.startEvent);
        end

        function tf = get.isStopFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.stopEvent);
        end

        function tf = get.isMarkFixedTime(ad)
            tf = strcmp(ad.zeroEvent, ad.markEvents);
        end

        function tf = get.isStartZero(ad)
            tf = strcmp(ad.zeroEvent, ad.startEvent) && ad.zeroOffset == ad.startOffset;
        end

        function tf = get.isStopZero(ad)
            tf = strcmp(ad.zeroEvent, ad.stopEvent) && ad.zeroOffset == ad.stopOffset;
        end

        function tf = get.isMarkZero(ad)
            tf = strcmp(ad.zeroEvent, ad.markEvents) & ad.markOffsets == ad.zeroOffset;
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
        
        % auto fill in default labels
        function str = get.startLabel(ad)
            if isempty(ad.startLabel)
                if ad.isStartFixedTime && ~ad.isStartZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', ad.startOffset);
                else
                    str = ad.buildLabel(ad.startEvent, ad.startOffset);
                end
            else
                str = ad.startLabel;
            end
        end
        
        function str = get.stopLabel(ad)
            if isempty(ad.stopLabel)
                if ad.isStopFixedTime && ~ad.isStopZero && ad.omitRedundantLabel 
                    str = ad.buildLabel('', ad.stopOffset);
                else
                    str = ad.buildLabel(ad.stopEvent, ad.stopOffset);
                end
            else
                str = ad.stopLabel;
            end
        end

        function str = get.zeroLabel(ad)
            if isempty(ad.zeroLabel)
%                 if ad.isStartZero
%                     str = '';
%                 elseif ad.isStopZero
%                     str = '';
%                 else
                    str = ad.buildLabel(ad.zeroEvent, ad.zeroOffset);
%                 end
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
                elseif strcmp(ad.markEvents{iMark}, ad.zeroEvent) && ad.omitRedundantLabel
                    % same as zero, just include the offset (e.g. '+100')
                    markLabels{iMark} = ad.buildLabel('', ad.markOffsets(i));
                else
                    % full label string
                    markLabels{iMark} = ad.buildLabel(ad.markEvents{iMark}, ...
                        ad.markOffsets(iMark));
                end
            end
        end
    end

    methods % constructor, manual event specification
        function ad = AlignDescriptor(varargin)
            ad.eventAbbrevLookup = ValueMap('KeyType', 'char', 'ValueType', 'char');
            
            ad.startEvent  = 'start';
            ad.stopEvent = 'end';
            ad.zeroEvent = 'start';
            
            % either use named property value pairs or string syntax
            if nargin == 1 && ischar(varargin{1});
                ad = ad.parseDescriptionString(varargin{1});
            else
                ad = structargs(ad, varargin);
            end
        end

        % return a list of event names this AlignDescriptor references, exluding 'start' and 'end'
        function eventList = getEventList(ad)
           eventList = unique({ad.startEvent ad.stopEvent ad.zeroEvent ...
                    ad.truncateBeforeEvents{:} ad.truncateAfterEvents{:} ad.invalidateEvents{:} ...
                    ad.markEvents{:} ad.intervalEvents{:}});
           eventList = setdiff(eventList, {'start', 'end'});
        end

        function ad = start(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.startEvent = eventName;
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
        end

        function ad = stop(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.stopEvent = eventName;
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
        end

        function ad = zero(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            ad.zeroEvent = eventName;
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
        end

        function ad = interval(ad, eventNameStart, offsetStart, eventNameStop, offsetStop, varargin)
            ad.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addParamValue('as', '', @ischar);
            p.addParamValue('color', 'k', @(x) true);
            p.addParamValue('conditionMatch', struct(), @(x) isstruct(x) && isscalar(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            color = p.Results.color;
            as = p.Results.as;
            conditionMatch = p.Results.conditionMatch;

            iInterval = size(ad.intervalEvents,1)+1;
            ad.intervalEvents{iInterval, 1} = eventNameStart; 
            ad.intervalEvents{iInterval, 2} = eventNameStop; 
            ad.intervalOffsets(iInterval, 1) = offsetStart;
            ad.intervalOffsets(iInterval, 2) = offsetStop;
            ad.intervalColors{iInterval} = color;
            ad.intervalInfo{iInterval} = p.Unmatched;
            ad.intervalConditionMatch{iInterval} = conditionMatch;
            ad.intervalLabels{iInterval} =as;
        end

        function ad = mark(ad, eventName, varargin)
            ad.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addOptional('offset', 0, @isscalar);
            p.addParamValue('as', '', @ischar);
            p.addParamValue('markData', [], @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            offset = p.Results.offset;

            iMark = length(ad.markEvents)+1;
            ad.markEvents{iMark} = eventName;
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
            
        end

        function ad = truncateBefore(ad, eventName, offset, varargin)
            ad.warnIfNoArgOut(nargout);
            ad.truncateBeforeEvents{end+1} = eventName;
            ad.truncateBeforeOffsets(end+1) = offset;
        end

        function ad = truncateAfter(ad, eventName, offset, varargin)
            ad.warnIfNoArgOut(nargout);
            ad.truncateAfterEvents{end+1} = eventName;
            ad.truncateAfterOffsets(end+1) = offset;
        end

        function ad = invalidateOverlap(ad, eventName, offset, varargin)
            ad.warnIfNoArgOut(nargout);
            ad.invalidateEvents{end+1} = eventName;
            ad.invalidateOffsets(end+1) = offset;
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
                'startOffset', ...
                'stopEvent', ...
                'stopOffset', ...
                'zeroEvent', ...
                'zeroOffset', ...
                'truncateBeforeEvents', ...
                'truncateBeforeOffsets', ...
                'truncateAfterEvents', ...
                'truncateAfterOffsets', ...
                'invalidateEvents', ...
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
                ad.zeroOffset = 0;
            end
        end
    end

    methods % String parsing and building, disp() 
        function ad = parseDescriptionString(ad, str)
            % parse the string in the format:
            % startEventName+offset : stopEventName+offset @ zeroEventName+offset, truncateBefore event+offset, truncateAfter event+offset, mark event+offset
  
            [startStopZero remain] = strtok(str, ',');

            % split the string into start : stop @ zero
            pat = '(?<start>[^:]+):(?<stop>[^@]+)(?<zero>@.*)?';
            info = regexp(startStopZero, pat, 'names', 'once');

            if isempty(info) || isempty(info.start) || isempty(info.stop)
                % Try parsing it as just a single event, which would be
                % come the zero
                try 
                    [ad.zeroEvent ad.zeroOffset ad.zeroLabel] = ad.parseEventOffsetString(startStopZero);
                    ad.zeroDefault = false;
                catch exc
                    % otherwise just fail
                    error('Error parsing align descriptor %s', string);
                end
            end

            if ad.zeroDefault
                % parse each of the start, stop, and zero strings
                [ad.startEvent ad.startOffset ad.startLabel] = ad.parseEventOffsetString(info.start, 'start');
                ad.startDefault = true;
                [ad.stopEvent ad.stopOffset ad.stopLabel] = ad.parseEventOffsetString(info.stop, 'stop');
                ad.stopDefault = true;

                if isfield(info, 'zero') && ~isempty(info.zero)
                    % zero specified explicitly
                    info.zero = strtrim(info.zero(2:end)); % remove the leading @
                    [ad.zeroEvent ad.zeroOffset ad.zeroLabel] = ad.parseEventOffsetString(info.zero, 'zero');
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
                [typeStr name offset label] = ad.parseTypedEventOffsetString(thisToken);

                % match typeStr against the known types:
                switch typeStr
                    case {'truncateAfter', 'after'}
                        ad.truncateAfterEvents{end+1} = name;
                        ad.truncateAfterOffsets(end+1) = offset;

                    case {'truncateBefore', 'before'}
                        ad.truncateBeforeEvents{end+1} = name;
                        ad.truncateBeforeOffsets(end+1) = offset;

                    case 'invalidate'
                        ad.invalidateEvents{end+1} = name;
                        ad.invalidateOffsets(end+1) = offset;

                    case 'mark'
                        ad = ad.mark(name, offset, 'as', label);

                    otherwise
                        error('Unknown descriptor keyword %s', typeStr);
                end    
            end
        end

        function [typeStr name offset label] = parseTypedEventOffsetString(ad, str)
            % looking for something like 'typeStr eventName + offset'
            pat = '(?<typeStr>\w+)\s*(?<offsetStr>.+)';
            str = strtrim(str);
            info = regexp(str, pat, 'names', 'once');

            if isempty(info)
                error('Error parsing offset string %s', str);
            end

            typeStr = info.typeStr;
            [name offset label] = ad.parseEventOffsetString(info.offsetStr, typeStr);
        end

        function [name offset label] = parseEventOffsetString(ad, str, errorName)
            % parse a string in a format 'eventName + offset' or 'eventName - offset'
            % or possibly 'eventName + offset as label'
            % errorName provides information for printing error message about what we're trying to parse

            str = strtrim(str);
            pat = '(?<event>\w+)\s*(?<offset>[+\-\.e\d\s]+)?\s*(as)?\s*(?<label>[\w\s\d+\-\.\<\>\?]+)?'; 
            info = regexp(str, pat, 'names', 'once');
            
            if isempty(info)
                error('Error parsing %s descriptor "%s"', errorName, str);
            end

            name = info.event;

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

        function str = buildLabel(ad, name, offset)
            if nargin < 3
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
                
            if offset == 0
                str = abbrev;
            else
                str = sprintf('%s%+.0f', abbrev, offset);
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
            if ad.startOffset == 0
                sStart = sprintf('%s as %s', ad.startEvent, ad.startLabel);
            else
                sStart = sprintf('%s%+.0f as %s', ad.startEvent, ad.startOffset, ad.startLabel);
            end

            if ad.stopOffset == 0
                sStop = sprintf('%s as %s', ad.stopEvent, ad.stopLabel);
            else
                sStop = sprintf('%s%+.0f as %s', ad.stopEvent, ad.stopOffset, ad.stopLabel);
            end

            if ad.zeroOffset == 0
                sZero = sprintf('%s as %s', ad.zeroEvent, ad.zeroLabel);
            else
                sZero= sprintf('%s%+.0f as %s', ad.zeroEvent, ad.zeroOffset, ad.zeroLabel);
            end

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

            for i = 1:length(ad.markEvents);
                tcprintf('inline', '\tmark {white}%s%+.0f{none} as {bright blue}%s\n', ad.markEvents{i}, ...
                    ad.markOffsets(i), ad.markLabels{i}); 
            end
            for i = 1:size(ad.intervalEvents, 1)
                tcprintf('inline', '\tinterval {white}%s%+.0f{none} : {white}%s%+0.f{none} as {bright blue}%s\n', ...
                    ad.intervalEvents{i, 1}, ad.intervalOffsets(i, 1), ...
                    ad.intervalEvents{i, 2}, ad.intervalOffsets(i, 2), ...
                    ad.intervalLabels{i}); 
                % TODO add condition match description
            end
            for i = 1:length(ad.truncateBeforeEvents);
                tcprintf('inline', '\ttruncate before {white}%s%+.0f{none}\n', ad.truncateBeforeEvents{i}, ...
                    ad.truncateBeforeOffsets(i)); 
            end
            for i = 1:length(ad.truncateAfterEvents);
                tcprintf('inline', '\ttruncate after {white}%s%+.0f{none}\n', ad.truncateAfterEvents{i}, ...
                    ad.truncateAfterOffsets(i)); 
            end
            for i = 1:length(ad.invalidateEvents);
                tcprintf('inline', '\tinvalidate overlap {white}%s%+.0f{none}\n', ad.invalidateEvents{i}, ...
                    ad.invalidateOffsets(i)); 
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
