function Q = addSegmentedEvents(Q, eventInfo, eventData)
	
    % strip off the low order byte if requested
    if eventInfo.lowOrderByte
        eventData.codeByte = bitand(eventData.code, 255);
    else
        eventData.codeByte = eventData.code;
    end
    
    % convert all of the events over
    if isfield(eventInfo, 'codes')
        codeNames = fieldnames(eventInfo.codes);
        codeValues = cellfun(@(fld) eventInfo.codes.(fld), codeNames);
    end
   
	prog = ProgressBar(length(Q), 'Adding segmented event data');
    for iq = 1:length(Q)
    	prog.update(iq);

		% grab the events that occurred within this trial (respecting the pre/post skips)
		% round is to ensure only ms precision is used, as there is some noise in determining trial
        % start and end boundaries and it is okay for events, since they'll likely be clustered
        % near the start boundary but not near the end boundary
        % CHANGING THIS TO FLOOR ON 2013-04-11 @djoshea
        eventsWithin = find(eventData.timestamp >= Q(iq).CerebusInfo.startTime & ...
                            eventData.timestamp < Q(iq).CerebusInfo.endTime);
                   
        Q(iq).evc.rawPulses = eventData.code(eventsWithin);
        
        % now cut off the beginning and ending pulses as requested
        eventsWithin = eventsWithin(eventInfo.skipPulsesStart+1:end-eventInfo.skipPulsesEnd);
        eventCodes = eventData.codeByte(eventsWithin);
        eventTimes = eventData.timestamp(eventsWithin) - Q(iq).CerebusInfo.startTime;

        Q(iq).evc.codes = eventCodes;
        Q(iq).evc.times = eventTimes;

        % for each code, assign a field into q.evc.codeName = [times list]
        if isfield(eventInfo, 'codes')
            for ic = 1:length(codeNames)
                Q(iq).evc.(codeNames{ic}) = eventTimes(eventCodes == codeValues(ic));
            end
        end
    end

	prog.finish(); 
end

