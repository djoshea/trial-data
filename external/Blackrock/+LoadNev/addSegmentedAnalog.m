function [Q, analogLookup] = addSegmentedAnalog(Q, analogInfo, nsxData)
% grab analog data according to analogLookup,
% subtract time offset, parse into channel groups,
% assign into struct

% build the lookup table which tells us where to find (in nsxData) each
% requested analog channel to speed things up a biq

    analogLookup = LoadNev.buildAnalogLookup(analogInfo, nsxData);

    if isempty(Q)
        return;
    end

    % we first want to split the data into multiple facets
    % grab the start and end inds of each trial
    CI = [Q.CerebusInfo];
    startTimes = [CI.startTime];
    endTimes = [CI.endTime];
    
    [timeCell, dataCell] = deal(cell(length(nsxData), 1));

    for insx = 1:length(nsxData)
        % find start and end inds
        
        prog = ProgressBar(length(Q), 'Finding segmented analog data times');
        [startInds, endInds] = deal(nan(length(Q), 1));
        for iq = 1:length(Q)
            temp = find(nsxData(insx).time >= startTimes(iq), 1, 'first');
            if ~isempty(temp)
                startInds(iq) = temp;
            end
            temp = find(nsxData(insx).time <= endTimes(iq), 1, 'last');
            if ~isempty(temp)
                endInds(iq) = temp;
            end
            prog.update(iq);
        end    
        prog.finish();
        
%         startInds = ceil((startTimes - nsxData(insx).timeStart) / nsxData(insx).timeSamplePeriod) + 1;
%         endInds =  floor((endTimes - nsxData(insx).timeStart) / nsxData(insx).timeSamplePeriod) + 1;

        % first remove trials outside of the boundaries of .time
        lastValidTrial = find(startInds < length(nsxData(insx).time), 1, 'last');
        % and then move the last trial's end to lie within .time as well
        if length(nsxData(insx).time) - endInds(lastValidTrial) < 0
            endInds(lastValidTrial) = length(nsxData(insx).time);
        end

        % move end inds to be at least 1 before the subsequent startInds
        for iR = 1:numel(startInds)-1
            endInds(iR) = min(startInds(iR+1)-1, endInds(iR));
        end

        startInds(lastValidTrial+1:end) = 0;
        endInds(lastValidTrial+1:end) = 0;

        % number of inds before first trial
        preTrialCols = startInds(1)-1;
        % inds devoted to ith trial
        trialCols = endInds - startInds + 1;
        trialCols(lastValidTrial+1:end) = 0;
        % inds devoted to gap between ith and i+1th trial
        interTrialCols = startInds(2:end) - endInds(1:end-1) - 1;
        interTrialCols(lastValidTrial:end) = 0;

        % inds after last trial
        if lastValidTrial == length(startInds)
            postTrialCols = length(nsxData(insx).time) - endInds(end);
        else
            postTrialCols = 0;
        end

        % now we want preTrialCols, then interleaved trialCols and [interTrialCols postTrialCols]
        colsCount = [preTrialCols reshape([trialCols'; [interTrialCols' postTrialCols]], 1, [])];

        % use these counts to split the giant data matrix into a cell
        debug('Splitting analog data into segments (may take some time)\n');
        timeCell{insx} = mat2cell(makerow(nsxData(insx).time), 1, colsCount);
        dataCell{insx} = mat2cell(nsxData(insx).data, size(nsxData(insx).data, 1), colsCount);
        indexIntoCellForEachQ = 2:2:length(colsCount);
    end

    prog = ProgressBar(length(Q), 'Adding segmented analog data');
    for iq = 1:length(Q)
        prog.update(iq);

        % precompute the segmented data for this trial
        %for insx = 1:length(nsxData)
        %    timeInds = nsxData(insx).time >= Q(iq).CerebusInfo.startTime & ...
        %               nsxData(insx).time <= Q(iq).CerebusInfo.endTime;
        %    data{insx} = nsxData(insx).data(:, timeInds);
        %    timeOffset{insx} = nsxData(insx).time(timeInds) - Q(iq).CerebusInfo.startTime;
        %end

        % look over each analog channel/channel set requested
        for ia = 1:length(analogLookup)
            lk = analogLookup(ia);
            insx = lk.nsxIndex;
            % check for time field in this channel group, add if missing
            if ~isfield(Q(iq), lk.groupName) || ~isfield(Q(iq).(lk.groupName), 'time') || isempty(Q(iq).(lk.groupName).time)
                Q(iq).(lk.groupName).time = timeCell{insx}{indexIntoCellForEachQ(iq)} - startTimes(iq);
            end

            % check for scaleFn and scaleLims field in this channel group, add if missing
            if ~isfield(Q(iq), lk.groupName) || ~isfield(Q(iq).(lk.groupName), 'scaleFn') || isempty(Q(iq).(lk.groupName).scaleFn)
                Q(iq).(lk.groupName).scaleFn = lk.scaleFn;
                Q(iq).(lk.groupName).scaleLims = lk.scaleLims;
            end

            % mark with sampling rate
            Q(iq).(lk.groupName).samplingFreq = nsxData(insx).samplingHz;

            % figure out how to assign the data (single or multiple channel assign with lookup?)
            if(lk.single)
                % single channel assign into group
                dataName = lk.name;
            else
                % multiple channels: assign into matrix, include lookup table
                dataName = 'data';
                Q(iq).(lk.groupName).lookup = lk.lookup;
            end

            % actually assign the data for this channel or channel set
            Q(iq).(lk.groupName).(dataName) = dataCell{insx}{indexIntoCellForEachQ(iq)}(lk.chInd, :);
        end
    end
    prog.finish();
end

