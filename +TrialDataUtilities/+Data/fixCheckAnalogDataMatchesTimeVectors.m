function [data, okayMask, messages] = fixCheckAnalogDataMatchesTimeVectors(data, dataFields, timeFields)
    % this is very similar to the function in TrialData but is externalized
    % to allow TrialDataInterfaces to call it directly to save time and
    % unnecessary warnings. Note that this EXPLICITLY assumes that no
    % channels share a time field, unless they are already in a group (with
    % one entry in dataField corresponding to it).

    nTotal = numel(dataFields);
    nTrials = numel(data);
    messages = cell(nTotal, 1);
    okayMask = true(nTotal, 1);
    
    prog = ProgressBar(nTotal, 'Checking sample count vs. times for analog channels');
    for iA = 1:nTotal
        prog.update(iA, 'Checking sample count vs. times for %s', dataFields{iA});

        dataField = dataFields{iA};
        timeField = timeFields{iA};

        okay = truevec(nTrials);
        [transpose_time, transpose_data] = falsevec(nTrials);

        for iT = 1:nTrials
            sz = size(data(iT).(timeField));
            
            if sz(2) > 1 && sz(1) == 1
                transpose_time(iT) = true;
                nTime = size(data(iT).(timeField), 2);
            else
                nTime = size(data(iT).(timeField), 1);
            end

            sz = size(data(iT).(dataField));
            if sz(2) > 1 && sz(1) == 1
                transpose_data(iT) = true;
                nData = size(data(iT).(dataField), 2);
                nChannelsThis = 1;
            else
                nData = size(data(iT).(dataField), 1);
                nChannelsThis = sz(2);
            end
            
            if prod(sz) > 0
                memoryClass = class(data(iT).(dataField));
            end

            okay(iT) = nTime == nData;
        end

        % do the actual transposing down here for efficiency
        if any(transpose_time)
            tr_time = cellfun(@transpose, {data(transpose_time).(timeField)}, 'UniformOutput', false);
            data = TrialDataUtilities.Data.assignIntoStructArray(data, timeField, tr_time, transpose_time);
        end

        if any(transpose_data)
            tr_data = cellfun(@transpose, {data(transpose_data).(dataField)}, 'UniformOutput', false);
            data = TrialDataUtilities.Data.assignIntoStructArray(data, timeField, tr_data, transpose_data);
        end

        if any(~okay)
            okayMask(iA) = false;
            messages{iA} = strjoin([messages{iA}, sprintf("%d trials have differing number of data samples in %s as timestamps in %s. Fixing by clearing data and time fields.", ...
                nnz(~okay), dataField, timeField)], "; ");

            if strcmp(memoryClass, 'char') %#ok<ISCHR>
                assert(nChannelsThis == 1, 'Enum channel groups not yet supported?');
                emptyVal = '';
            elseif strcmp(memoryClass, 'logical') %#ok<ISLOG>
                emptyVal = false(0, nChannelsThis);
            else
                emptyVal = zeros(0, nChannelsThis, memoryClass);
            end

            data = TrialDataUtilities.Data.assignIntoStructArray(data, timeField, zeros(0, 1), ~okay);
            data = TrialDataUtilities.Data.assignIntoStructArray(data, dataField, emptyVal, ~okay);
        end

        timeData = {data.(timeField)}';
        emptyMask = cellfun(@(t) numel(t) < 2, timeData);
        timeDelta = nanmedian(cellfun(@(x) nanmedian(diff(x)), timeData(~emptyMask)));

        % check sorted and no duplicates
        resort = truevec(nTrials);
        for iT = 1:nTrials
            timeThis = timeData{iT};
            timeThis = TrialDataUtilities.Data.removeSmallTimeErrors(timeThis, timeDelta, 0);
            resort(iT) = ~issorted(timeThis) || numel(timeThis) > numel(unique(timeThis));
        end
        if any(resort)
            messages{iA} = strjoin([messages{iA}, sprintf("%d trials have duplicate or non-monotonically increasing timestamps for %s. Fixing by deleting non-montonic samples.", nnz(resort), dataField)], "; ");
            okayMask(iA) = false;
            [timeSorted, timeSortIdx] = cellfun(@getTimeSorted, timeData(resort), 'UniformOutput', false);
            dataSorted = cellfun(@sortData, {data(resort).(dataField)}', timeSortIdx, 'UniformOutput', false);
            data = TrialDataUtilities.Data.assignIntoStructArray(data, dataField, dataSorted, resort);
            data = TrialDataUtilities.Data.assignIntoStructArray(data, timeField, timeSorted, resort);
        end

    end
    prog.finish();

    function [time, sortIdx] = getTimeSorted(time)
        time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, 0);
        [time, sortIdx] = unique(time, 'last');
    end

    function data = sortData(data, sortIdx)
        data = data(sortIdx, :, :, :, :, :, :, :, :);
    end

%             function [time, data] = resortTimeDedup(time, data)
%                 time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, 0);
%                 [time, idx] = unique(time, 'last');
%                 data = data(idx, :, :, :, :);
%             end
end