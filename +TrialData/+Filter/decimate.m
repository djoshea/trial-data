function [td, by] = decimate(td, channelName, varargin)
    % using data only within the the current alignment window,
    % compute the iti histogram for the spike unit
    assert(isa(td, 'TrialData'));
   
    p = inputParser();
    p.addParameter('by', [], @isscalar);
    p.addParameter('toHz', [], @isscalar);
    p.parse(varargin{:});
    
    prog = ProgressBar(td.nTrials, 'Decimating channel %s', channelName);
    
    by = [];
    if ~isempty(p.Results.by)
        by = p.Results.by;
    elseif ~isempty(p.Results.toHz)
        currentHz = td.timeUnitsPerSecond / td.getAnalogTimeDelta(channelName);
        by = round(currentHz / p.Results.toHz);
    end
        
    [data, times] = td.getAnalog(channelName);
    
    for i = 1:td.nTrials
        prog.update(i);
        if numel(data{i}) < 24
            continue;
        end
        data{i} = decimate(data{i}, by);
        times{i} = times{i}(1:by:end);
    end
    prog.finish();
    
    td = td.setAnalog(channelName, data, times);
end