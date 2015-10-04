function [time, data, tMask] = fixNonmonotonicTimeseries(time, data)
% [time, data, tMask] = fixNonmonotonicTimeseries(time, data)
%    time and data can be either vectors or cell arrays of vectors

    if nargin < 2
        data = [];
    end
    
    if iscell(time)
        if isempty(data)
            % time only, as cell
            [time, tMask] = cellfun(@fixSingleTime, time, 'UniformOutput', false);
        else
            % time and data as cell
            [time, data, tMask] = cellfun(@fixSingle, time, data, 'UniformOutput', false);
        end
    else
        if isempty(data)
            % time only as vector
            [time, tMask] = fixSingleTime(time);
        else
            % time and data as vectors
            [time, data, tMask] = fixSingle(time, data);
        end
    end
    
end
   
function [t, tMask] = fixSingleTime(t)
    if isempty(t)
        t = [];
        tMask = [];
    else
        diffT = diff(t);
        stuck = find(diffT(1:end-1) == 0 & diffT(2:end) == 2);
        t(stuck+1) = t(stuck+1) + 1;
        skip = find(diffT(1:end-1) == 2 & diffT(2:end) == 0);
        t(skip+1) = t(skip+1) - 1;

        tMask = [true; makecol(diff(t)>0)];
        t = t(tMask);
    end
end

function [t, d, tMask] = fixSingle(t, d)
    % d is data, t is time, both must be vectors
    [t, tMask] = fixSingleTime(t);
    d = d(tMask) ;
end
