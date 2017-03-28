function [tMin, tMax] = getValidTimeExtents(time, data)
% for cells time and data, return the min and max time for which
% data is present on each trial. If data is a a matrix, the first dim must
% be over trials. If data is a cell over trials, time runs over dim 1

[tMin, tMax] = deal(nan(size(time)));

if iscell(data)
    % different time vector for each
    for i = 1:numel(time)
        if ~isempty(data{i}) && ~isempty(time{i})
            mask = ~all(isnan(data{i}), 2);
            if any(mask)
                tMin(i) = min(time{i}(mask), [], 'omitnan');
                tMax(i) = max(time{i}(mask), [], 'omitnan');
            end
        end
    end
    
else
    % assume matrix data and vector time
    assert(ismatrix(data) && isvector(time));
    
    for i = 1:numel(time)
        mask = ~isnan(data(i, :));
        if any(mask)
            tMin(i) = time(find(mask, 1, 'first'));
            tMax(i) = time(find(mask, 1, 'last'));
        end
    end
end
    
