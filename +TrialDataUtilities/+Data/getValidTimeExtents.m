function [tMin, tMax, timeDelta, indMin, indMax] = getValidTimeExtents(time, data, varargin)
% for cells time and data, return the min and max time for which
% data is present on each trial. If data is a a matrix, the first dim must
% be over trials. If data is a cell over trials, time runs over dim 1
%
% if data is nTrials x nChannels cell, then origDelta is nChannels x 1 of 
% the time between samples for each channel

p = inputParser();
% these are used as a secondary guard to truncate data within tMin :
% tMax, when the input data includes padded edges to facilitate
% resampling
p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
p.parse(varargin{:});

[tMin, tMax, indMin, indMax] = deal(nan(size(time)));

if iscell(data)
    tMinExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMinExcludingPadding, size(time));
    tMaxExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMaxExcludingPadding, size(time));

    % different time vector for each
    timeDelta = nan(size(data));
    for i = 1:numel(data)
        if ~isempty(data{i}) && ~isempty(time{i})
            mask = ~all(isnan(data{i}), 2) & time{i} >= tMinExcludingPadding(i) & time{i} <= tMaxExcludingPadding(i);
            if any(mask)
                [tMin(i), indMin(i)] = min(time{i}(mask), [], 'omitnan');
                [tMax(i), indMax(i)] = max(time{i}(mask), [], 'omitnan');
            end
            if nargout > 2
                timeDelta(i) = nanmedian(diff(time{i}(mask))); % ignore the mask here
            end
        end
    end
    
    timeDelta = nanmedian(timeDelta, 1);
else
    tMinExcludingPadding = p.Results.tMinExcludingPadding;
    tMaxExcludingPadding = p.Results.tMaxExcludingPadding;

    % assume matrix data and vector time
    assert(ismatrix(data) && isvector(time));
    
    for i = 1:numel(time)
        mask = ~isnan(data(i, :)) & time >= tMinExcludingPadding & time <= tMaxExcludingPadding;
        if any(mask)
            indMin(i) = find(mask, 1, 'first');
            tMin(i) = time(indMin(i));
            indMax(i) = find(mask, 1, 'last');
            tMax(i) = time(indMax(i));
        end
    end
    
    mask = ~all(isnan(data), 1);
    if nargout > 2
        timeDelta = nanmedian(diff(time(mask)));
    end
end
