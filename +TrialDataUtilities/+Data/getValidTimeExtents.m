function [tMin, tMax, indMin, indMax] = getValidTimeExtents(time, data, varargin)
% for cells time and data, return the min and max time for which
% data is present on each trial. If data is a a matrix, the first dim must
% be over trials. If data is a cell over trials, time runs over dim 1
%
% if data is nTrials x nChannels cell, then origDelta is nChannels x 1 of 
% the time between samples for each channel
% if data is a matrix, it is nTrials x nTime
%

p = inputParser();
% these are used as a secondary guard to truncate data within tMin :
% tMax, when the input data includes padded edges to facilitate
% resampling
p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
p.parse(varargin{:});

if iscell(data)
    [tMin, tMax, indMin, indMax] = deal(nan(size(data)));
    
    tMinExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMinExcludingPadding, size(time));
    tMaxExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMaxExcludingPadding, size(time));

    % different time vector for each
    for i = 1:numel(data)
        if ~isempty(data{i}) && ~isempty(time{i})
            % need to make data{i} into matrix to detect if all NaNs
            mask = ~all(isnan(data{i}(:, :)), 2) & time{i} >= tMinExcludingPadding(i) & time{i} <= tMaxExcludingPadding(i);
            if any(mask)
                [tMin(i), indMin(i)] = min(time{i}(mask), [], 'omitnan');
                [tMax(i), indMax(i)] = max(time{i}(mask), [], 'omitnan');
            end
        end
    end
else
    [tMin, tMax, indMin, indMax] = deal(nan(size(data, 1), 1));
    
    tMinExcludingPadding = p.Results.tMinExcludingPadding;
    tMaxExcludingPadding = p.Results.tMaxExcludingPadding;

    % assume matrix data and vector time
    assert(ismatrix(data) && isvector(time));
    
    time = makecol(time);
    
    mask = ~all(isnan(data), 1);
    
    for i = 1:size(data, 1)
        mask = ~isnan(data(i, :))' & time >= tMinExcludingPadding & time <= tMaxExcludingPadding;
        if any(mask)
            indMin(i) = find(mask, 1, 'first');
            tMin(i) = time(indMin(i));
            indMax(i) = find(mask, 1, 'last');
            tMax(i) = time(indMax(i));
        end
    end
end
