function timeDelta = inferTimeDeltaFromSampleTimes(time, varargin)

p = inputParser();
p.addOptional('data', @(x) isempty(x) || ~isscalar(x) || ~islogical(x));
p.addParameter('ignoreNaNSamples', false, @islogical); % if true, nan samples will be ignored and the time jump will span non-nan samples
p.parse(varargin{:});

% time may be a vector, or cell array of time vectors.
% data may be a matrix whos first dim is over trials or a cell (nTrials x nChannels) whose contents have time on dim 1
data = p.Results.data;
ignoreNaN = p.Results.ignoreNaNSamples && ~isempty(data);

if iscell(time)
    % different time vector for each
    timeDelta = nan(size(time));
    if ignoreNaN
        for i = 1:numel(time)
            if ~isempty(data{i}) && ~isempty(time{i})
                mask = ~all(isnan(data{i}), 2);
                if any(mask)
                    if numel(time{i}(mask)) == 1
                        timeDelta(i) = NaN; % 0 means single sample
                    else
                        timeDelta(i) = median(diff(time{i}(mask)), 1, 'omitnan'); % ignore the mask here
                    end
                end
            end
        end
    else
        for i = 1:numel(time)
            if ~isempty(data{i}) && ~isempty(time{i})
                timeDelta(i) = median(diff(time{i}), 1, 'omitnan'); % ignore the mask here
            end
        end
    end
    
    timeDelta = median(timeDelta, 1, 'omitnan');
    
else
    % single time vector
    if ignoreNaN
        % assume matrix data and vector time
        assert(ismatrix(data) && isvector(time));
        time = makecol(time);
        mask = ~all(isnan(data), 1);
        timeDelta = median(diff(time(mask)), 1, 'omitnan');
    else
        timeDelta = median(diff(time), 1, 'omitnan');
    end
end
