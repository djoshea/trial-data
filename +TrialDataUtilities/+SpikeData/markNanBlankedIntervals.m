function [data, time] = markNanBlankedIntervals(blankingIntervals, data, time, varargin)

p = inputParser();
p.addParameter('padding', [0 0], @isvector);
p.parse(varargin{:});

pre = p.Results.padding(1);
post = p.Results.padding(2);

if iscell(data)
    nTrials = numel(data);
    
    for iT = 1:nTrials
        blank = blankingIntervals{iT};
        
        % padding means that data at time t requires that there be no blanking from t-pre to t+post 
        for iR = 1:size(blank, 1)
            mask = time{iR} >= blank(iR, 1) - pre & time{iR} <= blank(iR, 2) + post;
            data{iT}(mask, :) = NaN;
        end
    end
else
    nTrials = size(data, 1);
    
    if isempty(data)
        return;
    end
    
    for iT = 1:nTrials
        blank = blankingIntervals{iT};
        
        % padding means that data at time t requires that there be no blanking from t-pre to t+post 
        for iR = 1:size(blank, 1)
            mask = time >= blank(iR, 1) - pre & time <= blank(iR, 2) + post;
            data(iT, mask, :) = NaN;
        end
    end
end

end