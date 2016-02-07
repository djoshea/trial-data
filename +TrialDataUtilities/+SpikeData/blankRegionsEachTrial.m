function [data, masks] = blankRegionsEachTrial(data, times, blankingRegions, varargin)
% data is nTrials x 1 cell of data whose first dimension corresponds to
% each spike, like spike times or waveforms
%
% times is nTrials x 1 cell of spike times
%
% blanking regions is nTrials x 1 cell with nRegions x 2 matrices of time
% points inside
%
% remove all spikes within the indicated time windows

    p = inputParser();
    p.addParameter('mask', truevec(numel(times)), @islogical);
    p.parse(varargin{:});
    
    nTrials = numel(times);
    masks = cellvec(nTrials);
    trialMask = p.Results.mask;
    for iT = 1:nTrials
        masks{iT} = truevec(numel(times{iT}));
        if ~trialMask(iT), continue; end
        
        blank = blankingRegions{iT};
        if isempty(blank), continue; end
        assert(size(blank, 2) == 2)
        for iR = 1:size(blank, 1)
            if any(isnan(blank(iR, :))) || isempty(masks{iT}), continue, end
            masks{iT} = masks{iT} & (times{iT} < blank(iR, 1) | times{iT} > blank(iR, 2));
        end
        
        data{iT} = data{iT}(masks{iT}, :, :, :, :);
    end
end