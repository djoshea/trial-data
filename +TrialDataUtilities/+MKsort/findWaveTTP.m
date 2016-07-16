function [ttp, timeTrough, timePeak] = findWaveTTP(waves, waveTvec)
% based on code written by Matt Kaufman
% ttp is in units of waveTvec

    nW = size(waves, 1);
    [ttp, timeTrough, timePeak] = deal(nan(nW, 1));

    for iW = 1:nW
        wave = waves(iW, :);

        [~, mini] = min(wave);

        % The post peak is where the derivative first goes negative
        diffs = diff(wave);
        postPeak = mini + find(diffs(mini+1:end) < 0, 1);
        if ~isempty(postPeak)
          ttp(iW) = waveTvec(postPeak)-waveTvec(mini);
          timeTrough(iW) = waveTvec(mini);
          timePeak(iW) = waveTvec(postPeak);
        end
    end

end