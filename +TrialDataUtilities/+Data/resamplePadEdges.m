function [y, ty] = resamplePadEdges(x, tx, ty, binAlignmentMode, interpolateMode)
    % time runs along first dim of x
    % tx is timestamps associated with x, txReference is a timepoint that
    % will line up with the time sampling in y
    
    % our choice of first timestamp is important, as it sets the effective
    % times that will be sampled by the Y returned by resample, i.e. ty
    % will begin at tpre(1) but have interfal timeDeltaY instead.
    % So we need to take into account the binAlignmentMode here in order to
    % get data consistent with the binAlignmentMode.
    %
    % If we're centering the bins, then no adjustment is needed. But since
    % resampling assumes centered samples, we need to adjust for causal and
    % acausal.
    % Example 1:
    %   incoming is 1 ms, resample to 20 ms, causal bins:
    %   we treat the incoming data (0:100) as causal 1 ms binned, with
    %   centers -0.5:99.5. The outgoing data will be ty = 20:20:100 with centers
    %   10:20:90. So first, we want to interpolate x to get a centered x
    %   with centers at 0:100, which means interpolating x to 0.5:1:100.5.
    %   Then we want to pad the data backwards in time until it aligns with the outgoing
    %   centers (first timepoint is 10-some k * 20). Then we adjust the
    %   time vector 
    % 
    % Example 2
    %   incoming is 20 ms, causal bins, upsample to 1 ms.
    %   incoming data is 0:20:100, centers at -10:20:90. So we interp x to get
    %   x centered at 0:20:100. 
    
    if nargin < 5
        interpolateMode = 'linear';
    end
    
    x = makecol(x);
    tx = double(makecol(tx));
    ty = double(makecol(ty));

    % A trial with no samples in the requested window is a legitimate input, not an error: it
    % simply has no data there. resample() rejects an empty X, so answer directly with all-NaN
    % on the requested output time base -- the same thing the NaN-edge invalidation at the end
    % of this function produces for a trial that is only partially covered.
    %
    % Without this, ONE such trial among many aborts the whole extraction. That mixed case is
    % new: with a stim-locked alignment the window sits a fixed distance from the data edge, so
    % either every trial can supply it or none can, and the all-empty case is handled upstream.
    % With a movement-locked alignment (Estim<Fam>AlignMovePosthocCommon) the distance varies
    % per trial with RT, so a handful of trials come up empty while the rest are fine --
    % measured 2 of 70 on EstimReach P20180609_A baselineLate.
    if isempty(x) || isempty(tx)
        szOut = size(x);
        szOut(1) = numel(ty);
        y = nan(szOut, 'like', x);
        return;
    end

    timeDeltaX = median(diff(tx));
    timeDeltaY = median(diff(ty));

    % resample only works with doubles, so upcast it here and convert it
    % back later if needed
    castDouble = ~isa(x, 'double');
    if castDouble
        clsX = class(x);
        x = double(x);
    end
    
    % figure out padding needed
    [P, Q] = rat(timeDeltaX / timeDeltaY);
    maxPQ = max(P, Q);
    
    % resample uses a filter length of 20*max(P,Q) inside, with some zero edges added in
    % we overdo the padding here so that we can truncate it later
    pad = double(24*maxPQ);
    
    % first fill in NaN tails
    [x, idxFirst, idxLast] = TrialDataUtilities.Data.infillNanEdgesWithLastSample(x);
    
    % interpolate the incoming data to compensate for the bin centering
    % x will now have centered bins on tx
    if binAlignmentMode ~= BinAlignmentMode.Centered
        txInterpToCenter = tx - binAlignmentMode.getOffsetToBinCenter(timeDeltaX);
        x = interp1(tx, x, txInterpToCenter,  interpolateMode, 'extrap');
    end
    
    % then replicate-pad edges 
    x = padarray(x, pad, 'replicate', 'both');
    
    tpre = tx(1) + (-pad:-1)' .* timeDeltaX;
    isInt = @(x) ceil(x) == x;
    tpost = tx(end) + (1:pad)' .* timeDeltaX;
    
    timeFirst = nan(size(idxFirst));
    timeLast = nan(size(idxLast));
    mask = ~isnan(idxFirst) & ~isnan(idxLast);
    timeFirst(mask) = tx(idxFirst(mask));
    timeLast(mask) = tx(idxLast(mask));
    
    tx = [tpre; tx; tpost]; 
    
    % make into 2d matrix
    szX = size(x);
    x = x(:, :);
    colMask = ~all(isnan(x), 1);
    x = x(:, colMask);

    % Every channel all-NaN over this window leaves x with zero columns, and resample() rejects
    % an empty X. That is a legitimate input, not an error -- the trial simply has no data
    % here -- and the answer is all-NaN on the output time base, which is exactly what the
    % inflate below produces from a zero-column y. Short-circuit to it rather than resampling
    % nothing.
    %
    % ONE such trial among many otherwise aborts the whole extraction. The mixed case is new:
    % a stim-locked alignment puts the window a fixed distance from the data, so either every
    % trial covers it or none does, and the all-empty case is handled upstream. A
    % movement-locked alignment puts it at a per-trial-varying distance, so a few trials come
    % up empty while the rest are fine -- measured 2 of 70 on EstimReach P20180609_A
    % baselineLate under Estim<Fam>AlignMovePosthocCommon.
    if ~any(colMask)
        y = nan(numel(ty), 0, 'like', x);
        if castDouble
            y = cast(y, clsX);
        end
        y = TensorUtils.inflateMaskedTensor(y, 2, colMask);
        y = reshape(y, [size(y, 1), szX(2:end)]);
        return;
    end

    [P, Q] = rat(timeDeltaX / timeDeltaY, timeDeltaY*0.001);
    timeDeltaY = timeDeltaX * Q / P; % this will be the actually realized sampling rate
    [y] = resample(x, P, Q);
    
    % this is the effective time vector of the data in Y after resampling
    % to sampling rate timeDelta
    tyRaw = tx(1) + (0:timeDeltaY:(size(y, 1)-1)*timeDeltaY)';
    
    if binAlignmentMode ~= BinAlignmentMode.Centered
        % interpolate the resampled data to the correct bin alignment,
        % currently y is centered on ty
        tyInterpFromCenter = tyRaw - binAlignmentMode.getOffsetToBinCenter(timeDeltaY);
        y = interp1(tyInterpFromCenter, y, ty,  'linear');
    else
        % then we find the alignment between tyRaw and ty to grab the right
        % window
        [delta, idxFirst] = min(abs(tyRaw - ty(1)));
%         assert(delta < timeDeltaY / 4);
        y = y(idxFirst:idxFirst+numel(ty)-1, :);
    end
    
    % invalidate nan edges again based on time vectors
    for c = 1:size(y, 2)
        if ~isnan(timeFirst(c)) && ~isnan(timeLast(c))
            y(ty < timeFirst(c), c) = NaN;
            y(ty > timeLast(c), c) = NaN;
        end
    end
    
    if castDouble
        y = cast(y, clsX);
    end
    y = TensorUtils.inflateMaskedTensor(y, 2, colMask);
    szY = [size(y, 1), szX(2:end)];
    y = reshape(y, szY);
end

%% Old version

% function [y, ty] = resamplePadEdges(x, tx, txReference, timeDeltaX, timeDeltaY, interpMethod, binAlignmentMode, uniformlySampled)
%     % time runs along first dim of x
%     % tx is timestamps associated with x, txReference is a timepoint that
%     % will line up with the time sampling in y
%     
%     Q = TrialDataUtilities.Stats.ceiltol(timeDeltaY / timeDeltaX, min(timeDeltaY, timeDeltaX) / 1000);
% 
%     % resample uses a filter length of 20 inside, with some zero edges added in
%     P = 24*Q;
%     tpre = tx(1) + (-P:-1)' .* timeDeltaX;
%     isInt = @(x) ceil(x) == x;
%     
%     % do some offsetting so that resampling respects the binAlignmentMode
%     switch binAlignmentMode
%         case BinAlignmentMode.Causal
%             addToTy = timeDeltaY/2;
%         case BinAlignmentMode.Acausal
%             addToTy = -timeDeltaY/2;
%         case BinAlignmentMode.Centered
%             addToTy = 0;
%         otherwise
%             error('Unknown binAlignmentMode');
%     end
%     
%     idxFirst = find(isInt((tpre-txReference+addToTy) ./ timeDeltaY), 1);
%     tpre = tpre(idxFirst:end);
%     P = numel(tpre);
%     tpost = tx(end) + (1:P)' .* timeDeltaX;
%     tx = [tpre; tx; tpost];
%     x = padarray(x, P, 'replicate', 'both');
%     
%     castDouble = ~isa(x, 'double');
%     if castDouble
%         clsX = class(x);
%         x = double(x);
%     end
% 
%     % make into 2d matrix
%     szX = size(x);
%     x = x(:, :);
%     colMask = ~all(isnan(x), 1);
%     x = x(:, colMask);
% 
%     if uniformlySampled
%         [P, Q] = rat(timeDeltaX / timeDeltaY);
%         [y] = resample(x, P, Q);
%         ty = tx(1) + (0:timeDeltaY:(size(y, 1)-1)*timeDeltaY)';
%     else
%         [y,ty] = resample(x, tx, 1./timeDeltaY, interpMethod);
%     end
%     ty = ty+addToTy;
%     
%     if castDouble
%         y = cast(y, clsX);
%     end
%     y = TensorUtils.inflateMaskedTensor(y, 2, colMask);
%     szY = [size(y, 1), szX(2:end)];
%     y = reshape(y, szY);
%     
% end