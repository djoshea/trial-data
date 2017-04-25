function [y, ty] = resamplePadEdges(x, tx, txReference, timeDeltaX, timeDeltaY, interpMethod, binAlignmentMode, uniformlySampled)
    % time runs along first dim of x
    % tx is timestamps associated with x, txReference is a timepoint that
    % will line up with the time sampling in y
    
    Q = ceil(timeDeltaY / timeDeltaX);

    % resample uses a filter length of 20 inside, with some zero edges added in
    P = 24*Q;
    tpre = tx(1) + (-P:-1)' .* timeDeltaX;
    isInt = @(x) ceil(x) == x;
    
    % do some offsetting so that resampling respects the binAlignmentMode
    switch binAlignmentMode
        case BinAlignmentMode.Causal
            addToTy = timeDeltaY/2;
        case BinAlignmentMode.Acausal
            addToTy = -timeDeltaY/2;
        case BinAlignmentMode.Centered
            addToTy = 0;
        otherwise
            error('Unknown binAlignmentMode');
    end
    
    idxFirst = find(isInt((tpre-txReference+addToTy) ./ timeDeltaY), 1);
    tpre = tpre(idxFirst:end);
    P = numel(tpre);
    tpost = tx(end) + (1:P)' .* timeDeltaX;
    tx = [tpre; tx; tpost];
    x = padarray(x, P, 'replicate', 'both');
    
    castDouble = ~isa(x, 'double');
    if castDouble
        clsX = class(x);
        x = double(x);
    end

    % make into 2d matrix
    szX = size(x);
    x = x(:, :);
    colMask = ~all(isnan(x), 1);
    x = x(:, colMask);
    
    if uniformlySampled
        [P, Q] = rat(timeDeltaX / timeDeltaY);
        [y] = resample(x, P, Q);
        ty = tx(1) + (0:timeDeltaY:(size(y, 1)-1)*timeDeltaY)';
    else
        [y,ty] = resample(x, tx, 1./timeDeltaY, interpMethod);
    end
    ty = ty+addToTy;
    
    if castDouble
        y = cast(y, clsX);
    end
    y = TensorUtils.inflateMaskedTensor(y, 2, colMask);
    szY = [size(y, 1), szX(2:end)];
    y = reshape(y, szY);
    
end