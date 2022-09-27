function rgb = adjust_oklch(rgb, args)

    arguments
        rgb double
        args.shiftL = 0;
        args.gainL = 1;
        args.minL = 0;
        args.maxL = 1;

        args.shiftC = 0;
        args.gainC = 1;
        args.minC = 0;
        args.maxC = 0.35;
        
        args.shiftH = 0;
        args.minH = 0;
        args.maxH = 360;
    end
    
    if ~isnumeric(rgb)
        return;
    end
            
    % using the oklab color space, reduces the chroma + lighby dividing by factor
    was_permuted = false;
    if size(rgb,3) == 1 && size(rgb,2) == 3 || size(rgb,2) == 4
        rgb = permute(rgb,[1 3 2]);
        was_permuted = true;
    end

    if size(rgb, 3) == 4
        % has alpha channel
        alpha = rgb(:, :, 4);
        rgb = rgb(:, :, 1:3);
    else
        alpha = [];
    end

    oklch = TrialDataUtilities.Color.rgb2oklch(rgb);
    clamp = @(x, lo, hi) max(min(x, hi), lo);
    oklch(:, :, 1) = clamp(oklch(:, :, 1).*args.gainL + args.shiftL, args.minL, args.maxL);
    oklch(:, :, 2) = clamp(oklch(:, :, 2).*args.gainC + args.shiftC, args.minC, args.maxC);
    oklch(:, :, 3) = clamp(mod(oklch(:, :, 3) + args.shiftH, 360),   args.minH, args.maxH);

    rgb = TrialDataUtilities.Color.oklch2rgb(oklch);
    
    if ~isempty(alpha)
        rgb = cat(3, rgb, alpha);
    end
    if was_permuted
        rgb = permute(rgb, [1 3 2]);
    end

end