function rgb = computeWithAlpha(rgb, alpha, args)
    % alpha-blends color with background in linear sRGB space
    arguments
        rgb double
        alpha (1, 1) {mustBeInRange(alpha, 0, 1)}
        args.background (1, 3) = [1 1 1];
    end

    was_permuted = false;
    if size(rgb,3) == 1 && size(rgb,2) == 3
        rgb = permute(rgb,[1 3 2]);
        was_permuted = true;
    end

%     lin = rgb2lin(rgb);
    lin = rgb;
    lin = alpha*lin + (1-alpha)*shiftdim(args.background, -1);
    rgb = lin;
%     rgb = lin2rgb(lin);
    
    if was_permuted
        rgb = permute(rgb,[1 3 2]);
    end
end