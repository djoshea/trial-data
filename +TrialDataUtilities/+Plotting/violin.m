function [h, pts_lims] = violin(loc, vals, args, patchArgs)
    arguments
        loc (1, 1)
        vals (:, 1)
    
        args.orientation (1, 1) string = "vertical";
        args.support (1, 1) string = ""
        args.bandwidth (1, 1) = NaN;
        args.trimQuantiles (1, 2) = nan(1, 2);
        args.width (1, 1) = 1;
    
        args.FaceColor = [0.7 0.7 0.7];
        args.EdgeColor = [0.2 0.2 0.2];
    
        patchArgs.?matlab.graphics.primitive.Patch
    end
    
    ks_args = struct();
    if args.support == "minmax"
        ks_args.support = [min(vals) - eps(min(vals)), max(vals) + eps(max(vals))];
    elseif args.support ~= ""
        ks_args.support = char(args.support);
    end
    
    if ~isnan(args.bandwidth)
        ks_args.bandwidth = args.bandwidth;
    end

    ks_args_c = namedargs2cell(ks_args);
    [f, pts] = ksdensity(vals, ks_args_c{:});
    
    if any(~isnan(args.trimQuantiles))
        qtiles = quantile(vals, args.trimQuantiles);
        mask = pts < qtiles(1) | pts > qtiles(2);
        f(mask) = [];
        pts(mask) = [];
    end
    
    f = f';
    pts = pts';
    
    pts_lims = [min(pts), max(pts)];

    f = f / max(f) * args.width/2; %normalize
    XX = [f + loc; flipud(loc - f)];
    YY = [pts; flipud(pts)];
    if args.orientation == "horizontal"
        [XX, YY] = swap(XX, YY);
    end
    
    patchArgsC = namedargs2cell(patchArgs);
    h = fill(XX, YY, args.FaceColor, 'EdgeColor', args.EdgeColor, patchArgsC{:});

end

function [b, a] = swap(a, b)
end