function [hbars, hbaseline] = stem_indep_baseline(x, y, baselineValue, args)
    
arguments
        x (:, 1)
        y (:, 1)
        baselineValue = 0;
        args.width = [];
        args.color = 'k';
        args.drawBaseline (1, 1) logical = true;
        args.baselineArgs = struct();
        args.barArgs = struct();
    end
    
    baselineArgs = namedargs2cell(args.baselineArgs);
    barArgs = namedargs2cell(args.barArgs);
   
    if isempty(args.width)
        args.width = (x(2) - x(1)) * 0.8;
    end
    
    xb_lo = min(x, [], 'omitnan') - args.width;
    xb_hi = max(x, [], 'omitnan') + args.width;
    
    N = numel(x);
    x1 = x' - args.width/2;
    x2 = x' + args.width/2;
    y1 = repmat(baselineValue, 1, N);
    y2 = y';
    X = [x1; x1; x2; x2];
    Y = [y1; y2; y2; y1];
    
    washolding = ishold;
    
    hbars = fill(X, Y, args.color, 'EdgeColor', args.color, barArgs{:});
    
    if args.drawBaseline
        hold on;
        hbaseline = plot([xb_lo, xb_hi], [baselineValue baselineValue], '-', 'Color', args.color, baselineArgs{:});
    end
    if ~washolding
        hold off;
    end
    
    
end