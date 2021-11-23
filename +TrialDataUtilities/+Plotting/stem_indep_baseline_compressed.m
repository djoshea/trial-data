function [hbars, hbaseline] = stem_indep_baseline_compressed(x, y, baselineValue, args)
    % like stem_indep_baseline but bars are touching (much faster)
arguments
        x (:, 1)
        y (:, 1)
        baselineValue = 0;
        args.width = []; % only affects first and last bar
        args.color = 'k';
        args.drawBaseline (1, 1) logical = true;
        args.baselineArgs = struct();
        args.barArgs = struct();
    end
    
    baselineArgs = namedargs2cell(args.baselineArgs);
    barArgs = namedargs2cell(args.barArgs);
   
    if isempty(args.width)
        args.width = x(2) - x(1);
    end
    
    xb_lo = min(x, [], 'omitnan') - args.width;
    xb_hi = max(x, [], 'omitnan') + args.width;
    
    % row vectors 1 x N
    N = numel(x);
    x1 = [x(1) - args.width/2; 0.5 * (x(1:end-1) + x(2:end))]'; % left edges
    x2 = [0.5 * (x(1:end-1) + x(2:end)); x(end) + args.width/2]';
    y1 = repmat(baselineValue, 1, N);
    y2 = y';
    
    ymin = min(y1, y2);
    ymax = max(y1, y2);
    
    X = cat(2, [x1; x2], fliplr([x2; x1]));
    Y = cat(2, [ymax; ymax], fliplr([ymin; ymin]));
    
    X = cat(1, X(:), x1(1));
    Y = cat(1, Y(:), ymax(1));
    
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