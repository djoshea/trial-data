function [hl, hs] = errorshade(x, ym, ye, color, varargin)
% [ha] = shadeYInterval(x, y1, y2, varargin)
% shadeYInterval draws two lines on a plot and shades the area between those
% lines. ParamValues are same as for fill command (e.g. FaceColor,
% EdgeColor)
%

    p = inputParser();
    p.addParamValue('lineArgs', {}, @iscell);
    p.addParamValue('shadeArgs', {}, @iscell);
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('alpha', 1, @isscalar);
    p.parse(varargin{:});

    axh = p.Results.axh;
    
    y1 = ym - ye;
    y2 = ym + ye;
    
    if all(isnan(y1) | isnan(y2))
        hl = NaN;
        hs = NaN;
        return;
    end

    % plot the shaded area
    x = makerow(x);
    y1 = makerow(y1);
    y2 = makerow(y2);

    % desaturate the color for shading
    shadeColor = 1 - (1-color)*0.5;
    
    % need to split the vecs
    nanMask = isnan(y1) | isnan(y2);
    offset = 1;
    while(offset < numel(x))
        % find next non-nan sample
        newOffset = find(~nanMask(offset:end), 1, 'first');
        if isempty(newOffset)
            break;
        end
        
        offset = offset+newOffset-1;
        nextNaN = find(nanMask(offset:end), 1, 'first');
        if isempty(nextNaN)
            regionEnd = numel(x);
        else
            regionEnd = nextNaN+offset - 2;
        end
        
        regionStart = offset;
        mask = regionStart:regionEnd;
        
        [hs] = shadeSimple(axh, x(mask), y1(mask), y2(mask), 'FaceColor', shadeColor, ...
            'alpha', p.Results.alpha, p.Results.shadeArgs{:});
       
        offset = regionEnd + 1;
    end

    hold(axh, 'on');
    if p.Results.alpha < 1
        % use patchline for drawing translucent lines
        hl = TrialDataUtilities.Plotting.patchline(x, ym, ...
           'EdgeColor', color, 'EdgeAlpha', p.Results.alpha, ...
           p.Results.lineArgs{:});
    else
        % use plot for opaque lines
        hl = plot(x, ym, 'Color', color, 'Parent', axh, p.Results.lineArgs{:});
    end
    
end

function [ha] = shadeSimple(axh, x, y1, y2, varargin)

p = inputParser();
p.addParamValue('FaceColor', [0.8 0.8 1], @(x) true);
p.addParamValue('EdgeColor', 'none', @(x) true);
p.addParamValue('alpha', 1, @isscalar);
p.KeepUnmatched = true;
p.parse(varargin{:});

ha = fill([x, fliplr(x)], [y1, fliplr(y2)], p.Results.FaceColor, ...
    'EdgeColor', p.Results.EdgeColor, 'Parent', axh, 'FaceAlpha', p.Results.alpha, ...
    p.Unmatched);

% hide shading from legend
set(get(get(ha, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'Layer', 'top')

end
