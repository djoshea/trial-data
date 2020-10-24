function hs = errorshadeInterval(x, lo, hi, color, varargin)
% hs = errorshadeInterval(x, lo, hi, color, varargin)

    p = inputParser();
    p.addParameter('errorColor', [], @(x) true);
    p.addParameter('shadeArgs', {}, @iscell);
    p.addParameter('axh', [], @(x) true);
    p.addParameter('alpha', 1, @isscalar);
    p.addParameter('z', 0, @isscalar); % used for visual stacking on 2-d plots
    p.addParameter('Clipping', true, @islogical);
    p.parse(varargin{:}); 
    
    color = TrialDataUtilities.Plotting.convertColorsToMatrix(color);
    
    z = p.Results.z;

    if isempty(p.Results.axh)
        axh = newplot;
    else
        axh = p.Results.axh;
    end
   
    y1 = lo;
    y2 = hi;
    hs = [];
    
    if all(isnan(y1) | isnan(y2))   
        return;
    end

    % plot the shaded area
    x = makerow(x);
    y1 = makerow(y1);
    y2 = makerow(y2);

    % desaturate the color for shading if not translucent
    if isempty(p.Results.errorColor)
        if p.Results.alpha < 1
            shadeColor = color;
        else
            % move half the distance to white
            shadeColor = 1 - (1-color)*0.5;
        end
    else
        shadeColor = p.Results.errorColor;
    end
    
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
        
        hs = shadeSimple(axh, x(mask), y1(mask), y2(mask), z, 'FaceColor', shadeColor, ...
            'alpha', p.Results.alpha, 'Clipping', p.Results.Clipping, p.Results.shadeArgs{:});
        TrialDataUtilities.Plotting.hideInLegend(hs);
        
        offset = regionEnd + 1;
    end
end

function [ha] = shadeSimple(axh, x, y1, y2, z, varargin)

p = inputParser();
p.addParameter('FaceColor', [0.8 0.8 1], @(x) true);
p.addParameter('EdgeColor', 'none', @(x) true);
p.addParameter('alpha', 1, @isscalar);
p.addParameter('Clipping', true, @islogical);
p.KeepUnmatched = false;
p.parse(varargin{:});

faceColor = TrialDataUtilities.Plotting.convertColorsToMatrix(p.Results.FaceColor);
edgeColor = TrialDataUtilities.Plotting.convertColorsToMatrix(p.Results.EdgeColor);

xv = [x, fliplr(x)];
yv = [y1, fliplr(y2)];
zv = z * ones(size(xv));

ha = patch(xv, yv, zv, 'k', 'Parent', axh);
set(ha, 'FaceColor', faceColor, ...
    'EdgeColor', edgeColor, 'Parent', axh, 'FaceAlpha', p.Results.alpha, 'Clipping', p.Results.Clipping);

% hide shading from legend
set(get(get(ha, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

% set(axh, 'Layer', 'top')

end
