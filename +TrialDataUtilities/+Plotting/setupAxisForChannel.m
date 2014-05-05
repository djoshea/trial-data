function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('which', 'y', @(str) ismember(str, {'x', 'y', 'z'}));
    p.addParamValue('useAutoAxis', true, @islogical);
    p.addParamValue('style', 'tickBridge', @ischar);
    p.parse(varargin{:});
    
    which = p.Results.which;
    axh = p.Results.axh;
    
    if strcmp(p.Results.style, 'scaleBar')
        useScaleBar = true;
        label = channelDescriptor.name;
    elseif strcmp(p.Results.style, 'tickBridge')
        useScaleBar = false;
        label = channelDescriptor.getAxisLabelPrimary();
    else
        error('Invalid axis style %s', p.Results.style);
    end
    
    if p.Results.useAutoAxis
        au = AutoAxis(axh);
    end
    switch which
        case 'x'
            xlabel(axh, label);
            set(get(axh, 'XLabel'), 'Visible', 'on');
            if p.Results.useAutoAxis
                if useScaleBar
                    au.xUnits = channelDescriptor.unitsPrimary;
                    au.addAutoScaleBarX();
                else
                	au.addAutoAxisX();
                end
            end
        case 'y'
            ylabel(axh, label);
            set(get(axh, 'XLabel'), 'Visible', 'on');
            if p.Results.useAutoAxis
                if useScaleBar
                    au.yUnits = channelDescriptor.unitsPrimary;
                    au.addAutoScaleBarY();
                else
                    au.addAutoAxisY();
                end
            end
        case 'z'
            zlabel(axh, label);
    end

    if p.Results.useAutoAxis
        %au.update();
    end
end
