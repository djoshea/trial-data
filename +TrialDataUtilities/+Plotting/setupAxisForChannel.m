function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('which', 'y', @(str) ismember(str, {'x', 'y', 'z'}));
    p.addParamValue('useAutoAxis', true, @islogical);
    p.addParamValue('style', 'tickBridge', @ischar);
    p.parse(varargin{:});
    
    which = p.Results.which;
    axh = p.Results.axh;
    
    if isa(channelDescriptor, 'ChannelDescriptor')
        label = channelDescriptor.getAxisLabelPrimary();
        units = channelDescriptor.unitsPrimary;
    elseif isstruct(channelDescriptor)
        name = channelDescriptor.name;
        units = channelDescriptor.units;
        if isempty(units)
            label = name;
        else
            label = sprintf('%s (%s)', name, units);
        end
    end  
    
    useAutoAxis = p.Results.useAutoAxis;
    
    if strcmp(p.Results.style, 'scaleBar')
        useLabel = false;
        useScaleBar = true;
        useAutoAxis = true;
    elseif strcmp(p.Results.style, 'tickBridge')
        useLabel = true;
        useScaleBar = false;
    elseif strcmp(p.Results.style, 'none')
        useLabel = false;
        useScaleBar = false;
        useAutoAxis = false;
    else
        error('Invalid axis style %s', p.Results.style);
    end
    
    if strcmp(which, 'z')
        useAutoAxis = false;
    end
    
    if ~useAutoAxis || useLabel
        switch which
            case 'x'
                xlabel(axh, label);
                set(get(axh, 'XLabel'), 'Visible', 'on');
            case 'y'
                ylabel(axh, label);
                set(get(axh, 'YLabel'), 'Visible', 'on');
            case 'z'
                zlabel(axh, label);
                set(get(axh, 'ZLabel'), 'Visible', 'on');s
        end
    end
    
    if ~useAutoAxis
        return;
    end
    
    au = AutoAxis(axh);
    switch which
        case 'x'
            au.xUnits = units;
            if useScaleBar
                au.addAutoScaleBarX();
            else
                au.addAutoAxisX();
            end
        case 'y'
            au.yUnits = units;
            if useScaleBar
                au.addAutoScaleBarY();
            else
                au.addAutoAxisY();
            end
    end
    
    au.update();
end
