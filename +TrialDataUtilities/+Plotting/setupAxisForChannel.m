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
        useTickBridge = false;
    elseif strcmp(p.Results.style, 'tickBridge')
        useLabel = true;
        useScaleBar = false;
        useTickBridge = true;
    elseif strcmp(p.Results.style, 'tickBridgeScaleBar')
        useLabel = true;
        useScaleBar = true;
        useAutoAxis = true;
        useTickBridge = true;
    elseif strcmp(p.Results.style, 'none')
        useLabel = false;
        useScaleBar = false;
        useAutoAxis = false;
        useTickBridge = false;
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
                set(get(axh, 'XLabel'), 'Visible', 'on', 'Interpreter', 'none');
            case 'y'
                ylabel(axh, label);
                set(get(axh, 'YLabel'), 'Visible', 'on', 'Interpreter', 'none');
            case 'z'
                zlabel(axh, label);
                set(get(axh, 'ZLabel'), 'Visible', 'on', 'Interpreter', 'none');
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
            end
            if useTickBridge
                au.addAutoAxisX();
            end
        case 'y'
            au.yUnits = units;
            if useScaleBar
                au.addAutoScaleBarY();
            end
            if useTickBridge
                au.addAutoAxisY();
            end
    end
    
    au.update();
end
