function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('which', 'y', @(str) ismember(str, {'x', 'y', 'z'}));
    p.addParamValue('useAutoAxis', true, @islogical);
    p.addParamValue('style', 'tickBridge', @ischar);
    p.addParameter('label', '', @ischar);
    p.parse(varargin{:});
    
    which = p.Results.which;
    axh = p.Results.axh;
    
    hold(axh, 'on');
    
    if isa(channelDescriptor, 'ChannelDescriptor')
        units = channelDescriptor.unitsPrimary;
        if isempty(p.Results.label)
            label = channelDescriptor.getAxisLabelPrimary(); % includes units
        else
            label = sprintf('%s (%s)', p.Results.label, units);
        end
    elseif isstruct(channelDescriptor)
        if isempty(p.Results.label)
            name = channelDescriptor.name;
        else
            name = p.Results.label;
        end
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
    
%     au.update();
end
