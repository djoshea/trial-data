function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('which', 'y', @(str) ismember(str, {'x', 'y', 'z'}));
    p.addParamValue('useAutoAxis', true, @islogical);
    p.parse(varargin{:});
    
    which = p.Results.which;
    axh = p.Results.axh;
    
    label = channelDescriptor.getAxisLabelPrimary();
    
    if p.Results.useAutoAxis
        au = AutoAxis(axh);
    end
    switch which
        case 'x'
            xlabel(axh, label);
            if p.Results.useAutoAxis
                au.addAutoAxisX();
            end
        case 'y'
            ylabel(axh, label);
            if p.Results.useAutoAxis
                au.addAutoAxisY();
            end
        case 'z'
            zlabel(axh, label);
    end

    if p.Results.useAutoAxis
        %au.update();
        %au.installCallbacks();
    end
end
