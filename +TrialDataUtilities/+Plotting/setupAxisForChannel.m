function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('which', 'y', @(str) ismember(str, {'x', 'y', 'z'}));
    p.addParamValue('useAutoAxis', true, @islogical);
    p.addParamValue('scaleBar', false, @islogical);
    p.parse(varargin{:});
    
    which = p.Results.which;
    axh = p.Results.axh;
    
    if p.Results.scaleBar
        label = channelDescriptor.name;
    else
        label = channelDescriptor.getAxisLabelPrimary();
    end
    
    if p.Results.useAutoAxis
        au = AutoAxis(axh);
    end
    switch which
        case 'x'
            xlabel(axh, label);
            if p.Results.useAutoAxis
                if p.Results.scaleBar
                    au.addAutoScaleBarX();
                else
                	au.addAutoAxisX();
                end
            end
        case 'y'
            ylabel(axh, label);
            if p.Results.useAutoAxis
                if p.Results.scaleBar
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
        %au.installCallbacks();
    end
end
