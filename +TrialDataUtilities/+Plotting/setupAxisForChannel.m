function setupAxisForChannel(channelDescriptor, varargin)
    p = inputParser();
    p.addParamValue('axh', gca, @ishandle);
    p.addParamValue('xy', 'y', @(str) ismember(str, {'x', 'y'}));
    p.parse(varargin{:});
    
    xy = p.Results.xy;
    axh = p.Results.axh;
    
    label = channelDescriptor.getAxisLabelPrimary();
    
    au = AutoAxis(axh);
    if strcmp(xy, 'x')
        xlabel(axh, label);
        au.addAutoAxisX();
    else
        ylabel(axh, label);
        au.addAutoAxisY();
    end
    
    au.update();
    au.installCallbacks();
end