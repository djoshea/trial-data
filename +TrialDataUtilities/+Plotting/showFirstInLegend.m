function showFirstInLegend(h, name)
% showInLegend(h, names)
% show first object in h in default legend with specified name
% use legend(axh, 'show') to activate default legend

    if isempty(name)
        TrialDataUtilities.Plotting.hideInLegend(h);
    else
        h = h(TrialDataUtilities.Plotting.isGraphicsHandle(h));
        if ~isempty(h)
            TrialDataUtilities.Plotting.showInLegend(h(1), name);
            if numel(h) > 1
                TrialDataUtilities.Plotting.hideInLegend(h(2:end));
            end
        end
    end
end