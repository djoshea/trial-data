function showFirstInLegend(h, name)
% showInLegend(h, names)
% show first object in h in default legend with specified name
% use legend(axh, 'show') to activate default legend

    h = h(~isnan(h));
    if ~isempty(h)
        TrialDataUtilities.Plotting.showInLegend(h(1), name);
        if numel(h) > 1
            TrialDataUtilities.Plotting.hideInLegend(h(2:end));
        end
    end
end