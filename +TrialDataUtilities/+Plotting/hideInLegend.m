function hideInLegend(h)
    % prevent object h from appearing in legend by default
    for i = 1:numel(h)
        if ~TrialDataUtilities.Plotting.isGraphicsHandle(h(i)), continue; end
        if isprop(h(i), "Annotation")
            ann = get(h(i), 'Annotation');
            leg = get(ann, 'LegendInformation');
            set(leg, 'IconDisplayStyle', 'off');
        end
    end
end