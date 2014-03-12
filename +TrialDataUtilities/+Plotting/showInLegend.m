function showInLegend(h, names)
% showInLegend(h, names)
% show object h in legend with specified by default
% use legend(axh, 'show') to activate default legend
% names is either char (for scalar h), or cellstr

    if ~iscell(names)
        names = {names};
    end

    for i = 1:numel(h)
        ann = get(h(i), 'Annotation');
        leg = get(ann, 'LegendInformation');
        set(leg, 'IconDisplayStyle', 'on');
        
        set(h(i), 'DisplayName', names{i});
    end

end