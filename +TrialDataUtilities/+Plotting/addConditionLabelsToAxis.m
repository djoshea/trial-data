function addConditionLabelsToAxis(tdca, axh, varargin)

    au = AutoAxis(axh);
    labels = tdca.conditionNamesShort;
    colors = {tdca.conditionAppearances.Color};
    au.addColoredLabels(labels, colors, varargin{:});
    au.update();

end