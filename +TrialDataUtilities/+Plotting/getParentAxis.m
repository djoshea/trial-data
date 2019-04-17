function ax = getParentAxis(h)
    assert(isa(h, 'matlab.graphics.Graphics'));
    ax = h(1);
    while ~isempty(ax) && ~isa(ax, 'matlab.graphics.axis.Axes')
      ax = ax.Parent;
    end
    if isempty(ax) || ~isa(ax, 'matlab.graphics.axis.Axes')
        error('Could not determine parent axis');
    end
end