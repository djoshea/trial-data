function graphicsSortInPlaceAs(hvec)
    % given a list of Graphics handles, resort them visually in the order specified in hvec
    % (first object is on top), but stores the entire set at the position in axis.Children 
    % where the presently top-most object in hvec is stored. So this sets the relative order of the
    % items in hvec, but won't bring them on top of any other objects in the scene that were on top of all
    % hvec previously.
    
    if isempty(hvec)
        return;
    end
    
    assert(isa(hvec, 'matlab.graphics.Graphics'));
    
    mask = arrayfun(@(h) ~isa(h, 'matlab.graphics.GraphicsPlaceholder'), hvec);
    if ~any(mask)
        return;
        warning('No non-placeholders found');
    end
    hvec = hvec(mask);
    ax = TrialDataUtilities.Plotting.getParentAxis(hvec(1));
    [tf, idx] = ismember(hvec, ax.Children);
    
    assert(all(tf), 'Graphics objects must be in same axis');

    idxInsert = min(idx);
    childrenOther = ax.Children(setdiff(1:numel(ax.Children), idx));
    ax.Children = [childrenOther(1:idxInsert-1); hvec; childrenOther(idxInsert:end)];
end