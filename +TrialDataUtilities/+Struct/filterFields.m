function [sfilt, flds, fldsMask] = filterFields(s, fn)
    % fn should have signature tf = fn(value, fieldName)

    flds = fieldnames(s);
    sfilt = struct();
    fldsMask = false(numel(flds), 1);
    for iF = 1:numel(flds)
        if fn(s.(flds{iF}), flds{iF})
            sfilt.(flds{iF}) = s.(flds{iF});
            fldsMask(iF) = true;
        end
    end
    flds = flds(fldsMask);
end
