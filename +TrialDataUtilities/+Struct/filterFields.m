function sfilt = filterFields(s, fn)
    % fn should take tf = fn(value, fieldName)

    flds = fieldnames(s);
    sfilt = struct();
    for iF = 1:numel(flds)
        if fn(s.(flds{iF}), flds{iF})
            sfilt.(flds{iF}) = s.(flds{iF});
        end
    end
end
