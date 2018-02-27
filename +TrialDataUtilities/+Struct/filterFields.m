function [sfilt, flds, fldsMask] = filterFields(s, fn, nested)
    % fn should have signature tf = fn(value, fieldName)

    if nargin < 3
        nested = false;
    end
    
    flds = fieldnames(s);
    sfilt = struct();
    fldsMask = false(numel(flds), 1);
    for iF = 1:numel(flds)
        if nested && isstruct(s.(flds{iF}))
            sfilt.(flds{iF}) = TrialDataUtilities.Struct.filterFields(s.(flds{iF}), fn, nested);
            
        elseif fn(s.(flds{iF}), flds{iF})
            sfilt.(flds{iF}) = s.(flds{iF});
            fldsMask(iF) = true;
        end
    end
    flds = flds(fldsMask);
end
