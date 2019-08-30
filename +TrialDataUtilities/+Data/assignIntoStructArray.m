function S = assignIntoStructArray(S, fld, vals, idx, treatValsAsScalar)
% S = assignIntoStructArray(S, fld, vals, idx)
% for each element s in struct array S(idx), efficiently assigns
% s.(fld) = the corresponding element from vals.
%
% vals must either be a numel(S) vector or numel(S) x numel(fld)
% cell/matrix
% 
% If fld is a cellstr, assigns vals(iS, iFld) into each field in fld{:}

    if ~exist('idx', 'var') || isempty(idx)
        if isempty(S)
            % no existing struct, make it by figuring out how many vals are
            % given
%             if ~iscell(fld) || length(fld) == 1
%                 % make into vector if many fields
%                 vals = makecol(vals);
%             end
            idx = 1:size(vals, 1);
            if numel(vals) == 0
                idx = [];
            end
        else
            idx = 1:numel(S);
        end
    end
     
    if ~exist('treatValsAsScalar', 'var')
        treatValsAsScalar = false;
    end
    
    % cell wrap fld
    if ~iscell(fld)
        fld = {fld};
    end
    
    % cell wrap vals
    if ischar(vals) || isempty(vals)
        vals = {vals};
    elseif ~iscell(vals)
        vals = num2cell(vals);
    end
    
    if treatValsAsScalar && ~isscalar(vals) % this is for the case where we need to assign a scalar cell wrapped into each element, which would otherwise be treated as nElements x 1 cell array with one value for each element
        vals = {vals};
    end
    
    % scalar expand vals
    if isscalar(vals)
        vals = repmat(vals, nnz(idx), numel(fld));
    elseif isvector(vals)
        if size(vals, 1) == numel(fld)
            vals = vals';
        end
        vals = repmat(vals, 1, numel(fld));
    end 
        
    if numel(idx) == size(vals, 1)
        vals = vals(idx, :);
    elseif islogical(idx) && nnz(idx) == size(vals, 1)
        % correct size already
    else
        error('size(vals, 1) must match either numel(idx) or nnz(idx)');
    end
    
    assert(numel(fld) == size(vals, 2), ...
        'size(vals, 2) must match numel(fields)');
    
    createdS = isempty(S);    
    
    for iFld = 1:length(fld)
        [S(idx).(fld{iFld})] = deal(vals{:, iFld});
    end
    
    if createdS
        S = makecol(S);
    end

end
