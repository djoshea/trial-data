function newClass = cellDetermineCommonClass(data, origClass)
    if nargin > 1
        newClass = origClass;
    else
        newClass = '';
    end
    nonEmpty = ~cellfun(@isempty, data);
    if ~any(nonEmpty)
        newClass = origClass;
    else
        % support recursive cell of cell of .. of class(data)
        classes = string(cellfun(@getCellElementClass, data(nonEmpty), 'UniformOutput', false));
        
        newClass = TrialDataUtilities.Data.determineCommonClass(classes{:}, newClass);
    end
end

function cl = getCellElementClass(dcell)
    if ~iscell(dcell)
        cl = class(dcell);
    else
        cl = '';
        for i = 1:numel(dcell)
            if ~isempty(dcell{i})
                cl = class(dcell{i});
                return;
            end
        end
    end
end