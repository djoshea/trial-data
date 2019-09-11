function intervalCell = removeOverlappingIntervals(varargin)
    % intervalCell = removeOverlappingIntervals(intervalCell, [intervalCell2])
    % 
    % intervalCell2 is optional, will be combined with intervalCell if
    % specified
    %
    % both intervalCell are cells of matching size with nInterval x 2 matrices inside
    % 
    % simplifies the blanking regions to be non-overlapping
   
    if numel(varargin) > 1
        % combine the intervals in all windows
        intervalCell = varargin{1};
        
        for iArg = 2:numel(varargin)
            intervalCell2 = varargin{iArg};
            for iT = 1:numel(intervalCell)
                if isempty(intervalCell{iT})
                    intervalCell{iT} = intervalCell2{iT};
                elseif ~isempty(intervalCell2{iT})
                    intervalCell{iT} = cat(1, intervalCell{iT}, intervalCell2{iT});
                end
            end
        end
    else
        intervalCell = varargin{1};
    end
    
    for iT = 1:numel(intervalCell)
        mat = intervalCell{iT};
       
        maskRemove = falsevec(size(mat, 1));
        newMat = mat;
        for i = 1:size(mat, 1)
            for j = i+1:size(mat, 1)
                if mat(i, 2) >= mat(j, 1) || ... % i overlaps left of j
                   mat(i, 2) >= mat(j, 1) % i overlaps right of j
                    % i and j overlap
                    newMat(i, 1) = min(mat(i,1), mat(j,1));
                    newMat(i, 2) = max(mat(i,2), mat(j,2));
                    maskRemove(j) = true;
                end
            end
        end
        
        intervalCell{iT} = newMat(~maskRemove, :);
    end
end