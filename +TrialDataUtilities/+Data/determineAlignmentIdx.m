function [idx_take, idx_insert, coords_aligned, coords_valid] = determineAlignmentIdx(coord_set, args)
% this function determines how to take slices from each element of coord_set so as to align the coordinates
% in a concatenated output. This function assumes even, shared spacing of coord_set vectors
    arguments
        coord_set (:, 1) cell % cell array of indices
        args.valid_masks (:, 1) cell = {};
        args.output_range (1, 2) = [-Inf Inf];
        args.overlap_only (1, 1) logical = false;
    end

    nC = numel(coord_set);
    dx = coord_set{1}(2) - coord_set{1}(1);

    if ~isempty(args.valid_masks)
        for iC = 1:nC
            coord_set{iC}(~args.valid_masks{iC}) = NaN;
        end
    end

    if ~args.overlap_only
        lo = max(args.output_range(1), min(cellfun(@(x) min(x, [], 'omitnan'), coord_set)));
        hi = min(args.output_range(2), max(cellfun(@(x) max(x, [], 'omitnan'), coord_set)));
    else
        lo = max(args.output_range(1), max(cellfun(@(x) min(x, [], 'omitnan'), coord_set)));
        hi = min(args.output_range(2), min(cellfun(@(x) max(x, [], 'omitnan'), coord_set)));
    end
    coords_aligned = (lo:dx:hi)';
    if isempty(coords_aligned)
        error('No overlapping region found in coord_set');
    end
    coords_valid = false(numel(coords_aligned), nC);

    [idx_take, idx_insert] = cellvec(nC);
    for iC = 1:nC
        [tf, where] = ismembertol(coord_set{iC}, coords_aligned, dx/10, DataScale=1);
        idx_take{iC} = find(tf);
        idx_insert{iC} = where(tf);
        coords_valid(where(tf), iC) = true;
    end
end

%% Testing
%coord_set = { 1:5, 3:10, 4:8};
%[idx_take, idx_insert, coords_aligned, coords_valid] = determineAlignmentIdx(coord_set, args)



