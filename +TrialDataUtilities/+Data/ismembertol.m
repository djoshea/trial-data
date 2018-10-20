function [lia, locb] = ismembertol(a, b, tol)
% same as uniquetol with a fallback for versions prior to R2015a

    if iscategorical(a)
        if ~iscategorical(b)
            b = categorical(b);
        end
        [lia, locb] = ismember(a, b);
        return
    end

    if isfloat(a) && isfloat(b)
        if nargin < 3
            tol = 1e-6;
        end

        if true || verLessThan('matlab', '8.5')
            equiv = abs(bsxfun(@minus, makecol(a), makerow(b))) < tol;
            lia = any(equiv, 2);
            locb = nan(size(a));
            [~, locb(lia)] = max(equiv(lia, :), [], 2);
        else
            [lia, locb] = builtin('ismembertol', makecol(a(:)), makecol(b(:)), tol);
        end
    else
        % need one arg to be double if classes don't match
        if (islogical(b) || isnumeric(b)) && ~isa(b, 'double')
            b = double(b);
        end
        [lia, locb] = ismember(a, b); % need one arg to be double
    end
end