function [lia, locb] = ismembertol(a, b, tol)
% same as uniquetol with a fallback for versions prior to R2015a

    if isnumeric(a) && isnumeric(b)
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
        [lia, locb] = ismember(a, b);
    end
end