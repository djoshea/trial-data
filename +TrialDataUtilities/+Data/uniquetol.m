function [c, ia, ic] = uniquetol(a, tol)
% same as uniquetol with a fallback for versions prior to R2015a
    if nargin < 3
        tol = 1e-6 * nanmax(abs(a(:)));
    end

    if isempty(a)
        c = [];
        ia = [];
        ic = [];
        return;
    end
    if verLessThan('matlab', '8.5')
        c = builtin('_mergesimpts',makecol(a(:)),tol);
        if nargout > 1
            error('Multiple outputs not supported for earlier Matlab versions');
        end
    else
        [c, ia, ic] = builtin('uniquetol', makecol(a(:)), tol);
    end
end