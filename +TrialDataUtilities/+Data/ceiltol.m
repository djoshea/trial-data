function a = ceiltol(a, tol)
    if nargin < 2
        tol = 1e-6;
    end
    a = ceil(a - tol);
end