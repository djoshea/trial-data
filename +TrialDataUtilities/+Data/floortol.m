function a = floortol(a, tol)
    if nargin < 2
        tol = 1e-6;
    end
    a = floor(a + tol);
end