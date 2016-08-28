function cvec = desaturate(c, gamma)
    if nargin < 2
        gamma = 0.2;
    end
    c = TrialDataUtilities.Color.toRGB(c);
    cvec = 1 - ((1-c) .* (1-gamma));
end
