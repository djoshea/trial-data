function cmap = saturate(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    cmap = TrialDataUtilities.Color.desaturate(cmap, -amount);
end
