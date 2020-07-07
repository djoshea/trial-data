function cmap = brighten(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    
    cmap = TrialDataUtilities.Color.darken(cmap, -amount);
end