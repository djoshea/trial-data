function cmap = darken(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    
    Kn = 18;
    lab = TrialDataUtilities.Color.convert('RGB->Lab', cmap);

    lab(:, 1) = lab(:, 1) - Kn * amount;

    cmap = TrialDataUtilities.Color.convert('Lab->RGB', lab);

    % old version
    %     cmap = TrialDataUtilities.Color.adjustHSL(cmap, 'multL', factor);
end