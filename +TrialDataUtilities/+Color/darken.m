function cmap = darken(cmap, factor)

    cmap = TrialDataUtilities.Color.adjustHSL(cmap, 'multL', factor);

end