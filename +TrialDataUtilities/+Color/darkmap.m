function cmap = darkmap(color, n, varargin)
    p = inputParser();
    p.addParameter('dark', [0.2 0.2 0.2], @isvector);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    color = TrialDataUtilities.Color.toRGB(color);
    dark = TrialDataUtilities.Color.toRGB(p.Results.dark);
    cmap = TrialDataUtilities.Color.blendcolors([dark; color], n, p.Unmatched);
end
