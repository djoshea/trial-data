function cmap = flare_reversed(nColors)
% from Seaborn, (c) Michael Waskom.  (c.f. seaborn/seaborn/cm.py)

if nargin < 1
    nColors = [];
end
cmap = flipud(TrialDataUtilities.Colormaps.flare(nColors));

if nargout == 0
    colormap(cmap);
end

end