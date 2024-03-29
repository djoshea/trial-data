function cmap = BlueWhiteRed(n)

if nargin < 1
    n = 301;
end

red = [185 39 50] / 255;
% white = [247 246 246] / 255;	
white = [1 1 1];
blue = [41 113 177] / 255;
cmap_internal = TrialDataUtilities.Color.blendcolors([blue; white; red], n);

if nargout == 0
    clim = caxis;
    cmax = max(abs(clim));
    colormap(cmap_internal);
    caxis([-cmax cmax]);
else
    cmap = cmap_internal;
end

end
    