function lch = oklab2oklch(lab)
  % OKLAB2LCH Convert colors stored in OKLab to OKLCH color space 
  %
  % oklch = oklab2oklch(oklab)
  %
  % Inputs:
  %   oklab  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   lch  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2rgb, rgb2oklab, lin2oklab

  was_permuted = false;
  if size(lab,3) == 1 && size(lab,2) == 3
    lab = permute(lab,[1 3 2]);
    was_permuted = true;
  end

  % L = L, a = C cosd(h), b = C sind(h)
  lch = cat(3, ...
    lab(:, :, 1), ...
    sqrt(lab(:, :, 2).^2 + lab(:, :, 3).^2), ...
    atan2d(lab(:, :, 3), lab(:, :, 2)));

  if was_permuted
    lch = permute(lch, [1 3 2]);
  end

end
