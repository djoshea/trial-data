function lab = oklch2oklab(lch)
  % OKLCH2LAB Convert colors stored in OKLCH to OKLab color space 
  %
  % oklab = oklch2lab(oklch)
  %
  % Inputs:
  %   oklch  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   lab  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2rgb, rgb2oklab, lin2oklab

  was_permuted = false;
  if size(lch,3) == 1 && size(lch,2) == 3
    lch = permute(lch,[1 3 2]);
    was_permuted = true;
  end

  % L = L, a = C cosd(h), b = C sind(h)
  lab = cat(3, ...
    lch(:, :, 1), ...
    lch(:, :, 2) .* cosd(lch(:, :, 3)), ...
    lch(:, :, 2) .* sind(lch(:, :, 3)));

  if was_permuted
    lab = permute(lab,[1 3 2]);
  end

end
