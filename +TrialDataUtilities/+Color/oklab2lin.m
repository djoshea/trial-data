function lin = oklab2lin(ok)
  % OKLAB2LIN Convert colors stored in OKLab color space to linear sRGB color space
  %
  % lin = oklab2lin(ok)
  %
  % Inputs:
  %   ok  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   lin  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2rgb, rgb2oklab, lin2oklab
  %
  % @misc{gptoolbox,
  %   title = {{gptoolbox}: Geometry Processing Toolbox},
  %   author = {Alec Jacobson and others},
  %   note = {http://github.com/alecjacobson/gptoolbox},
  %   year = {2021},
  % }
  % updated with new 2021-01-25 matrices from https://bottosson.github.io/posts/oklab/

  was_permuted = false;
  if size(ok,3) == 1 && size(ok,2) == 3
    ok = permute(ok,[1 3 2]);
    was_permuted = true;
  end

  lms = cat(3, ...
    ok(:,:,1)+0.3963377774*ok(:,:,2)+0.2158037573*ok(:,:,3), ...
    ok(:,:,1)-0.1055613458*ok(:,:,2)-0.0638541728*ok(:,:,3), ...
    ok(:,:,1)-0.0894841775*ok(:,:,2)-1.2914855480*ok(:,:,3));
  lms = lms.^3;
  lin = cat(3, ...
    +4.0767416621*lms(:,:,1) - 3.3077115913*lms(:,:,2) + 0.2309699292*lms(:,:,3), ...
    -1.2684380046*lms(:,:,1) + 2.6097574011*lms(:,:,2) - 0.3413193965*lms(:,:,3), ...
    -0.0041960863*lms(:,:,1) - 0.7034186147*lms(:,:,2) + 1.7076147010*lms(:,:,3));

  if was_permuted
    lin = permute(lin,[1 3 2]);
  end
end