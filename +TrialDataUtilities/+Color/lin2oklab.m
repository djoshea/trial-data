function ok = lin2oklab(lin)
  % LIN2OKLAB Convert colors stored in linear sRGB color space to OKLab color
  % space
  %
  % ok = lin2oklab(lin)
  %
  % Inputs:
  %   lin  #h by #w by 3 image of colors or #lin by 3 list of colors
  % Outputs:
  %   ok  converted ok colors (same size as lin)
  %
  % See also: rgb2lin, lin2rgb, oklab2rgb, rgb2oklab, oklab2lin
  %
  % @misc{gptoolbox,
  %   title = {{gptoolbox}: Geometry Processing Toolbox},
  %   author = {Alec Jacobson and others},
  %   note = {http://github.com/alecjacobson/gptoolbox},
  %   year = {2021},
  % }
  % upodated with new matrices from https://bottosson.github.io/posts/oklab/

  was_permuted = false;
  if size(lin,3) == 1 && size(lin,2) == 3
    lin = permute(lin,[1 3 2]);
    was_permuted = true;
  end

  lms = cat(3, ...
    0.4122214708*lin(:,:,1) + 0.5363325363*lin(:,:,2) + 0.0514459929*lin(:,:,3), ...
    0.2119034982*lin(:,:,1) + 0.6806995451*lin(:,:,2) + 0.1073969566*lin(:,:,3), ...
    0.0883024619*lin(:,:,1) + 0.2817188376*lin(:,:,2) + 0.6299787005*lin(:,:,3));
  lms = lms.^(1/3);
  ok = cat(3, ...
    0.2104542553*lms(:,:,1) + 0.7936177850*lms(:,:,2) - 0.0040720468*lms(:,:,3), ...
    1.9779984951*lms(:,:,1) - 2.4285922050*lms(:,:,2) + 0.4505937099*lms(:,:,3), ...
    0.0259040371*lms(:,:,1) + 0.7827717662*lms(:,:,2) - 0.8086757660*lms(:,:,3));

  if was_permuted
    ok = permute(ok,[1 3 2]);
  end
end