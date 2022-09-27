function rgb = oklch2rgb(ok)
  % OKLCH2RGB Convert colors stored in OKLab color space to RGB color space
  %
  % rgb = oklch2rgb(ok)
  %
  % Inputs:
  %   ok  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   rgb  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2lin , rgb2oklab, lin2oklab
  %
  % @misc{gptoolbox,
  %   title = {{gptoolbox}: Geometry Processing Toolbox},
  %   author = {Alec Jacobson and others},
  %   note = {http://github.com/alecjacobson/gptoolbox},
  %   year = {2021},
  % }

  lab = TrialDataUtilities.Color.oklch2oklab(ok);
  rgb = TrialDataUtilities.Color.oklab2rgb(lab);
end