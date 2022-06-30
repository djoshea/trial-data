function h = plot_cols_with_colormap(x, y, cmap, varargin)
   nTraces = size(y, 2);
   assert(size(cmap, 1) == nTraces, 'cmap must have size(y, 2) rows');

   h = plot(x, y, varargin{:});

   for iH = 1:nTraces
       h(iH).Color = cmap(iH, :);
   end
end