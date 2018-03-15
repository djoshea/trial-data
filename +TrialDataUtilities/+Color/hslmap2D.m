function img = hslmap2D(varargin)
    % vary hue along first axis, then saturation
    
    p = inputParser();
    % to match hslmap's defaults (too much purple in default)
    p.addParameter('hue', TrialDataUtilities.Data.circspace(-0.01*360, 0.8*360, 360), @isvector);
    p.addParameter('lum', [], @isvector);
    p.addParameter('sat', [], @isvector);
    p.parse(varargin{:});
    
    hue = makecol(p.Results.hue);
    lum = makecol(p.Results.lum);
    sat = makecol(p.Results.sat);
    
    if isempty(lum) && isempty(sat)
        % default is lum and sat increase
        lum = linspace(0.1, 0.6, 100);
        sat = linspace(0.1, 1, 100);
    elseif isempty(lum)
        lum = 0.6;
    elseif isempty(sat)
        sat = 0.6;
    end
    
    if isscalar(lum)
        lum = repmat(lum, numel(sat), 1);
    end
    if isscalar(sat)
        sat = repmat(sat, numel(lum), 1);
    end
    
    [H, L] = ndgrid(hue, lum);
    [~, S] = ndgrid(hue, sat);
    
    hslImg = cat(3, H, S, L);
    
    img = TrialDataUtilities.Color.convert2D('HSL->RGB', hslImg);
end