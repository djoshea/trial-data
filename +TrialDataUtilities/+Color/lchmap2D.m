function img = lchmap2D(varargin)
    % vary hue along first axis, then saturation
    
    p = inputParser();
    p.addParameter('hue', TrialDataUtilities.Data.circspace(0, 360, 360), @isscalar);
    p.addParameter('lum', [], @isvector);
    p.addParameter('chroma', [], @isvector);
    p.parse(varargin{:});
    
    hue = makecol(p.Results.hue);
    lum = makecol(p.Results.lum);
    chroma = makecol(p.Results.chroma);
    
    if isempty(lum) && isempty(chroma)
        % default is lum and sat increase
        lum = linspace(10, 80, 100);
        chroma = linspace(10, 80, 100);
    elseif isempty(lum)
        lum = 80;
    elseif isempty(chroma)
        chroma = 80;
    end
    
    if isscalar(lum)
        lum = repmat(lum, numel(chroma), 1);
    end
    if isscalar(chroma)
        chroma = repmat(chroma, numel(lum), 1);
    end
    
    [H, L] = ndgrid(hue, lum);
    [~, C] = ndgrid(hue, chroma);
    
    lchImg = cat(3, L, C, H);
    
    img = TrialDataUtilities.Color.convert2D('LCH->RGB', lchImg);
end