function img = phaseAmplitude2DMap(varargin)

p = inputParser();
p.addOptional('nPhase', 360, @isscalar);
p.addOptional('nAmplitude', 60, @isscalar);
p.addParameter('varyChroma', true, @islogical); % if false, no loss of saturation
p.addParameter('varyLuminance', true, @islogical); % if false, everything fades to gray instead of black
p.addParameter('chromaGains', [0 1], @isvector); % either limits for linspace or nAmplitude x 1 vector
p.addParameter('luminanceGains', [0 1], @isvector); % either limits for linspace or nAmplitude x 1 vector
p.addParameter('keepInLCH', false, @islogical); % keep in LCH color space
p.parse(varargin{:});

nPhase = p.Results.nPhase;
nAmplitude = p.Results.nAmplitude;

% uses Chad Green's cmocean's phase map, originally developed by created by Kristen Thyng
phase = TrialDataUtilities.Color.cmocean('phase', nPhase);

phaseLCH = TrialDataUtilities.Color.convert('RGB->LCH', phase);

phaseLCH_start = phaseLCH;

phaseLCH = repmat(phaseLCH, [nAmplitude 1]);

chromaGains = p.Results.chromaGains;
if numel(chromaGains) == 2
    chromaGains = linspace(chromaGains(1), chromaGains(2), nAmplitude);
else
    assert(numel(chromaGains) == nAmplitude);
end

luminanceGains = p.Results.luminanceGains;
if numel(luminanceGains) == 2
    luminanceGains = linspace(luminanceGains(1), luminanceGains(2), nAmplitude);
else
    assert(numel(luminanceGains) == nAmplitude);
end

gain = linspace(0, 1, nAmplitude);
for iA = 1:nAmplitude
    offset = (iA-1)*nPhase;
    if p.Results.varyLuminance
        phaseLCH((1:nPhase) + offset, 1) = phaseLCH_start(:, 1) * luminanceGains(iA);
    end
    if p.Results.varyChroma
        phaseLCH((1:nPhase) + offset, 2) = phaseLCH_start(:, 2) * chromaGains(iA);
    end
end

% we'll use the first dim to vary luminance
if ~p.Results.keepInLCH
    phase = TrialDataUtilities.Color.convert('LCH->RGB', phaseLCH);
else
    phase = phaseLCH;
end

img = reshape(phase, [nPhase, nAmplitude, 3]);

end
