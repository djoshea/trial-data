function cmap = phaseAmplitude2DMapEvalAt(phaseValues, amplitudeValues, varargin)
% evaluates a 2d hue and chroma/luminance map with phase from 0:2*pi and
% amplitude from 0 to 1

p = inputParser();
p.addParameter('degrees', false, @islogical);
p.addParameter('phaseShift', 0, @isscalar);
p.addParameter('amplitudeLimits', [0 1], @isvector);
p.KeepUnmatched = true;
p.parse(varargin{:});

phaseShift = p.Results.phaseShift;
if p.Results.degrees
    phaseShift = deg2rad(phaseShift);
    phaseValues = deg2rad(phaseValues);
end

% uses Chad Green's cmocean's phase map, originally developed by created by Kristen Thyng
mapImg = TrialDataUtilities.Color.phaseAmplitude2DMap(360, 100, p.Unmatched, 'keepInLCH', false);
phaseLims = [0 2*pi];
evalPhaseAt = mod(phaseValues + phaseShift, 2*pi);

cmap = TrialDataUtilities.Color.evalColorMap2DAt(mapImg, evalPhaseAt, amplitudeValues, phaseLims, p.Results.amplitudeLimits);

