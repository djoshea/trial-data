function cmap = phaseMapEvalAt(phaseValues, varargin)

p = inputParser();
p.addParameter('degrees', false, @islogical);
p.addParameter('phaseShift', 0, @isscalar);
p.parse(varargin{:});

phaseShift = p.Results.phaseShift;
if p.Results.degrees
    phaseShift = deg2rad(phaseShift);
    phaseValues = deg2rad(phaseValues);
end

% uses Chad Green's cmocean's phase map, originally developed by created by Kristen Thyng
N = 360;
map = TrialDataUtilities.Color.cmocean('phase', N);
cmapLims = [0 2*pi];
evalAt = mod(phaseValues + phaseShift, 2*pi);

cmap = TrialDataUtilities.Color.evalColorMapAt(map, evalAt, cmapLims);

