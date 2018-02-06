function map = phaseMap(N, varargin)
% uses Chad Green's cmocean's phase map, originally developed by created by Kristen Thyng

p = inputParser();
p.addParameter('phaseShift', 0, @isscalar);
p.parse(varargin{:});

if nargin < 1
    N = 256;
end
shift = floor(mod(p.Results.phaseShift, 2*pi)/(2*pi) * N);

map = TrialDataUtilities.Color.cmocean('phase', N);

map = circshift(map, shift, 1);

