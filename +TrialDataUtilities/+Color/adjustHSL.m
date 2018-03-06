function map = adjustHSL(map, varargin)

    p = inputParser();
    p.addParameter('addL', 0, @isscalar);
    p.addParameter('multL', 1, @isscalar);
    p.addParameter('addS', 0, @isscalar);
    p.addParameter('multS', 1, @isscalar);
    p.addParameter('addH', 0, @isscalar); % degrees
    p.parse(varargin{:});

    mapHSL = TrialDataUtilities.Color.convert('RGB->HSL', map);

    mapHSL(:, 1) = mod(mapHSL(:, 1) + p.Results.addH, 360);
    mapHSL(:, 2) = mapHSL(:, 2) * p.Results.multS + p.Results.addS;
    mapHSL(:, 3) = mapHSL(:, 3) * p.Results.multL + p.Results.addL;

    map = TrialDataUtilities.Color.convert('HSL->RGB', mapHSL);

end