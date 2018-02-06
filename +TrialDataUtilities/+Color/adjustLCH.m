function map = adjustLCH(map, varargin)

    p = inputParser();
    p.addParameter('addL', 0, @isscalar);
    p.addParameter('multL', 1, @isscalar);
    p.addParameter('addC', 0, @isscalar);
    p.addParameter('multC', 1, @isscalar);
    p.addParameter('addH', 0, @isscalar); % degrees
    p.parse(varargin{:});

    mapLCH = TrialDataUtilities.Color.convert('RGB->LCH', map);

    mapLCH(:, 1) = mapLCH(:, 1) * p.Results.multL + p.Results.addL;
    mapLCH(:, 2) = mapLCH(:, 2) * p.Results.multC + p.Results.addC;
    mapLCH(:, 3) = mod(mapLCH(:, 3) + p.Results.addH, 360);

    map = TrialDataUtilities.Color.convert('LCH->RGB', mapLCH);

end