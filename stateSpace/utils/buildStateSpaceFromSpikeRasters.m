function [pset mask] = buildStateSpaceFromSpikeRasters(srCell, varargin)

p = inputParser();
p.addParamValue('minTrialsByNonEmptyCondition', 6, @isscalar);
p.addParamValue('minMeanHz', 5, @isscalar);
p.KeepUnmatched = true;
p.parse(varargin{:});

% filter for SpikeRasters

isValidFn = @(sr) isa(sr, 'SpikeRaster') && ...
    sr.minTrialsByNonEmptyCondition >= p.Results.minTrialsByNonEmptyCondition && ...
    sr.getMeanHz() >= p.Results.minMeanHz;

% all rasters in a given row must be valid to use
mask = all(cellfun(isValidFn, srCell), 2);
srCell = srCell(mask, :);

pset = PopulationTrajectorySet(p.Unmatched);
pset.addUnitFromSpikeRaster(srCell);

end