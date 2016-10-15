function [xc, edges] = computeCrossCorrelation(td, unitA, unitB, varargin)

p = inputParser();
p.addParameter('edges', 0:0.1:20, @isvector);
p.addParameter('showPlot', false, @islogical);
p.parse(varargin{:});

timesA = td.getSpikeTimes(unitA, 'combine', true);
timesB = td.getSpikeTimes(unitB, 'combine', true);

edges = p.Results.edges;

xc = zeros(numel(edges)-1, 1);
for t = 1:td.nTrials
    tA = timesA{t};
    tB = timesB{t};

    dist = pdist2(tA, tB, 'cityblock');
    
    xc = xc + histcounts(dist(:), edges)';
end

if p.Results.showPlot
    xcn = xc ./ sum(xc);
    stairs(edges(1:end-1), xcn, 'LineWidth', 1);
    xlabel('{\Delta}Time (ms)');
    ylabel('Density');
    AutoAxis.replace();
end
