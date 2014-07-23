function plotISIHistogram(td, unit, varargin)
    [freq, binEdges] = TrialData.Spike.computeISIHistogram(td, unit, varargin{:});
    
    freq = freq ./ sum(freq);

    stairs(binEdges, freq);
    ylabel(sprintf('ISI (%s)', td.timeUnitName));
    xlabel('Proportion');
    title(sprintf('ISI Histogram for %s', unit));
    
    au = AutoAxis(gca);
    au.replace();
end