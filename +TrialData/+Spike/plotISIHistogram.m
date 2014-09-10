function plotISIHistogram(td, unit, varargin)
    p = inputParser();
    p.addParamValue('binWidth', 0.5, @isscalar);
    p.addParamValue('maxISI', 30, @isscalar);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    [freq, binEdges] = TrialData.Spike.computeISIHistogram(td, unit, p.Results);
    
    %freq = freq ./ sum(freq);

    TrialDataUtilities.Plotting.plotHistogram(binEdges, freq, p.Unmatched);
    xlabel(sprintf('ISI (%s)', td.timeUnitName));
    %ylabel('Proportion');
    ylabel('Spike count');
    %title(sprintf('ISI Histogram for %s', unit));
    
    au = AutoAxis(gca);
    au.replace();
end