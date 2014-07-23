function [freq, binEdges] = computeISIHistogram(td, unit, varargin)
    % using data only within the the current alignment window,
    % compute the iti histogram for the spike unit
    assert(isa(td, 'TrialData'));
    
    p = inputParser();
    p.addParamValue('binWidth', 0.5, @isscalar);
    p.addParamValue('maxISI', 40, @isscalar);
    p.parse(varargin{:});

    timesCell = td.getSpikeTimes(unit);
    isiVec = cellfun(@(times) diff(times), timesCell, 'UniformOutput', false);

    isiCat = cat(1, isiVec{:});

    binEdges = 0:p.Results.binWidth:p.Results.maxISI;
    freq = histc(isiCat, binEdges);
end