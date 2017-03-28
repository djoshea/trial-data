function [freq, binEdges] = computeISIHistogram(td, unit, varargin)
    % using data only within the the current alignment window,
    % compute the iti histogram for the spike unit
    assert(isa(td, 'TrialData'));
    
    p = inputParser();
    p.addParamValue('binWidth', 0.5, @isscalar);
    p.addParamValue('maxISI', 30, @isscalar);
    p.parse(varargin{:});

    timesCell = td.getSpikeTimes(unit);
    mask = cellfun(@(times) numel(times) > 1, timesCell);
    isiVec = cellfun(@(times) makecol(diff(times)), timesCell(mask), 'UniformOutput', false);

    isiCat = cat(1, isiVec{:});
    
    if isempty(isiCat)
        error('No spike pairs found for unit %s', unit);
    end

    binEdges = 0:p.Results.binWidth:p.Results.maxISI;
    freq = histc(isiCat, binEdges);
end