function [dataCell, timeCell] = resampleDataCellAtSpecificTimes(dataCell, timeCell, sampleAtCell, varargin)
% this is a shorthand for calling resampleTensorInTime in a loop with some
% double checks to handle truncation when padding is used

    p = inputParser();
   p.addParameter('interpolateMethod', 'linear', @ischar); % pchip, linear, spline
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('extrapolate', false); 
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
    
    % these are used as a secondary guard to truncate data within tMin :
    % tMax, when the input data includes padded edges to facilitate
    % resampling
    p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
    p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
    p.parse(varargin{:});
    
    tMinExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMinExcludingPadding, size(dataCell));
    tMaxExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMaxExcludingPadding, size(dataCell));
   
    timeCell = TensorUtils.singletonExpandToSize(timeCell, size(dataCell));
    
    if ~isempty(p.Results.origDelta)
        tol = p.Results.orig / 1000;
    else
        tol = 1e-8;
    end
    
    if ~iscell(sampleAtCell)
        sampleAtCell = repmat({sampleAtCell}, numel(dataCell), 1);
    end
    
    %prog = ProgressBar(numel(dataCell), 'Resampling data');
    for iD = 1:numel(dataCell)
        %prog.update(iD);
        if isempty(dataCell{iD}) || size(dataCell{iD}, 1) == 1, continue; end % test for empty or single-sample
        [d, t] = TrialDataUtilities.Data.resampleTensorAtSpecificTimes(dataCell{iD}, 1, timeCell{iD}, sampleAtCell{iD}, ...
            'origDelt', p.Results.origDelta, ...
            'binAlignmentMode', p.Results.binAlignmentMode, ...
            'interpolateMethod', p.Results.interpolateMethod, ...
            'extrapolate', p.Results.extrapolate);
        
        mask = t >= tMinExcludingPadding(iD) - tol & t <= tMaxExcludingPadding(iD) + tol;
        dataCell{iD} = d(mask, :, :, :);
        timeCell{iD} = t(mask);
    end
    %prog.finish();
end

