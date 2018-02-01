function [data, timeNew] = resampleTensorAtSpecificTimes(data, timeDim, time, sampleAtTimes, varargin)
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar); % pchip, linear, spline
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('extrapolate', false); 
%     p.addParameter('binAlignmentModeOutput', [], @(x) isempty(x) || isa(x, 'BinAlignmentMode'));
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
%     p.addParameter('newDelta', 0, @isscalar); % if the output data should be resampled with the same bin alignment mode 
    p.parse(varargin{:});
    
    assert(isvector(time));
    assert(isvector(sampleAtTimes));
    nTime = numel(time);
    assert(size(data, timeDim) == nTime);
    
    origDelta = p.Results.origDelta;
    if isempty(origDelta)
        origDelta = nanmedian(diff(time));
    end
%     newDelta = p.Results.newDelta; 
    
    % figure out new time vector
    interpolateMethod = p.Results.interpolateMethod;
    binAlignmentMode = p.Results.binAlignmentMode;
%     if ~isempty(p.Results.binAlignmentModeOutput)
%         binAlignmentModeOutput = p.Results.binAlignmentModeOutput;
%     else
%         binAlignmentModeOutput = binAlignmentMode;
%     end
    
    if p.Results.extrapolate
        interpArgs = {'extrap'};
    else
        interpArgs = {};
    end
    
    nDimsOrig = ndims(data);
    timeNew = sampleAtTimes;
    
    % first shift the original timestamps to their centers
    % if causal bins, e.g. 20 ms, want to interp data to centers at
    % +10, +30 and then label as 20, 40. offset is -10
    interpInputTimes =  time + binAlignmentMode.getOffsetToBinCenter(origDelta);
%     interpOutputTimes = sampleAtTimes + binAlignmentModeOutput.getOffsetToBinCenter(newDelta);
    interpOutputTimes = timeNew;
    data = interp1(interpInputTimes, data, interpOutputTimes, interpolateMethod, interpArgs{:});
            
    assert(size(data, 1) == numel(timeNew));

    data = TensorUtils.unshiftdimToFirstDim(data, timeDim, nDimsOrig);
end