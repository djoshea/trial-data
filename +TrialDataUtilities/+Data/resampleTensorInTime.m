function [data, timeNew] = resampleTensorInTime(data, timeDim, time, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.parse(varargin{:});
    
    assert(isvector(time));
    nTime = numel(time);
    assert(size(data, timeDim) == nTime);
    
    % figure out new time vector
    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    interpolateMethod = p.Results.interpolateMethod;
   
    tMinRaw = min(time);
    tMaxRaw = max(time);
    if isempty(timeDelta)
        timeDelta = nanmedian(diff(x), time);
    end

    tMin = timeReference + floor((tMinRaw - timeReference) / timeDelta) * timeDelta;
    tMax = timeReference + ceil((tMaxRaw - timeReference) / timeDelta) * timeDelta;
   
    timeNew = makecol(tMin:timeDelta:tMax);
    nDimsOrig = ndims(data);
    data = TensorUtils.shiftdimToFirstDim(data, timeDim);
    data = interp1(time, data, timeNew, interpolateMethod);
    data = TensorUtils.unshiftdimToFirstDim(data, timeDim, nDimsOrig);
    
end