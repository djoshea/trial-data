function [data, timeNew] = resampleTensorInTime(data, timeDim, time, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('resampleMethod', 'filter', @ischar); % valid modes are filter, average, repeat , interp   
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('uniformlySampled', false, @islogical); % can speed things up if you know it's arleady uniform
    p.parse(varargin{:});
    
    if isempty(data)
        timeNew = [];
        return;
    end
    
    assert(isvector(time));
    nTime = numel(time);
    assert(size(data, timeDim) == nTime);
    
    % figure out new time vector
    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    interpolateMethod = p.Results.interpolateMethod;
   
    origDelta = p.Results.origDelta;
    if isempty(origDelta)
        origDelta = nanmedian(diff(time));
    end
    if isempty(timeDelta)
        timeDelta = origDelta;
    end

    tMin = p.Results.tMin;
    tMax = p.Results.tMax;
    [tMinCalc, tMaxCalc] = p.Results.binAlignmentMode.getTimeLimitsForRebinning(min(time), max(time), origDelta, timeDelta, timeReference);
    if isempty(tMin), tMin = tMinCalc; end
    if isempty(tMax), tMax = tMaxCalc; end
    
    timeNew = (tMin:timeDelta:tMax)';
    nDimsOrig = ndims(data);
    data = TensorUtils.shiftdimToFirstDim(data, timeDim);
    
    % build time vector for the original that starts at the appropriate
    % tMin so that we end up with the right samples
    timeUniform = (tMin:origDelta:tMax)';
    
    switch p.Results.resampleMethod
        case 'filter'
            if timeDelta == origDelta
                data = interp1(time, data, timeNew, interpolateMethod);
            else
                % use resampling
                [data, ty] = TrialDataUtilities.Data.resamplePadEdges(data, time, p.Results.timeReference, origDelta, timeDelta, ...
                    interpolateMethod, p.Results.binAlignmentMode, p.Results.uniformlySampled);
                mask = ty >= tMin & ty <= tMax;
                data = data(mask, :, :, :);
                assert(size(data, 1) == numel(timeNew));
            end
            
        case 'repeat'
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta < origDelta
                % upsample via repelem
                P = origDelta/timeDelta;
                assert(ceil(P) == P, 'New sampling interval must be an integer multiple of old sampling interval');
                data = repelem(data, P, 1);
                
                % trim extra copies at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif timeDelta == origDelta
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
            
        case 'average'
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta > origDelta
                % downsample via averaging
                P = timeDelta/origDelta;
                assert(ceil(P) == P, 'New sampling interval must be an integer multiple of old sampling interval');
                data = blockproc(data, [P 1], @(block) mean(block.data));
                
                % trim average from partial blocks at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif timeDelta == origDelta
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
                
        case 'interp'
            data = interp1(time, data, timeNew, interpolateMethod);
            
        otherwise
            error('Unknown resampleMethod %s', p.Results.resampleMethod)
            
    end
    
    data = TensorUtils.unshiftdimToFirstDim(data, timeDim, nDimsOrig);
end