function [data, timeNew] = resampleTensorInTime(data, timeDim, time, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp   
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x)); % min time for data valid
    p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x)); % max time for data valid

    p.addParameter('tMinOutput', [], @(x) isempty(x) || isscalar(x)); % min time for data that is output, if < tMin will be padded with NaNs 
    p.addParameter('tMaxOutput', [], @(x) isempty(x) || isscalar(x)); % max time for data that is output, if > tMax will be padded with NaNs
    
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
    binAlignmentMode = p.Results.binAlignmentMode;
   
    origDelta = p.Results.origDelta;
    if isempty(origDelta)
        origDelta = median(diff(time), 'omitnan');
    end
    if isempty(timeDelta)
        timeDelta = origDelta;
    end
    
    time = TrialDataUtilities.Data.removeSmallTimeErrors(time, timeDelta, timeReference);

    origDelta = round(origDelta, 5, 'significant');
    timeDelta = round(timeDelta, 5, 'significant');
    
    tMin = p.Results.tMin;
    tMax = p.Results.tMax;
    [tMinCalc, tMaxCalc] = p.Results.binAlignmentMode.getTimeLimitsForRebinning(min(time), max(time), origDelta, timeDelta, timeReference);
    if isempty(tMin)
        tMin = tMinCalc;
    else
        tMin = TrialDataUtilities.Data.removeSmallTimeErrors(tMin, timeDelta, timeReference);
    end
    if isempty(tMax)
        tMax = tMaxCalc; 
    else
        tMax = TrialDataUtilities.Data.removeSmallTimeErrors(tMax, timeDelta, timeReference);
    end
    
    timeNew = (tMin:timeDelta:tMax)';
    nDimsOrig = ndims(data);
    
    if origDelta == 0 && tMin == tMax
        % single sample special case
        return;
    end
    
    data = TensorUtils.shiftdimToFirstDim(data, timeDim);
    deltaIsChanging = ~TrialDataUtilities.Stats.isequaltol(timeDelta, origDelta, origDelta / 1000);
    
    % do this to avoid off by one errors when generating time vectors
    if ~deltaIsChanging
        origDelta = timeDelta;
    end
    
    % build time vector for the original that starts at the appropriate
    % tMin so that we end up with the right samples
    timeUniform = (tMin:origDelta:tMax)';
    
    origClass = class(data);
    if ~ismember(class(data), {'single', 'double'})
        data = single(data);
    end
    
    switch p.Results.resampleMethod
        case 'filter'
            if ~p.Results.uniformlySampled
                % sample to uniform grid
                data = interp1(time, data, timeUniform, interpolateMethod, 'extrap');
                time = timeUniform;
            end

            if deltaIsChanging
                % use resampling
                [data] = TrialDataUtilities.Data.resamplePadEdges(data, time, timeNew, p.Results.binAlignmentMode, interpolateMethod);   
            end

            assert(size(data, 1) == numel(timeNew));
            
        case 'repeat'
            % TODO NEED TO TAKE BINALIGNMENTMODE INTO ACCOUNT
            
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta < origDelta
                % upsample via repelem
                P = origDelta/timeDelta;
                assert(ceil(P) == P, 'New sampling interval must be an integer multiple of old sampling interval');
                data = repelem(data, P, 1);
                
                % trim extra copies at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif ~deltaIsChanging
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
            
        case 'average'
            % TODO NEED TO TAKE BINALIGNMENTMODE INTO ACCOUNT
            
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta > origDelta
                % downsample via averaging
                P = timeDelta/origDelta;
                assert(TrialDataUtilities.Data.isequaltol(round(P), P, P / 1000), 'New sampling interval must be an integer multiple of old sampling interval');
                P = ceil(P);
                data = blockproc(data, [P 1], @(block) mean(block.data));
                
                % trim average from partial blocks at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif ~deltaIsChanging
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
                
        case 'interp'
            % first shift the original timestamps to their centers
            % if causal bins, e.g. 20 ms, want to interp data to centers at
            % +10, +30 and then label as 20, 40. offset is -10
            
            interpInputTimes =  time + binAlignmentMode.getOffsetToBinCenter(origDelta);
            interpOutputTimes = timeNew + binAlignmentMode.getOffsetToBinCenter(timeDelta);
            data = interp1(interpInputTimes, data, interpOutputTimes,  'linear', 'extrap');
            
%             interpToTime = timeNew + p.Results.binAlignmentMode.getOffsetToBinCenter(timeDelta);
%             data = interp1(time, data, interpToTime, interpolateMethod);
%
        case 'nearest'
            interpInputTimes =  time + binAlignmentMode.getOffsetToBinCenter(origDelta);
            interpOutputTimes = timeNew + binAlignmentMode.getOffsetToBinCenter(timeDelta);
            data = interp1(interpInputTimes, data, interpOutputTimes, 'nearest', 'extrap');

        otherwise
            error('Unknown resampleMethod %s', p.Results.resampleMethod)
    end
    
    if ~strcmp(origClass, class(data))
        data = cast(data, origClass);
    end

    tMinOutput = p.Results.tMinOutput;
    tMaxOutput = p.Results.tMaxOutput;
    output_specified = false;
    if isempty(tMinOutput) || isnan(tMinOutput)
        tMinOutput = tMin;
    else
        output_specified = true;
    end
    if isempty(tMaxOutput) || isnan(tMaxOutput)
        tMaxOutput = tMax;
    else
        output_specified = true;
    end
    if output_specified
        tMinOutput = TrialDataUtilities.Data.removeSmallTimeErrors(tMinOutput, timeDelta, timeReference);
        tMaxOutput = TrialDataUtilities.Data.removeSmallTimeErrors(tMaxOutput, timeDelta, timeReference);
%         assert(tMinOutput <= tMin);
%         assert(tMaxOutput >= tMax);

        % There are two distinct time grids in play here, and they are decided by different rules:
        %
        %   timeNew (built above from tMin/tMax, which come from min(time)/max(time) via
        %     getTimeLimitsForRebinning) is the VALIDITY grid: the output bins that the actual input
        %     samples can fully support. The origDelta inset applied there is essential -- e.g. a
        %     Centered 10 ms bin labeled t needs coverage over t +/- 5, but with 1 ms input samples
        %     the outermost sample only certifies +/- 0.5, so the first supportable label lies
        %     newDelta/2 - origDelta/2 inside the data extent.
        %
        %   timeOutput (here) is the FRAMING grid: which labels the caller asked us to return.
        %     tMinOutput/tMaxOutput arrive as desired output sample labels already on the timeDelta
        %     grid (e.g. getAnalog passes the alignment window from getTimeStartStopEachTrial), NOT
        %     as a data extent to be probed for support.
        %
        % So we must NOT push tMinOutput/tMaxOutput through getTimeLimitsForRebinning: that would
        % apply the data-extent inset a second time, to a quantity that is not a data extent. It
        % silently drops one bin from each edge of every window whenever origDelta < timeDelta, and
        % it inverts (-> NaN, -> a crash below) any window narrower than newDelta - origDelta. E.g.
        % a requested window of [810 820] with origDelta = 1, timeDelta = 10, Centered yields
        % tMin = ceil((810 - 0.5 + 5)/10)*10 = 820 > tMax = floor((820 + 0.5 - 5)/10)*10 = 810.
        %
        % Instead we only snap the requested labels onto the bin lattice (multiples of timeDelta
        % offset from timeReference) so that timeOutput and timeNew line up bin-for-bin for the
        % assignment below. Rounding is inward so we never widen beyond what was requested; for a
        % window already on the grid (essentially always) this is a no-op.
        %
        % Validity is still enforced, just in the right place: data is preallocated all-NaN over
        % timeOutput and only the labels present in timeNew are written into it. Requesting a frame
        % wider than the input samples support therefore yields NaN at the unsupported edges rather
        % than partial or extrapolated values -- which is exactly the intent under
        % AlignDescriptor.IGNORE, and matches the NaN edges the filtered spike path produces.
        tol = timeDelta / 1000;
        tMinOutput = timeReference + TrialDataUtilities.Stats.ceiltol((tMinOutput - timeReference) / timeDelta, tol) * timeDelta;
        tMaxOutput = timeReference + TrialDataUtilities.Stats.floortol((tMaxOutput - timeReference) / timeDelta, tol) * timeDelta;

        if isnan(tMinOutput) || isnan(tMaxOutput) || tMinOutput > tMaxOutput
            % requested window does not contain a single label on the bin grid: no output samples
            data = data(1:0, :, :, :, :, :);
            timeNew = zeros(0, 1);
        else
            timeOutput = (tMinOutput:timeDelta:tMaxOutput)';

            data_preExpand = data;
            sz_out = size(data);
            sz_out(1) = numel(timeOutput);
            data = nan(sz_out, like=data_preExpand);

            mask_take = timeNew >= tMinOutput - tol & timeNew <= tMaxOutput + tol;
            data_preExpand = data_preExpand(mask_take, :, :, :, :);
            timeNew = timeNew(mask_take);

            % timeNew may be empty if the supportable bins lie entirely outside the requested
            % window, in which case the all-NaN preallocation above is the correct output
            if ~isempty(timeNew)
                mask_assign = timeOutput >= timeNew(1) - tol & timeOutput <= timeNew(end) + tol;
                assert(nnz(mask_assign) == numel(timeNew), 'Output time grid does not align with resampled time grid');
                data(mask_assign, :, :, :, :, :) = data_preExpand;
            end
            timeNew = timeOutput;
        end
    end

    data = TensorUtils.unshiftdimToFirstDim(data, timeDim, nDimsOrig);
end