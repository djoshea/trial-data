function [dataCell, timeCell] = resampleDataCellInTime(dataCell, timeCell, varargin)
% this is a shorthand for calling resampleTensorInTime in a loop with some
% double checks to handle truncation when padding is used

    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
    p.addParameter('resampleMethod', 'filter', @isstringlike); % valid modes are filter, average, repeat , interp   
    p.addParameter('uniformlySampled', false, @islogical); % can speed things up if you know it's arleady uniform
    p.addParameter('progress', false, @slogical);
    % these are used as a secondary guard to truncate data within tMin :
    % tMax, when the input data includes padded edges to facilitate
    % resampling
    p.addParameter('tMinExcludingPadding', -Inf, @ismatrix);
    p.addParameter('tMaxExcludingPadding', Inf, @ismatrix);
    p.addParameter('expandToTimeMinMax', false, @islogical); % if false, time will ony be truncated to within tMin/tMax
    p.parse(varargin{:});

    
    progress = p.Results.progress;
    tMinExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMinExcludingPadding, size(dataCell));
    tMaxExcludingPadding = TensorUtils.singletonExpandToSize(p.Results.tMaxExcludingPadding, size(dataCell));
   
    timeCell = TensorUtils.singletonExpandToSize(timeCell, size(dataCell));
    tol = p.Results.timeDelta / 1000;
    
    if progress, prog = ProgressBar(numel(dataCell), 'Resampling data'); end
    for iD = 1:numel(dataCell)
        if isempty(dataCell{iD}) || size(dataCell{iD}, 1) == 1, continue; end % test for empty or single-sample

        if p.Results.expandToTimeMinMax
            extra_args = {'tMinOutput', tMinExcludingPadding(iD), 'tMaxOutput', tMaxExcludingPadding(iD)};
        else
            extra_args = {};
        end

        [d, t] = TrialDataUtilities.Data.resampleTensorInTime(dataCell{iD}, 1, timeCell{iD}, ...
            'origDelta', p.Results.origDelta, ...
            'timeDelta', p.Results.timeDelta, 'timeReference', p.Results.timeReference, ...
            'binAlignmentMode', p.Results.binAlignmentMode, ...
            'interpolateMethod', p.Results.interpolateMethod, ...
            'resampleMethod', p.Results.resampleMethod, ...
            'uniformlySampled', p.Results.uniformlySampled, ...
            extra_args{:}); % this is new, not sure if correct
        
        % if tmin and tmax are larger than the window.
        mask = t >= tMinExcludingPadding(iD) - tol & t <= tMaxExcludingPadding(iD) + tol;
        this_data = d(mask, :, :, :);
        this_time = t(mask);

        dataCell{iD} = this_data;
        timeCell{iD} = this_time;
        if progress, prog.update(iD); end
    end
    if progress, prog.finish(); end
end

