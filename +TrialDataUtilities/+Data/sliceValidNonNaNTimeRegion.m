function [varargout] = sliceValidNonNaNTimeRegion(varargin)
% [out1, out2, ..., maskOrMaskCell] = sliceValidNonNaNTimeRegion(out1, out2, )
% Given a set of either:
%   1. cell arrays with the same length (N) containing vectors
%        of the same lengths across each argument (T_i) 
%   2. matrices with the same size along (N x T x D) or be vectors with length
%      T consistent across arguments
% Find the time points where at least one non-nan value is present across
% all of the inputs, slice each input to that time region, and return those
% sliced matrices.

isCell = cellfun(@iscell, varargin);

if all(isCell)
    useCell = true;
    
    % check sizes
    n = cellfun(@numel, varargin);
    assert(all(n == n(1)), 'All cell inputs must have same numel');
    
elseif all(~isCell)
    useCell = false;
    
    isvec = cellfun(@isvector, varargin);
    n = cellfun(@(x) size(x, 1), varargin);
    t = cellfun(@(x) size(x, 2), varargin);
    nel = cellfun(@numel, varargin);
    t(isvec) = nel(isvec);
    
    if any(~isvec)
        assert(all(n(isvec) == max(n(isvec))), 'All matrix inputs must have same numrows');
    end
    assert(all(t == t(1)), 'All matrix and vector inputs must have same number of timepoints');

else
    error('All inputs must be cells or matrices');
end
   
varargout = cell(numel(varargin), 1);
nArg = numel(varargin);

if useCell
    nTrials = numel(varargin{1});
    for iArg = 1:nArg+1
        varargout{iArg} = cell(nTrials, 1);
    end
    
    for iTrial = 1:numel(varargin{1})
        % find the valid region for this trial
        nanmask = isnan(varargin{1}(iTrial));
        for iArg = 2:numel(varargin)
            nanmask = nanmask & isnan(varargin{iArg}{iTrial});
        end
        i1 = find(~nanmask, 1, 'first');
        i2 = find(~nanmask, 1, 'last');
        
        % select the valid region for this trial
        if ~isempty(i1) && ~isempty(i2)
            for iArg = 1:nArg
                varargout{iArg}{iTrial} = varargin{iArg}{iTrial}(i1:i2);
            end
        end
        
        m = falsevec(size(nanmask));
        m(i1:i2) = true;
        varargout{nArg+1} = m;
    end
else
    if any(~isvec)
        indMat = find(~isvec, 1);
        nTrials = size(varargin{indMat}, 1);
        nTime = size(varargin{indMat}, 2);
    else
        nTrials = 1;
        nTime = numel(varargin{1});
    end
    
    % find the valid window
    keepmask = true(nTime, 1);
    for iArg = 1:nArg
%         if isvec(iArg)
%             keepmask = keepmask & makecol(~isnan(varargin{iArg}));
%         else
            if ndims(varargin{iArg}) == 3
                validPortion = all(~isnan(varargin{iArg}), 3);
            else
                validPortion = ~isnan(varargin{iArg});
            end
            keepmask = keepmask & any(validPortion, 1)';
%         end
    end
    
    i1 = find(keepmask, 1, 'first');
    i2 = find(keepmask, 1, 'last');
    
    % select the valid region
    if ~isempty(i1) && ~isempty(i2)
        for iArg = 1:nArg
            if isvec(iArg)
                varargout{iArg} = varargin{iArg}(i1:i2);
            else
                varargout{iArg} = varargin{iArg}(:, i1:i2, :);
            end
        end
    else
        % no valid window found
        warning('No valid window found across matrices');
        for iArg = 1:nArg
            if isvec(iArg)
                varargout{iArg} = [];
            else
                varargout{iArg} = nan(nTrials, 0);
            end
        end
    end
    
    m = false(nTime, 1);
    m(i1:i2) = true;
    varargout{nArg+1} = m;
end