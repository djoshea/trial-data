function [delays, aligned] = findDelaysMultiAnalogTensor(data,varargin)
% Estimates delays between sets of simultaneously sampled signals.
%   delays = finddelaysMultiAnalogTensor(x,varargin)
%
%   inputs: 
%     X is a 3 dimensional tensor as returned by .getMultiAnalogAsTensor
%       dims are trials (R) x time (T) x channels (C)
%
%   returns:
%     delays is R x 1 set of delays used to slide each T x C matrix to align
%
%   parameters:
%     alignTo: T x C matrix to align each trace to. By default is
%       nanmean over trials
%
%     maxLag: maximum shift to detect, default T/2
%
%   also parameters for removeDelaysMultiAnalogTensor:
%     fillMode: how to fill the "exposed" edges of the signal matrices,
%       either a scalar like NaN or 0 or 'hold' to extend the edge values
%
%   based off Matlab's finddelay

    assert(ndims(data) <= 3);

    % R x T x C --> T x C x R
    x = permute(data, [2 3 1]);

    T = size(x, 1);
    C = size(x, 2);
    R = size(x, 3);

    p = inputParser();
    p.addOptional('maxLag', floor(T / 2), @isscalar);
    p.addParameter('alignTo', mean(x, 3, 'omitnan'), @isnumeric);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    y = p.Results.alignTo;
    maxLag = p.Results.maxLag;

    % The largest maximum window size determines the size of the 
    % cross-correlation vector/matrix c.
    % Preallocate normalized cross-correlation vector/matrix c.
    c_normalized = zeros(2*maxLag+1,C,R);
    lags = -maxLag:maxLag;
    nLags = size(c_normalized, 1);
    index_max = zeros(1,1,R);

    % Compute absolute values of normalized cross-correlations between x and
    % all columns of y: function XCORR does not take into account special case
    % when either x or y is all zeros, so we don't use its normalization option
    % 'coeff'. Values of normalized cross-correlations computed for a lag of
    % zero are stored in the middle row of c at index i = max_maxlag+1 (c has
    % an odd number of rows).

    for r = 1:R
        for c = 1:C
            c_normalized(:,c,r) = absnanxcorrNorm(x(:,c,r), y(:,c), maxLag);
        end
    end

    % sum the cross correlations across channels
    c_normalized_sum = TensorUtils.squeezeDims(sum(c_normalized, 2), 2); % nLags x R;

    % Find indices of lags resulting in the largest absolute values of
    % normalized cross-correlations: to deal with periodic signals, seek the
    % lowest (in absolute value) lag giving the largest absolute value of
    % normalized cross-correlation.
    % Find lowest positive or zero indices of lags (negative delays) giving the
    % largest absolute values of normalized cross-correlations. 
    [max_c_pos,index_max_pos] = max(c_normalized_sum(maxLag+1:end,:),[],1);    
    % Find lowest negative indices of lags (positive delays) giving the largest
    % absolute values of normalized cross-correlations. 
    [max_c_neg,index_max_neg] = max(flipud(c_normalized_sum(1:maxLag,:)),[],1);

    max_c = nan(R, 1);

    if isempty(max_c_neg)
        % Case where MAXLAG is all zeros.
        index_max = maxLag + index_max_pos;
    else
        for r=1:R
            if max_c_pos(r)>max_c_neg(r)
                % The estimated lag is positive or zero.
                index_max(r) = maxLag + index_max_pos(r);
                max_c(r) = max_c_pos(r);
            elseif max_c_pos(r)<max_c_neg(r)
                % The estimated lag is negative.
                index_max(r) = maxLag + 1 - index_max_neg(r);
                max_c(r) = max_c_neg(r);
            elseif max_c_pos(r)==max_c_neg(r)
                if index_max_pos(r)<=index_max_neg(r)
                    % The estimated lag is positive or zero.
                    index_max(r) = max_maxlag + index_max_pos(r);
                    max_c(r) = max_c_pos(r);
                else
                    % The estimated lag is negative.
                    index_max(r) = max_maxlag + 1 - index_max_neg(r);
                    max_c(r) = max_c_neg(r);
                end 
            end   
        end
    end

    delays = makecol(lags(index_max));

    % Set to zeros estimated delays for which the normalized cross-correlation
    % values are below a given threshold (spurious peaks due to FFT roundoff
    % errors).
    for i=1:R
        if max_c(r)<1e-8
            delays(r) = 0;
            if maxlag(r)~=0
                warning(message('signal:finddelay:noSignificantCorrelationVector', r));
            end
        end    
    end

    if nargout > 1
        aligned = TrialDataUtilities.Data.removeDelaysMultiAnalogTensor(data, delays, p.Unmatched);
    end

end

function c = absnanxcorrNorm(x, y, maxlag)
    cxx = nansum(abs(x).^2);
    cyy = nansum(abs(y).^2);

    x = (x - nanmean(x)) / nanstd(x);
    x(isnan(x)) = 0;
    
    c = xcorr(x, y, maxlag);
    
    if cxx==0 || cyy==0
        c(:) = 0;
    else
        c = c / sqrt(cxx*cyy);
    end
end
