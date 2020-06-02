function y = varWeightedCorrected(x, w, dim, asFrequencyWeights)
    % varWeightedCorrected(x, w, dim, asFrequencyWeights)
    % corrected version of the weighted variance.
    % if asFrequencyWeights is true: assumes weights represent integer number of cases [default]
    % if asFrequencyWeights is true: assumes weights represent reliability
    % always ignores NaNs as missing data. 
    %
    % c.f. https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    
    if nargin < 4
        asFrequencyWeights = true; 
    end
    
    % make w same size as x
    w = abs(w);
    if isvector(w)
        if size(w, dim) ~= size(x, dim)
            w = TensorUtils.orientVectorAlongDim(w, dim);
        end
        w = w .* ones(size(x));
    end
    assert(isequal(size(x), size(w)));

    % Count up non-NaN weights at non-NaN elements
    mask = isnan(x) | isnan(x);
    w(mask) = 0;
    x(mask) = 0;

    mu = sum(w .* x, dim, 'omitnan') ./ sum(w, dim, 'omitnan');
    numer = sum(w.*(x-mu).^2, dim, 'omitnan');
          
    if asFrequencyWeights
        denom = sum(w, dim) - 1;
    else
        % as reliability weights
        v1 = sum(w, dim, 'omitnan');
        v2 = sum(w.^2, dim, 'omitnan');
        denom = v1 - (v2 ./ v1);
    end
    y = numer ./ denom;
end