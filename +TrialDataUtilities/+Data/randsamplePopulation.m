function [samples, idx] = randsamplePopulation(varargin)
% works like Matlab's randsample, with a few simplifications:
% - always treats the population as a population, even if scalar, no automatic translation to 1:n
% - reduces the number of samples to numel(population) if sampling without replacement
% - no weighting supported

if isa(varargin{1}, 'RandStream')
    args = varargin(2:end);
    s = varargin{1};
else
    s = [];
    args = varargin;
end

p = inputParser();
p.addRequired('population', @isvector);
p.addRequired('numSamples', @isscalar);
p.addOptional('replace', false, @islogical);
p.parse(args{:});

population = p.Results.population;
n = numel(population);
k = p.Results.numSamples;
replace = p.Results.replace;

if replace
    if isempty(s)
        idx = randi(n, k, 1);
    else
        idx = randi(s, n, k, 1);
    end
else
    if k > n
        % too many samples requested, just shuffle the full population
        % (i.e. set k == n)
        if isempty(s)
            idx = randperm(n);
        else
            idx = randperm(s,n);
        end
    else
        % this is taken directly from randsample code:
        % If the sample is a sizable fraction of the population,
        % just randomize the whole population (which involves a full
        % sort of n random values), and take the first k.
        if 4*k > n
            if isempty(s)
                rp = randperm(n);
            else
                rp = randperm(s,n);
            end
            idx = rp(1:k);

        % If the sample is a small fraction of the population, a full sort
        % is wasteful.  Repeatedly sample with replacement until there are
        % k unique values.
        else
            x = zeros(1,n); % flags
            sumx = 0;
            while sumx < floor(k) % prevent infinite loop when 0<k<1
                if isempty(s)
                    x(randi(n,1,k-sumx)) = 1; % sample w/replacement
                else
                    x(randi(s,n,1,k-sumx)) = 1; % sample w/replacement
                end
                sumx = sum(x); % count how many unique elements so far
            end
            idx = find(x > 0);
            if isempty(s)
                idx = idx(randperm(k));
            else
                idx = idx(randperm(s,k));
            end
        end
    end
end

samples = population(idx);

    