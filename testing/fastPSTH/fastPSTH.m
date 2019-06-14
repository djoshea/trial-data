function [cond_avg_CxTxU, n_trials_CxT, cond_sem_CxTxU, single_trial_RxTxU] = fastPSTH(spikes_RxU, conditions_R, t_start_delta_stop, filterFn, filterWindow, varargin)
% fancy version involving accum array, turns out to be quite slow and memory inefficient

% rates = fastPSTH(spikes_RxU, t_start_delta_stop, filterFn, filterWindow, varargin)
%
% dimensions: R trials, U units, T time bins, C conditions
%
% Inputs:
% 
%   spikes_RxU: cell of spike times by trial (R) by unit (U)
% 
%   t_start_delta_stop specifies the output time vector as [start delta stop] which will sample the rate at start:delta:stop
%
%   filterFn is a function f(tspike - tvecTimepoint) which returns the filter output for a spike at tspike (it MUST be vectorized and return output the same shape as the input)
% 
%   filterWindow a 2 vector (pre post) e.g. [-5 5] is the effective function width
%     where the function can be assumed to be zero outside that window

p = inputParser();
p.addParameter('t_start_stop_valid_R', [], @ismatrix);
p.addParameter('batch_size', 1e7, @islogical);
p.parse(varargin{:});

t_start_stop_valid_R = p.Results.t_start_stop_valid_R;
batch_size = p.Results.batch_size; % in spikes

R = size(spikes_RxU, 1);
U = size(spikes_RxU, 2);
C = max(conditions_R);

% combine all the spikes into 1 big list, keeping track of the trial and unit identity
[spikes_cat, sp_subs_ru] = catCellWhichSubs(1, spikes_RxU);
nsp = numel(spikes_cat);

tstart = t_start_delta_stop(1);
tdelta = t_start_delta_stop(2);
tstop = t_start_delta_stop(3);
tvec = tstart:tdelta:tstop;
T = numel(tvec);

filterWindowDeltas = [floor(filterWindow(1) * tdelta) ceil(filterWindow(2) * tdelta)];
F = filterWindowDeltas(2) - filterWindowDeltas(1) + 1;

% a spike a tspike needs to be evaluated at the first tvec(:) >= tspike - filterWindow(1) thru tvec(:) <= tspike + filterWindow(2)
if ~isempty(t_start_stop_valid_R)
    trial_valid = tvec >= t_start_stop_valid_R(:,1) + filterWindow(2) & tvec >= t_start_stop_valid_R(:, 2) + filterWindow(1);
else
    trial_valid = true(R, T);
end

% figure out the uniquified filter evals, where we evaluate it for the full filter window
% for each unique offset from the start of the binning
spike_bin = floor((spikes_cat - tstart + 1) ./ tdelta);
spike_filter_window_start = spike_bin + filterWindowDeltas(1); % 
delta_from_tvec_entries = spikes_cat - spike_bin ./ tdelta;

[unique_deltas, ~, which_delta_ind] = unique(delta_from_tvec_entries);
filter_tvec = filterWindowDeltas(1):filterWindowDeltas(2);
offsets = filter_tvec + unique_deltas;
filter_vals = single(filterFn(offsets)); 

[cond_sum, cond_ssq] = deal(zeros(C, T, U));

computeSingleTrials = true;
if computeSingleTrials
    % compute single trial rates, then aggregate into conditions, this requires more memory
    
    nbatch = ceil(nsp / batch_size);
    
    for ibatch = 1:nbatch
        fprintf('batch %d\n', ibatch);
        if ibatch == nbatch
            inds = nbatch*(ibatch-1) + 1 : nsp;
        else
            inds = (1:batch_size) + nbatch*(ibatch-1);
        end
        % generate 2-d matrix of filter_vals to match spikes_cat at each value of the filter
        vals = filter_vals(which_delta_ind(inds), :);

        % assemble trials x time x units subscripts to accumulate these into
        uind = repmat(sp_subs_ru(inds, 2), 1, F);
        tind = spike_filter_window_start(inds) + filter_tvec;
        rind = repmat(sp_subs_ru(inds, 1), 1, F);

        % remove filter timepoints that would land outside of tvec, accumulate filter coefficients into single trial rates
        within_bounds = tind >= 1 & tind <= T;
        vals = vals(within_bounds(:));
        subs = [rind(within_bounds(:)), tind(within_bounds(:)), uind(within_bounds(:))];
        single_trial = accumarray(subs, vals(:), [R T U]);

        % now accumulate for each condition the values and squared values
        [cind, tind, uind] = ndgrid(conditions_R, 1:T, 1:U);
        subs = [cind(:), tind(:), uind(:)];
        cond_sum = cond_sum + accumarray(subs, single_trial(:), [C T U]);
        cond_ssq = cond_ssq + accumarray(subs, single_trial(:).^2, [C T U]);
    end
    
else
    % skip holding onto single trial and just accumulate in a loop

    for iR = 1:R
        mask = sp_subs_ru(:, 1) == iR;
        vals = filter_vals(which_delta_ind(mask), :);
        
        % assemble time x units subscripts to accumulate this trial into
        uind = repmat(sp_subs_ru(mask, 2), 1, F);
        tind = spike_filter_window_start(mask) + filter_tvec;
        
        % remove filter timepoints that would land outside of tvec, accumulate filter coefficients into single trial rates
        within_bounds = tind >= 1 & tind <= T;    

        vals = vals(within_bounds(:));
        subs = [tind(within_bounds(:)), uind(within_bounds(:))];

        this_trial = shiftdim(accumarray(subs, vals(:), [T U]), -1); % 1 x T x U
        
        this_cond = conditions_R(iR);
        cond_sum(this_cond, :, :) = cond_sum(this_cond, :, :) + this_trial;
        cond_ssq(this_cond, :, :) = cond_ssq(this_cond, :, :) + this_trial.^2;
    end
    
    single_trial = [];
end

% count valid trials for each condition at each time bin (C x T)
[cind, tind] = ndgrid(conditions_R, 1:T);
subs = [cind(:), tind(:)];
n_trials = accumarray(subs, trial_valid(:), [C T]);

% compute condition averages
cond_avg = cond_sum ./ n_trials;

% compute standard error from sample standard deviation, num trials, and c_4 correction factor
% https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation#Estimating_the_standard_deviation_of_the_mean
ssd = sqrt( (cond_ssq - cond_sum.^2./n_trials) ./ (n_trials - 1) );
c4 = @(n) sqrt(2 ./ (n-1)) .* gamma(n/2) ./ gamma((n-1)/2); % correction factor for unbiased estimate
cond_sem = ssd ./ c4(n_trials) ./ sqrt(n_trials);

% rename outputs
cond_avg_CxTxU = cond_avg;
n_trials_CxT = n_trials;
cond_sem_CxTxU = cond_sem;
single_trial_RxTxU = single_trial;

end

function [out, subs] = catCellWhichSubs(dim, inCell)
    % works like cat, but returns a vector indicating which of the
    % inputs each element of out came from
    sz = size(inCell);
    out = cat(dim, inCell{:});
    nel = cellfun(@numel, inCell);
    ind = cell2mat(arrayfun(@(N, idx) idx*ones(N, 1), nel(:), (1:numel(nel))', 'UniformOutput', false));
    subs = TensorUtils.ind2subAsMat(sz, ind);
end

