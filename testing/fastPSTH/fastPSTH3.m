function [cond_avg_CxTxU, n_trials_CxT, cond_sem_CxTxU, single_trial_RxTxU] = fastPSTH3(spikes_RxU, conditions_R, ...
    t_start_delta_stop, filterFn, filterWindow, varargin)
% even faster version, simply iterates once over spikes and keeps track of which unit, trial we're on
%
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
p.addParameter('trial_start_stop_valid', [], @ismatrix);
p.addParameter('quiet', false, @islogical);
p.parse(varargin{:});

t_start_stop_valid_R = p.Results.trial_start_stop_valid;
quiet = p.Results.quiet;

R = size(spikes_RxU, 1);
U = size(spikes_RxU, 2);
C = max(conditions_R);

if ~quiet
    prog = ProgressBar(R, 'Computing PSTHs for %d trials', R);
end

% combine all the spikes into 1 big list, keeping track of the trial and unit identity
spikes_UxR = spikes_RxU';
clear spikes_RxU;
spikes_cat = cat(1, spikes_UxR{:});
nsp_UxR = zeros(R, U, 'uint32');
for i = 1:numel(spikes_UxR)
    nsp_UxR(i) = uint32(numel(spikes_UxR{i}));
end
nsp = sum(nsp_UxR(:), 'native');

[Uind, Rind] = ndgrid(1:U, 1:R);
sp_ind_start_UxR = cumsum(nsp_UxR(:));

tstart = t_start_delta_stop(1);
tdelta = t_start_delta_stop(2);
tstop = t_start_delta_stop(3);
tvec = tstart:tdelta:tstop;
T = numel(tvec);
gain = single(1000 / tdelta); % convert from ms to spikes/sec

% a spike a tspike needs to be evaluated at the first tvec(:) >= tspike - filterWindow(1) thru tvec(:) <= tspike + filterWindow(2)
if ~isempty(t_start_stop_valid_R)
    trial_valid = tvec >= t_start_stop_valid_R(:,1) + filterWindow(2) & tvec >= t_start_stop_valid_R(:, 2) + filterWindow(1);
else
    trial_valid = true(R, T);
end

% figure out the uniquified filter evals, where we evaluate it for the full filter window
% for each unique offset from the start of the binning
filterWindowDeltas = [floor(filterWindow(1) * tdelta) ceil(filterWindow(2) * tdelta)];
spike_bin = floor((spikes_cat - tstart + 1) ./ tdelta);
spike_filter_window_start = spike_bin + filterWindowDeltas(1);
delta_from_tvec_entries = spikes_cat - spike_bin ./ tdelta;

[unique_deltas, ~, which_delta_ind] = uniquetol(delta_from_tvec_entries, tdelta / 100);
filter_tvec = filterWindowDeltas(1):filterWindowDeltas(2);
offsets = filter_tvec + unique_deltas;
filter_vals = single(filterFn(offsets) .* gain);
F = numel(filter_tvec);

[cond_sum, cond_ssq] = deal(zeros(C, T, U, 'single'));
single_trial = deal(zeros(R, T, U, 'single'));
this_trial = zeros(1, T, U, 'single');

u_r_ind = 1; % u, r linear ind
r = 1;
u = 1;
for i = 1:nsp % itrerate over spikes
    for f = 1:F
        tind = spike_filter_window_start(i) + f - 1;
        % remove filter timepoints that would land outside of tvec, accumulate filter coefficients into single trial rates
        if tind >= 1 && tind <= T
            this_trial(1, tind, u) = this_trial(1, tind, u) + filter_vals(which_delta_ind(i), f);
        end
    end
    
    if u_r_ind < U*R && i >= sp_ind_start_UxR(u_r_ind+1)
        % move to the next trial or unit within trial
        u_r_ind = u_r_ind + 1;
        u = Uind(u_r_ind);
        newr = Rind(u_r_ind);
        
        if newr > r
            % moving to new trial, accumulate this trial into condition sums 
            
            if ~quiet, prog.update(r); end
            this_cond = conditions_R(r);
            single_trial(r, :, :) = this_trial;
            cond_sum(this_cond, :, :) = cond_sum(this_cond, :, :) + this_trial;
            cond_ssq(this_cond, :, :) = cond_ssq(this_cond, :, :) + this_trial.^2;
            this_trial(:) = 0;
        end
        r = newr;
    end
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

if ~quiet, prog.finish(); end

end

