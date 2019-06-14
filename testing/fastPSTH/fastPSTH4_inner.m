function [single_trial, cond_sum, cond_ssq] = fastPSTH4_inner(R, C, T, U, nsp_UxR, filter_vals, which_delta_ind, spike_filter_window_start, conditions_R)

R = uint32(R);
C = uint32(C);
T = uint32(T);
U = uint32(U);
filter_vals = single(filter_vals);
which_delta_ind = single(which_delta_ind);
spike_filter_window_start = uint32(spike_filter_window_start);

nsp_UxR = uint32(nsp_UxR);
nsp = sum(nsp_UxR(:), 'native');
assert(uint32(numel(which_delta_ind)) == nsp);
assert(uint32(numel(conditions_R)) == R);

[Uind, Rind] = ndgrid(1:U, 1:R);
sp_ind_start_UxR = cumsum(nsp_UxR(:));

F = size(filter_vals, 2);
[cond_sum, cond_ssq] = deal(zeros(C, T, U, 'single'));
single_trial = deal(zeros(R, T, U, 'single'));

u_r_ind = uint32(1); % u, r linear ind
r = uint32(1);
u = uint32(1);

for i = 1:nsp % itrerate over spikes
    vals = filter_vals(which_delta_ind(i), :);
    tind = spike_filter_window_start(i) + (uint32(0):uint32(F-1));

     % remove filter timepoints that would land outside of tvec, accumulate filter coefficients into single trial rates
    within_bounds = tind >= uint32(1) & tind <= uint32(T);    

    single_trial(r, tind(within_bounds), u) = single_trial(r, tind(within_bounds), u) + vals(within_bounds);

    if u_r_ind < U*R && i >= sp_ind_start_UxR(u_r_ind+1)
        % move to the next trial or unit within trial
        u_r_ind = u_r_ind + 1;
        u = Uind(u_r_ind);
        newr = Rind(u_r_ind);
        
        if newr > r
%             fprintf('trial %u / %u\r', r, R);
            % moving to new trial, accumulate this trial into condition sums 
            this_cond = conditions_R(r);
            cond_sum(this_cond, :, :) = cond_sum(this_cond, :, :) + single_trial(r, :, :);
            cond_ssq(this_cond, :, :) = cond_ssq(this_cond, :, :) + single_trial(r, :, :).^2;
            r = newr;
        end
    end
end

end