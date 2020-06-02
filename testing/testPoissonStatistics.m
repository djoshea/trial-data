n = 1000;
r = 20;
p = r / n;

trials = 20;

rng(1);
spike_occurred = rand(trials, n) < p;
times = arrayfun(@(r) find(spike_occurred(r, :))', (1:trials)', 'UniformOutput', false);

%%
trials = size(times, 1);
td_test = TrialData.buildEmptyWithTrialDurations(repmat(n, trials, 1));

td_test = TrialDataConditionAlign(td_test);
td_test = td_test.addSpikeChannel('unit1_1', times);

%%

sf = NonOverlappingCausalSpikeBinFilter(10);

[m, tvec, se, st] = td_test.getPSTH('unit1_1', 'spikeFilter', sf);

[~, ~, sep, stp] = td_test.getPSTH('unit1_1', 'spikeFilter', sf, 'assumePoissonStatistics', true);

clf;
plot(tvec, st, 'k-');
hold on;
plot(tvec, stp, 'r-');
hold off;

%%

global_mean = mean(x(:, 50:70), 'all');

var(x(:, 50:70) / 100, 0, 'all') 

clf;
histogram(x(:, 50:70), 0:300, 'normalization', 'pdf')
hold on;

plot(0:300, poisspdf(0:300, global_mean), 'k-');

