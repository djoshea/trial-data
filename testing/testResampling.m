
N = 5;
c = 2;

psetR = pset.filterBases(1:N);
dm = psetR.dataMean{1};

drSet = psetR.computeDataMeanResampleTrialsWithinConditions(100, 0);
dr = drSet{1};
drperc = quantile(dr, [0.05 0.95], 4);

%%
i = 5;

clf
plot(squeeze(dr(i, c, :, :)), '-', 'Color', [0.5 0.5 0.5]);
hold on
plot(squeeze(dm(i, c, :))', 'b-', 'LineWidth', 2);
plot(squeeze(drperc(i, c, :, :)), 'r-', 'LineWidth', 2);

