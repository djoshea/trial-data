%% testOutsideOfTrialModes
% Characterize what the spike filter / analog resampling return at OUT-OF-TRIAL edges under each
% AlignDescriptor.outsideOfTrialMode (TRUNCATE / IGNORE / INVALIDATE).
%
% Minimal synthetic td: 1 trial, TrialStart=0, TrialEnd=100 ms, an event `start` at t=0, and one
% spike channel. Aligned to [start-100, start+100] so the window extends 100 ms BEFORE the trial
% (i.e. [-100, 0) is out-of-trial). Filtered with a 10 ms rectangular (causal) filter.
%
% Expected:
%   IGNORE     -> full [-100,+100] window kept; out-of-trial [-100,0) holds valid 0 rates (NOT NaN)
%   TRUNCATE   -> window clipped toward the trial extent [0,100] (also loses the causal edge)
%   INVALIDATE -> trial dropped (padded window exceeds the trial)
%
% This grounds the premise behind the PopulationTrajectorySet IGNORE fix: spike bins are a valid 0
% at out-of-trial edges (so the neural window is already the full nominal window), whereas analog
% resampling returns NaN there (so its window gets trimmed unless IGNORE is honored at the pset
% per-trial-limits stage).

td = TrialData.buildEmptyWithTrialStartTrialEnd(0, 100);          % 1 trial, [0,100] ms
td = TrialDataConditionAlign(td);
td = td.addSpikeChannel('unit1_1', {[10;20;30;55;90]});          % nTrials x 1 cell of spike times (ms)
td = td.addEvent('start', 0);                                    % event at t=0
td = td.start('start', -100).stop('start', +100).zero('start'); % window [-100, +100]

sf = RectangularCausalSpikeFilter(10);                           % 10 ms boxcar (causal)

%% Spike rates across the three modes
for mode = ["Truncate", "Ignore", "Invalidate"]
    tdm = td.("setOutsideOfTrial" + mode)();
    [rates, tvec] = tdm.getSpikeRateFilteredAsMatrix('unit1_1', 'spikeFilter', sf);
    fprintf('--- spikes %s --- tvec[%g..%g] n=%d\n', mode, min(tvec), max(tvec), numel(tvec));
    disp([tvec(:)'; rates]);
end

%% Analog counterpart (out-of-trial edge should come back as NaN, not 0)
% The spike-vs-analog asymmetry the PopulationTrajectorySet fix reconciles. NOTE: adjust the
% addAnalogChannel call to match its signature on your build if this errors.
% tda = td.addAnalogChannel('ana1', (0:100)', sin((0:100)'/10));   % samples over the trial [0,100]
% for mode = ["Truncate", "Ignore", "Invalidate"]
%     tdm = tda.("setOutsideOfTrial" + mode)();
%     [mat, tvec] = tdm.getAnalogAsMatrix('ana1', 'timeDelta', sf.timeDelta, 'binAlignmentMode', sf.binAlignmentMode);
%     fprintf('--- analog %s --- tvec[%g..%g] n=%d nNaN=%d\n', mode, min(tvec), max(tvec), numel(tvec), nnz(isnan(mat)));
%     disp([tvec(:)'; mat]);
% end

%% Pset-level regression check (run BEFORE and AFTER the PopulationTrajectorySet IGNORE fix)
% Build an IGNORE pset from a single td and inspect tvecDataMean{1}. After the fix, an ANALOG pset
% keeps the full nominal window (with NaN edge bins); a SPIKE pset is unchanged (already valid-0
% edges). Uncomment once the analog channel above is wired.
% psS = PopulationTrajectorySetBuilder.fromMultipleTrialData({td.setOutsideOfTrialIgnore()}, spikeFilter=sf);
% psS = psS.setConditionDescriptor(ConditionDescriptor());
% psS.minTrialsForTrialAveraging = 1;
% fprintf('spike pset tvecDataMean{1} = [%g .. %g], n=%d\n', ...
%     psS.tvecDataMean{1}(1), psS.tvecDataMean{1}(end), psS.nTimeDataMean(1));
