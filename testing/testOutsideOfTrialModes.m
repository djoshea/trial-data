%% testOutsideOfTrialModes
% Characterize + check what the filtered spike rate matrix returns at OUT-OF-TRIAL edges under each
% AlignDescriptor.outsideOfTrialMode (TRUNCATE / IGNORE / INVALIDATE), and verify the
% markNanOutsideTrial flag on getSpikeRateFilteredAsMatrix.
%
% Minimal synthetic td: 1 trial, TrialStart=0, TrialEnd=100 ms, an event `start` at t=0, one spike
% channel. Aligned to [start-100, start+100] so the window extends 100 ms BEFORE the trial (i.e.
% [-100, 0) is out-of-trial). Filtered with a 10 ms rectangular (causal) filter.
%
% Expected (with the default markNanOutsideTrial=true):
%   IGNORE     -> full [-100,+100] window kept; out-of-trial [-100,0) is NaN (was a fabricated 0)
%   IGNORE (markNanOutsideTrial=false) -> historical behavior: out-of-trial [-100,0) is 0, no NaN
%   TRUNCATE   -> window clipped toward the trial extent [0,100] (also loses the causal edge)
%   INVALIDATE -> trial dropped (padded window precedes the trial) -> empty

P = {'FAILED', 'PASS'};

td = TrialData.buildEmptyWithTrialStartTrialEnd(0, 100);          % 1 trial, [0,100] ms
td = TrialDataConditionAlign(td);
td = td.addSpikeChannel('unit1_1', {[10;20;30;55;90]});          % nTrials x 1 cell of spike times (ms)
td = td.addEvent('start', 0);                                    % event at t=0
td = td.start('start', -100).stop('start', +100).zero('start'); % window [-100, +100]; trial covers [0,100]

sf = RectangularCausalSpikeFilter(10);                           % 10 ms boxcar (causal)

% ---- IGNORE, default (markNanOutsideTrial=true): out-of-trial [-100,0) must be NaN ----
tdI = td.setOutsideOfTrialIgnore();
[rates, tvec] = tdI.getSpikeRateFilteredAsMatrix('unit1_1', 'spikeFilter', sf);
oo = tvec < 0; in = tvec >= 0;                                   % startTruncated = max(0,-100) = 0
fprintf('\n== IGNORE default == tvec[%g..%g] n=%d | out-of-trial %d/%d NaN | in-trial %d/%d NaN\n', ...
    min(tvec), max(tvec), numel(tvec), nnz(isnan(rates(oo))), nnz(oo), nnz(isnan(rates(in))), nnz(in));
fprintf('[%s] IGNORE default: out-of-trial [-100,0) all NaN\n', P{1 + all(isnan(rates(oo)))});
fprintf('[%s] IGNORE default: in-trial [0,100] all finite\n',   P{1 + ~any(isnan(rates(in)))});
fprintf('[%s] IGNORE default: window reaches -100\n',           P{1 + (abs(min(tvec) + 100) < 1e-9)});

% ---- IGNORE, markNanOutsideTrial=false: reproduce the historical 0-fill ----
[r0, t0] = tdI.getSpikeRateFilteredAsMatrix('unit1_1', 'spikeFilter', sf, 'markNanOutsideTrial', false);
fprintf('[%s] IGNORE false: no NaN anywhere\n',            P{1 + ~any(isnan(r0(:)))});
fprintf('[%s] IGNORE false: out-of-trial [-100,0) all 0\n', P{1 + all(r0(t0 < 0) == 0)});

% ---- TRUNCATE: window clipped toward the trial, no out-of-trial region, all finite ----
tdT = td.setOutsideOfTrialTruncate();
[rT, tT] = tdT.getSpikeRateFilteredAsMatrix('unit1_1', 'spikeFilter', sf);
fprintf('\n== TRUNCATE == tvec[%g..%g] n=%d\n', min(tT), max(tT), numel(tT));
fprintf('[%s] TRUNCATE: window clipped to >= 0\n', P{1 + (min(tT) >= 0)});
fprintf('[%s] TRUNCATE: no NaN\n',                 P{1 + ~any(isnan(rT(:)))});

% ---- INVALIDATE: window precedes the trial, so the trial is dropped ----
tdV = td.setOutsideOfTrialInvalidate();
[rV, tV] = tdV.getSpikeRateFilteredAsMatrix('unit1_1', 'spikeFilter', sf);
dropped = isempty(rV) || all(isnan(rV(:)));
if isempty(tV)
    fprintf('\n== INVALIDATE == tvec empty | rates empty=%d\n', isempty(rV));
else
    fprintf('\n== INVALIDATE == tvec[%g..%g] n=%d | rates %d/%d NaN\n', min(tV), max(tV), numel(tV), nnz(isnan(rV(:))), numel(rV));
end
fprintf('[%s] INVALIDATE: trial dropped (empty or all-NaN)\n', P{1 + dropped});

%% Follow-up (deferred): analog counterpart + pset-level regression
% Analog is already NaN-outside-trial by construction; getAnalogAsMatrix under IGNORE spans the align
% window (expandToTimeMinMax='auto', commit 58658bd). A neural (spike) IGNORE pset built from this td now
% shows NaN at the out-of-trial edge bins (matching analog) with tvecDataMean still the full nominal
% window (the per-trial-limits widen, commit 1f80518).
