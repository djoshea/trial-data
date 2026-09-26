function test_degenerate_segments_through_resample_chain()
% Drive DEGENERATE per-trial segments through the analog resample / embed chain and report
% everything that breaks, instead of discovering them one dataset at a time.
%
% WHY THIS EXISTS. A window pinned to a STIM event sits a fixed distance from each trial's
% data, so within a pool either every trial can fill it or none can -- and the all-empty case
% is handled by callers. A window pinned to MOVEMENT (e.g. EstimReach's
% Estim<Fam>AlignMovePosthocCommon) sits at a distance that varies per trial with RT, so a
% handful of trials come back with zero samples, one sample, or all-NaN while the rest are
% fine. That MIXED case reached this chain for the first time in Sep 2026 and broke it in
% three separate places, each found only by running the next dataset:
%
%   1. resamplePadEdges           -- resample() rejects an empty X
%   2. resamplePadEdges           -- makecol turns a 1 x nCh segment into nCh x 1, silently
%                                    reinterpreting channels as timepoints
%   3. inferTimeDeltaFromSampleTimes -- median(diff(t)) on one sample is empty and will not
%                                    assign into a scalar slot
%
% Each is a legitimate input meaning "this trial has no usable data in this window", and the
% correct answer is all-NaN on the requested output time base, not an error.
%
% Run: test_degenerate_segments_through_resample_chain

    nCh = 10;          % EstimReach's 'hand' group size -- the count that made case 2 visible
    dtIn = 1;
    dtOut = 10;
    nGood = 40;

    fprintf('Degenerate segments through the resample/embed chain (nCh = %d)\n\n', nCh);

    cases = {
        'zero samples',        zeros(0, nCh),                  zeros(0, 1)
        'one sample',          rand(1, nCh),                   0
        'two samples',         rand(2, nCh),                   [0; dtIn]
        'all-NaN, many rows',  nan(nGood, nCh),                (0:nGood-1)' * dtIn
        'one all-NaN channel', makeOneNaNCol(nGood, nCh),      (0:nGood-1)' * dtIn
        'normal',              rand(nGood, nCh),               (0:nGood-1)' * dtIn
    };

    ty = (0:dtOut:(nGood-1)*dtIn)';
    nFail = 0;

    %% 1. resamplePadEdges directly, one segment at a time
    fprintf('--- resamplePadEdges ---\n');
    for i = 1:size(cases, 1)
        [name, x, tx] = cases{i, :};
        try
            y = TrialDataUtilities.Data.resamplePadEdges(x, tx, ty, BinAlignmentMode.Causal, 'linear');
            ok = isequal(size(y), [numel(ty), nCh]);
            fprintf('  %-22s -> %-12s %s\n', name, mat2str(size(y)), passfail(ok, ...
                sprintf('expected %s', mat2str([numel(ty), nCh]))));
            nFail = nFail + ~ok;
        catch ME
            fprintf('  %-22s -> ERROR %s: %s\n', name, ME.identifier, oneline(ME.message));
            nFail = nFail + 1;
        end
    end

    %% 2. resampleDataCellInTime -- the MIXED cell, which is the real situation
    fprintf('\n--- resampleDataCellInTime (all cases mixed in one cell, as a real pool is) ---\n');
    dataCell = cases(:, 2);
    timeCell = cases(:, 3);
    try
        [d, t] = TrialDataUtilities.Data.resampleDataCellInTime(dataCell, timeCell, ...
            'timeDelta', dtOut, 'binAlignmentMode', BinAlignmentMode.Causal);
        fprintf('  returned %d segments, %d time vectors\n', numel(d), numel(t));
        for i = 1:numel(d)
            fprintf('    %-22s -> %s\n', cases{i, 1}, mat2str(size(d{i})));
        end
    catch ME
        fprintf('  ERROR %s: %s\n', ME.identifier, oneline(ME.message));
        nFail = nFail + 1;
    end

    %% 3. inferTimeDeltaFromSampleTimes -- both branches, since only one had the guard
    fprintf('\n--- inferTimeDeltaFromSampleTimes ---\n');
    for ig = [false true]
        try
            td = TrialDataUtilities.Data.inferTimeDeltaFromSampleTimes(timeCell, dataCell, ...
                'ignoreNaNSamples', ig);
            fprintf('  ignoreNaNSamples=%d -> %g %s\n', ig, td, passfail(isscalar(td) && ~isempty(td), 'expected a scalar'));
            nFail = nFail + ~(isscalar(td) && ~isempty(td));
        catch ME
            fprintf('  ignoreNaNSamples=%d -> ERROR %s: %s\n', ig, ME.identifier, oneline(ME.message));
            nFail = nFail + 1;
        end
    end

    %% 4. embedTimeseriesInMatrix -- the tensor assembly the exporter actually calls
    fprintf('\n--- embedTimeseriesInMatrix (mixed pool) ---\n');
    try
        [mat, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell, timeCell);
        fprintf('  -> %s over tvec %g..%g (%d bins) %s\n', mat2str(size(mat)), ...
            tvec(1), tvec(end), numel(tvec), passfail(size(mat, 1) == numel(dataCell), ...
            'expected one row per segment'));
        nFail = nFail + ~(size(mat, 1) == numel(dataCell));
    catch ME
        fprintf('  ERROR %s: %s\n', ME.identifier, oneline(ME.message));
        nFail = nFail + 1;
    end

    %% 5. embedTimeseriesInMatrix where EVERY segment is degenerate.
    % Taken from a real failing call on EstimReach P20180609_C: 40 segments, 39 of one sample
    % and 1 of none. Every segment being single-sample makes the pooled origDelta NaN, so
    % there is no common time vector at all -- tvec and tMinGlobal come back empty while tMin
    % stays full size, and the index arithmetic subtracts one from the other. Distinct from
    % the mixed pool above, which always has some trial to define the time base.
    fprintf('\n--- embedTimeseriesInMatrix (EVERY segment degenerate) ---\n');
    nAll = 40;
    allDegData = [repmat({rand(1, nCh)}, nAll - 1, 1); {zeros(0, nCh)}];
    allDegTime = [repmat({0}, nAll - 1, 1); {zeros(0, 1)}];
    try
        [mat2, tvec2] = TrialDataUtilities.Data.embedTimeseriesInMatrix(allDegData, allDegTime);
        ok = size(mat2, 1) == nAll;
        fprintf('  -> %s, tvec %d bins %s\n', mat2str(size(mat2)), numel(tvec2), ...
            passfail(ok, sprintf('expected %d rows', nAll)));
        nFail = nFail + ~ok;
    catch ME
        fprintf('  ERROR %s: %s\n', ME.identifier, oneline(ME.message));
        nFail = nFail + 1;
    end

    fprintf('\n%s (%d failures)\n', ternary(nFail == 0, 'ALL PASS', 'FAILURES PRESENT'), nFail);
    assert(nFail == 0, 'test_degenerate_segments_through_resample_chain: %d failures', nFail);
end

function x = makeOneNaNCol(n, nCh)
    x = rand(n, nCh);
    x(:, 2) = NaN;
end

function s = passfail(tf, why)
    if tf, s = 'OK'; else, s = sprintf('<-- WRONG SHAPE, %s', why); end
end

function s = oneline(s)
    s = regexprep(s, '\s+', ' ');
    if strlength(s) > 110, s = extractBefore(s, 110) + "..."; end
end

function v = ternary(tf, a, b)
    if tf, v = a; else, v = b; end
end
