function dataSpliced = interpolateLinear(dataCat, joinIdx, ...
    method, fitWindow, ignoreWindow, varargin)
    % fit window is the number of points to either side that will be used
    % to fit the interpolating line. ignoreWindow is the number of points
    % that will be discarded at the edges and ignored by the interpolation
    % process. the fit window begins outside of the ignore window, so the
    % interpolating line will span ignoreWindow points to either
    % side of the splice.
    
    p = inputParser;
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    sz = size(dataCat);
    dataSpliced = nan(sz);
    nTraj = prod(sz(3:end));

    for c = 1:nTraj
        dataPre = dataCat(:, 1:joinIdx(c)-1, c);
        dataPost = dataCat(:, joinIdx(c):end, c);
        
        nOverwritePre = fitWindow + ignoreWindow;
        nOverwritePost = fitWindow + ignoreWindow;
    
        % take fitWindow points to either side for fitting
        % excluding the last ignoreWindow points to either side
        dataOrig = cat(2, dataPre(:, end-ignoreWindow-fitWindow+1:end-ignoreWindow), ...
            dataPost(:, ignoreWindow+1:ignoreWindow+fitWindow));
    
        [dataOrig, ndimsOrig] = TensorUtils.shiftdimToFirstDim(dataOrig, 2);
        
        % build a time vector from +/- fitWindow+ignoreWindow
        % with a gap of +/- ignore Window in the middle
        % e.g. for fitWindow = 2, ignoreWindow = 1, want [1 2 5 6]
        tOrig = cat(2, 1 : fitWindow, fitWindow + 2*ignoreWindow + 1 : 2*fitWindow + 2*ignoreWindow);
        tEval = 1 : 2*fitWindow + 2*ignoreWindow;
        
        dataInterp = interp1(tOrig, dataOrig, tEval, method);
        
        % unshift permutation
        dataInterp = TensorUtils.unshiftdimToFirstDim(dataInterp, 2, ndimsOrig);
        
        % reinsert into dataSpliced
        dataPreSmoothed = dataPre(:, :);
        dataPreSmoothed(:, end-nOverwritePre+1:end) = dataInterp(:, 1:nOverwritePre);

        dataPostSmoothed = dataPost(:, :);
        dataPostSmoothed(:, 1:nOverwritePost) = dataInterp(:, nOverwritePre+1:end);

        dataSpliced(:, :, c) = cat(2, dataPreSmoothed, dataPostSmoothed);

        if p.Results.showPlot
            % plot traces and traces ends in red
            plot3(dataPre(1, :), dataPre(2, :), dataPre(3, :), 'k.');
            hold on;
            plot3(dataPre(1, end), dataPre(2, end), dataPre(3, end), 'ro', 'MarkerFaceColor', 'r');

            plot3(dataPost(1, :), dataPost(2, :), dataPost(3, :), 'k.');
            plot3(dataPost(1, 1), dataPost(2, 1), dataPost(3, 1), 'ro', 'MarkerFaceColor', 'r');

            % plot last timepoints retained for splicing in green
            plot3(dataPre(1, end), dataPre(2, end), dataPre(3, end), 'go', 'MarkerFaceColor', 'g');
            plot3(dataPost(1, 1), dataPost(2, 1), dataPost(3, 1), 'go', 'MarkerFaceColor', 'g');

            % plot splice results
            plot3(dataSpliced(1, :), dataSpliced(2, :), dataSpliced(3, :), 'r');

            set(findall(gca, 'Type', 'line'), 'Clipping', 'off');
        end
    end
end