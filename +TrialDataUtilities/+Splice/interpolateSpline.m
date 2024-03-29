function dataSpliced = interpolateSpline(dataCat, joinIdx, ...
    splineFitWindow, splineIgnoreWindow, varargin)
    
    p = inputParser;
    p.addParameter('showPlot', false, @islogical);
    p.parse(varargin{:});
    
    sz = size(dataCat);
    dataSpliced = nan(sz);
    nTraj = prod(sz(3:end));

    for c = 1:nTraj
        dataPre = dataCat(:, 1:joinIdx(c)-1, c);
        dataPost = dataCat(:, joinIdx(c):end, c);
        
        nOverwritePre = splineFitWindow + splineIgnoreWindow;
        nOverwritePost = splineFitWindow + splineIgnoreWindow;
    
        % take splineFitWindow points to either side for fitting
        % excluding the last splineIgnoreWindow points to either side
        dataFitSpline = cat(2, dataPre(:, end-splineIgnoreWindow-splineFitWindow+1:end-splineIgnoreWindow), ...
            dataPost(:, splineIgnoreWindow+1:splineIgnoreWindow+splineFitWindow));

        basisMask = all(~isnan(dataFitSpline), 2);
        
        if any(basisMask)
            dataFitSpline = dataFitSpline(basisMask, :);

            % use technique from cscvn but use smoothing spline rather than
            % interpolation
            if size(dataFitSpline, 1)==1
                dt = 0;
            else
                dt = sum((diff(dataFitSpline').^2).'); 
            end
            tFitSpline = cumsum([0,dt.^(1/4)]);
            tEvalSpline = linspace(min(tFitSpline), max(tFitSpline), nOverwritePre + nOverwritePost);

            % cs = csaps(tFitSpline, dataFitSpline, 'variational');
            cs = csaps(tFitSpline, dataFitSpline(:, :));

            % evaluate on even time grid that matches original data but
            % excludes numTimepointsDrop (which can be negative, such that we
            % add additional points)
            dataEvalSpline = fnval(cs, tEvalSpline);

            dataPreSmoothed = dataPre(:, :);
            dataPreSmoothed(basisMask, end-nOverwritePre+1:end) = dataEvalSpline(:, 1:nOverwritePre);

            dataPostSmoothed = dataPost(:, :);
            dataPostSmoothed(basisMask, 1:nOverwritePost) = dataEvalSpline(:, nOverwritePre+1:end);

            dataSpliced(:, :, c) = cat(2, dataPreSmoothed, dataPostSmoothed);
        else
            % nothing to do here, just copy in the NaNs
            dataSpliced(:, :, c) = cat(2, dataPre, dataPost);
            
        end

        if p.Results.showPlot
            if c == 1
                figure(); clf;
            end
            
            % plot traces and traces ends in red
            plot3(dataPre(1, :), dataPre(2, :), dataPre(3, :), 'k.');
            hold on;
            plot3(dataPre(1, end), dataPre(2, end), dataPre(3, end), 'ro', 'MarkerFaceColor', 'r');

            plot3(dataPost(1, :), dataPost(2, :), dataPost(3, :), 'k.');
            plot3(dataPost(1, 1), dataPost(2, 1), dataPost(3, 1), 'ro', 'MarkerFaceColor', 'r');

            % plot last timepoints retained for splicing in green
            plot3(dataPre(1, end-splineIgnoreWindow), dataPre(2, end-splineIgnoreWindow), dataPre(3, end-splineIgnoreWindow), 'go', 'MarkerFaceColor', 'g');
            plot3(dataPost(1, splineIgnoreWindow+1), dataPost(2, splineIgnoreWindow+1), dataPost(3, splineIgnoreWindow+1), 'go', 'MarkerFaceColor', 'g');

            % plot splice results
            plot3(dataSpliced(1, :, c), dataSpliced(2, :, c), dataSpliced(3, :, c), 'r');

            set(findall(gca, 'Type', 'line'), 'Clipping', 'off');
        end
    end
end