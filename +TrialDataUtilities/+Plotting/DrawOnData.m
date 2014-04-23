classdef DrawOnData
    methods(Static)
        function h = plotMark(axh, dMark, app, alpha, markerSize, varargin)
            % plot a set of mark points onto data
            % dMark is ? x D x ? where D is the dimensionality (2 or 3)
            import TrialDataUtilities.Plotting.patchcircle;
            import TrialDataUtilities.Plotting.patchsphere;
            
            p = inputParser;
            p.addParamValue('useTranslucentMark3d', false, @islogical);
            p.parse(varargin{:});
            
            % plot a single mark on the data
            D = size(dMark, 2);
            
            flatten = @(x) x(:);
            
            if D == 1 || D == 2
                if alpha < 1
                    h = patchcircle(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), markerSize, ...
                        'z', 0.1, 'FaceAlpha', alpha, 'FaceColor', app.Color, 'axh', axh);
                else   
                    zvals = 0.1 * ones(size(flatten(dMark(:, 1, :))));
                    h  = plot3(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), zvals, ...
                        'o', 'MarkerEdgeColor', 'none',  'MarkerFaceColor', app.Color, ...
                        'MarkerSize', markerSize, 'Parent', axh);
                end
            elseif D == 3
                if alpha < 1 && p.Results.useTranslucentMark3d
                    h = patchsphere(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), flatten(dMark(:, 3, :)), ...
                        markerSize, 'FaceColor', app.Color, 'FaceAlpha', alpha);
                else
                    h  = plot3(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), flatten(dMark(:, 3, :)), ...
                        'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', app.Color, ...
                        'Parent', axh, 'MarkerSize', markerSize);
                end
            else
                error('Invalid Dimensionality of data');
            end
        end
        
        function h = plotMarkOnRaster(axh, timesMat, app, alpha, varargin)
            % plot a set of mark points onto data displayed in a raster
            % timeMat is nOccur x nTrials matrix
            p = inputParser();
            p.addParamValue('xOffset', 0, @isscalar);
            p.addParamValue('yOffset', 0, @isscalar);
            p.addParamValue('asTick', true, @islogical);
            p.addParamValue('tickHeight', 1, @isscalar);
            p.addParamValue('tickWidth', 2, @isscalar);
            p.parse(varargin{:});

            import TrialDataUtilities.Plotting.patchcircle;
            
            nOccur = size(timesMat, 1);
            nTrials = size(timesMat, 2);
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            
            if ~p.Results.asTick
                Y = yOffset + repmat(0:-1:-(nTrials-1), nOccur, 1);
                X = xOffset + timesMat;
                X = X(:);
                Y = Y(:);
                h = patchcircle(X, Y, p.Results.tickHeight, 'z', 0.1, ...
                    'FaceAlpha', alpha, ...
                    'FaceColor', app.Color, 'axh', axh);
            else
                timesByTrial = mat2cell(timesMat, nOccur, ones(nTrials, 1))';
                h = TrialDataUtilities.Plotting.drawTickRaster(timesByTrial, ...
                    'axh', axh, 'color', app.Color, 'lineWidth', p.Results.tickWidth, ...
                    'xOffset', xOffset, 'yOffset', yOffset, 'tickHeight', p.Results.tickHeight);
            end
            
        end
        
        function h = plotInterval(axh, data, D, app, thickness, alpha)
            % plot a single interval of data
            % data is max(D,2) x T where D is dimensionality, T is number of
            % time points. D must be provided to discriminate D==1 (which
            % uses errorshade) from D==2 (which uses tubeplot)
            
            import TrialDataUtilities.Plotting.errorshade;
            
            nOccur = numel(data);
            h = nan(nOccur, 1);
            if D == 1 || D == 2
                % scale tube height in y dimension from points to data units
                %[~, yd] = TrialDataUtilities.Plotting.getPointsToAxisDataScaling(axh);
    
                for iOccur = 1:numel(data)
                    if isempty(data{iOccur}), continue; end
                    x = data{iOccur}(:, 1);
                    y = data{iOccur}(:, 2);
                  
                    if alpha < 1
                       h(iOccur) = TrialDataUtilities.Plotting.patchline(x, y, ...
                           'EdgeColor', app.Color, 'EdgeAlpha', alpha, ...
                           'LineWidth', thickness);
                    else
                        zvals = 0.09 * ones(size(x,1), 1);
                        h(iOccur) = plot3(axh, x, y, zvals, '-', ...
                            'Color', app.Color, 'LineWidth', thickness);
                    end
                end
                
            elseif D == 3
                % scale tube height in y dimension from points to data units
                %[~, yd] = TrialDataUtilities.Plotting.getPointsToAxisDataScaling(axh);
    
                for iOccur = 1:numel(data)
                    if isempty(data{iOccur}), continue; end
                    x = data{iOccur}(:, 1);
                    y = data{iOccur}(:, 2);
                    z = data{iOccur}(:, 3);
                  
                    if alpha < 1
                       h(iOccur) = TrialDataUtilities.Plotting.patchline(x, y, z, ...
                           'EdgeColor', app.Color, 'EdgeAlpha', alpha, ...
                           'LineWidth', thickness);
                    else
                        h(iOccur) = plot3(axh, x, y, z, '-', ...
                            'Color', app.Color, 'LineWidth', thickness);
                    end
                end
            end
                
        end

        function h = plotIntervalOnRaster(axh, intStart, intStop, app, alpha, varargin)
            % plot intervals on a tick raster
            % intStart/Stop are nOccur x nTrials
            p = inputParser();
            p.addParamValue('xOffset', 0, @isscalar);
            p.addParamValue('yOffset', 0, @isscalar);
            p.addParamValue('intervalHeight', 1, @isscalar);
            p.parse(varargin{:});

            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            
            nOccur = size(intStart, 1);
            nTrials = size(intStart, 2);
            yStart = yOffset + repmat(0:-1:-(nTrials-1), nOccur, 1);
            yStop = yStart - p.Results.intervalHeight;

            xStart = xOffset + intStart;
            xStop = xOffset + intStop;

            color = app.Color;
            if alpha < 1
                color = AppearanceSpec.brightenColor(color, alpha);
            end

            h = TrialDataUtilities.Plotting.patchrectangle(xStart, yStart, xStop, yStop, ...
               'FaceColor', color, 'z', -0.09, 'axh', axh);
        end
        
        function dMean = interpMarkLocation(time, data, tMean)
            % time is a time vector
            % data is a T x D set of traces
            % tMean must be vector with length nOccurrences
            % dMean will be nOccurrences x max(2,D)
            % where tMean will be inserted in row 1 if D == 1
            
            D = size(data, 2);
            
            % interpolate to find the mean location
            % dMean will be nOccurrences x D x N            
            dMean = interp1(time, data, tMean(:), 'nearest');
            
            % insert time as column 1 if D == 1
            if D == 1
                dMean = cat(2, makecol(tMean), dMean);
            end
        end
        
        function dCell = sliceIntervalLocations(time, data, tMin, tMax)
            % time is a time vector
            % data is a T x D set of traces
            % tMin and tMax must be vectors with length nOccurrences
            % dCell is nOccurrences x 1 cell of T x max(D, 2) data points 
            % that lie within each interval's occurrence
            
            T = size(data, 1);
            D = size(data, 2);
            nOccur = numel(tMin);
            
            assert(numel(time) == T);
            assert(numel(tMin) == numel(tMax));
            dCell = cell(nOccur, 1);
            
            for iOccur = 1:nOccur
                tMask = time >= tMin(iOccur) & time <= tMax(iOccur);
                
                % slice will be Tslice x D
                dCell{iOccur} = data(tMask, :);
                if D == 1
                    dCell{iOccur} = cat(2, makecol(time(tMask)), dCell{iOccur});
                end
            end
        end
    end  
end
