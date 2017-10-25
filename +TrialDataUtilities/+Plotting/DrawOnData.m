classdef DrawOnData
    methods(Static)
        function h = plotMark(axh, dMark, app, alpha, markerSize, varargin)
            % plot a set of mark points onto data
            % dMark is ? x D x ? where D is the dimensionality (2 or 3)
           % import TrialDataUtilities.Plotting.patchcircle;
            %import TrialDataUtilities.Plotting.patchsphere;
            
            p = inputParser;
            p.addParameter('useTranslucentMark3d', true, @islogical);
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('zOffset', 0, @isscalar);
            p.addParameter('clipping', 'on', @ischar);
            p.addParameter('frontLayer', true, @islogical);
            p.CaseSensitive = false;
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            zOffset = p.Results.zOffset;
            
            % plot a single mark on the data
            D = size(dMark, 2);
            
            flatten = @(x) x(:);
            
            if D == 1 || D == 2
                zvals = 1 * ones(size(flatten(dMark(:, 1, :))));
                h  = plot3(flatten(dMark(:, 1, :)) + xOffset, flatten(dMark(:, 2, :)) + yOffset, zvals, ...
                    'o', 'MarkerEdgeColor', 'none',  'MarkerFaceColor', app.Color, ...
                    'MarkerSize', markerSize, 'Parent', axh, 'Clipping', p.Results.clipping);
                if alpha < 1
                    TrialDataUtilities.Plotting.setMarkerOpacity(h, alpha, 0);
                end
                
            elseif D == 3
%                 if alpha < 1 && p.Results.useTranslucentMark3d
%                     h = patchcircle(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), markerSize, ...
%                         'z', 1, 'alpha', alpha, 'color', app.Color, 'axh', axh);
% %                     h = patchsphere(flatten(dMark(:, 1, :)), flatten(dMark(:, 2, :)), flatten(dMark(:, 3, :)), ...
% %                         markerSize, 'FaceColor', app.Color, 'FaceAlpha', alpha);
%                 else
%                 h  = plot3(flatten(dMark(:, 1, :)) + xOffset, flatten(dMark(:, 2, :)) + yOffset, flatten(dMark(:, 3, :)) + zOffset, ...
%                     'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', app.Color, ...
%                     'Parent', axh, 'MarkerSize', markerSize, 'Clipping', p.Results.clipping);
                
                h = scatter3(flatten(dMark(:, 1, :)) + xOffset, flatten(dMark(:, 2, :)) + yOffset, flatten(dMark(:, 3, :)) + zOffset, markerSize, ...
                    'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', app.Color, ...
                    'Parent', axh, 'MarkerFaceAlpha', alpha, 'Clipping', p.Results.clipping);
                
%                 if alpha < 1 && p.Results.useTranslucentMark3d
%                     TrialDataUtilities.Plotting.setMarkerOpacity(h, alpha, 0);
%                 end
            else
                error('Invalid Dimensionality of data');
            end
            
            if p.Results.frontLayer
                TrialDataUtilities.Plotting.setMarkerLayerFront(h);
            end
        end
        
        function h = plotMarkOnRaster(axh, timesMat, app, alpha, varargin)
            % plot a set of mark points onto data displayed in a raster
            % timeMat is nOccur x nTrials matrix
            p = inputParser();
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('asTick', true, @islogical);
            p.addParameter('tickHeight', 1, @isscalar);
            p.addParameter('tickWidth', 2, @isscalar);
            p.parse(varargin{:});

           % import TrialDataUtilities.Plotting.patchcircle;
            
            nOccur = size(timesMat, 1);
            nTrials = size(timesMat, 2);
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            
            if ~p.Results.asTick
                Y = yOffset + repmat(0:-1:-(nTrials-1), nOccur, 1);
                X = xOffset + timesMat;
                X = X(:);
                Y = Y(:);
                h = plot3(X, Y, p.Results.tickHeight, 'z', 0.1, 'o', ...
                    'color', app.Color, 'Parent', axh);
                TrialDataUtilities.Plotting.setMarkerOpacity(h, alpha, 0);
            else
                timesByTrial = mat2cell(timesMat, nOccur, ones(nTrials, 1))';
                h = TrialDataUtilities.Plotting.drawTickRaster(timesByTrial, ...
                    'axh', axh, 'color', app.Color, 'lineWidth', p.Results.tickWidth, 'alpha', alpha, ...
                    'xOffset', xOffset, 'yOffset', yOffset, 'tickHeight', p.Results.tickHeight);
            end
            
        end
        
        function h = plotInterval(axh, data, D, app, thickness, alpha, varargin)
            % plot a single interval of data
            % data is max(D,2) x T where D is dimensionality, T is number of
            % time points. D must be provided to discriminate D==1 (which
            % uses errorshade) from D==2 (which uses tubeplot)
            
            import TrialDataUtilities.Plotting.errorshade;
            
            p = inputParser();
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('zOffset', 0, @isscalar);
            p.addParameter('style', 'line', @ischar);
            p.addParameter('clipping', 'on', @ischar);
            p.CaseSensitive = false;
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            zOffset = p.Results.zOffset;
            
            nOccur = numel(data);
            h = nan(nOccur, 1);
            if D == 1 || D == 2
                % scale tube height in y dimension from points to data units
                %[~, yd] = TrialDataUtilities.Plotting.getPointsToAxisDataScaling(axh);
    
                for iOccur = 1:numel(data)
                    if isempty(data{iOccur}), continue; end
                    x = data{iOccur}(:, 1);
                    y = data{iOccur}(:, 2);
                  
                    if strcmp(p.Results.style, 'line')
                        if alpha < 1
                           h(iOccur) = TrialDataUtilities.Plotting.patchline(x + xOffset, y + yOffset, ...
                               'EdgeColor', app.Color, 'EdgeAlpha', alpha, ...
                               'LineWidth', thickness, 'Clipping', p.Results.clipping);
                        else
                            zvals = 0.09 * ones(size(x,1), 1);
                            h(iOccur) = plot3(axh, x + xOffset, y + yOffset, zvals, '-', ...
                                'Color', app.Color, 'LineWidth', thickness, 'Clipping', p.Results.clipping);
                        end
                    elseif strcmp(p.Results.style, 'stairs')
                        h(iOccur) = TrialDataUtilities.Plotting.stairs(x + xOffset, y + yOffset, ...
                               'Color', app.Color, 'LineWidth', thickness, 'Clipping', p.Results.clipping);
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
                       h(iOccur) = TrialDataUtilities.Plotting.patchline(x + xOffset, y + yOffset, z + zOffset, ...
                           'EdgeColor', app.Color, 'EdgeAlpha', alpha, ...
                           'LineWidth', thickness, 'Clipping', p.Results.clipping);
                    else
                        h(iOccur) = plot3(axh, x + xOffset, y + yOffset, z + zOffset, '-', ...
                            'Color', app.Color, 'LineWidth', thickness, 'Clipping', p.Results.clipping);
                    end
                end
            end
                
        end

        function h = plotIntervalOnRaster(axh, intStart, intStop, app, alpha, varargin)
            % plot intervals on a tick raster
            % intStart/Stop are nOccur x nTrials
            p = inputParser();
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            p.addParameter('intervalHeight', 1, @isscalar);
            p.addParameter('intervalMinWidth', NaN, @isscalar); % if specified, ensure that each interval is drawn at least this wide
            p.parse(varargin{:});

            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            
            nOccur = size(intStart, 1);
            nTrials = size(intStart, 2);
            yStart = yOffset + repmat(0:-1:-(nTrials-1), nOccur, 1);
            yStop = yStart - p.Results.intervalHeight;
            
            xStart = xOffset + intStart;
            xStop = xOffset + intStop;
            
            % enforce minimum interval width for display (to ensure
            % visibility)
            if ~isnan(p.Results.intervalMinWidth)
                minWidth = p.Results.intervalMinWidth;
                width = xStop - xStart;
                xStop(width < minWidth) = xStart(width < minWidth) + minWidth;
            end

            color = app.Color;
%             if alpha < 1
%                 color = AppearanceSpec.brightenColor(color, alpha);
%             end

            %fprintf('Interval pre patchrect: %s', get(gcf, 'Renderer'));
           h = TrialDataUtilities.Plotting.patchrectangle(xStart, yStart, xStop, yStop, ...
               'FaceColor', color, 'FaceAlpha', alpha, 'z', -0.09, 'axh', axh);
%            if ~isempty(h)
%                h.FaceAlpha = alpha;
%            end
           %fprintf(', post patchrect: %s\n', get(gcf, 'Renderer'));
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
            [time, data] = TrialDataUtilities.Data.fixNonmonotonicTimeseries(time, data);
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
                if ~isempty(dCell{iOccur}) && D == 1
                    dCell{iOccur} = cat(2, makecol(time(tMask)), dCell{iOccur});
                end
            end
        end
    end  
end
