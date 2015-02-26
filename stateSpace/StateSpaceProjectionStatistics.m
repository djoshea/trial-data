classdef StateSpaceProjectionStatistics

    properties
        basisNamesSource
        basisNamesProj
        
        nBasesSource
        nBasesProj
        
        % nBases x nBases covariance matrix
        corrSource
        covSource
        latent
        explained

        covMarginalized
        covMarginalizedNames
        latentMarginalized % variance in each basis explained by each type of variance in covMarginalizedNames
    end

    methods
        function plotExplained(proj, varargin)
            p = inputParser();
            p.addParameter('threshold', 90, @issclalar);
            p.parse(varargin{:});
            threshold = p.Results.threshold;
            
            axh = newplot;
            
            cumExplained = cumsum(proj.explained) * 100;
            plot(axh, 1:length(cumExplained), cumExplained, 'x--', ...
                'Color', [0.7 0.7 0.7], 'MarkerEdgeColor', 'r', ...
                'LineWidth', 2);
            box off;
            xlabel(axh, 'Basis');
            ylabel(axh, 'Cumulative % variance explained')

            if(cumExplained(end) >= threshold)
                hold on
                xl = get(gca, 'XLim');
                plot(axh, xl, [threshold threshold], '--', 'Color', 0.8*ones(3,1)); 

                % find the first pc that explains threshold variance
                indCross = find(cumExplained > threshold, 1, 'first');
                if ~isempty(indCross)
                    title(axh, sprintf('%d bases explain %.1f%%  of variance', indCross, threshold));
                end

                yl = get(gca, 'YLim');
                yl(2) = 100;
                ylim(axh, yl);
                hold(axh, 'off');
            end
            
            au = AutoAxis.replace(axh);
        end
        
        function plotCovSource(proj, varargin)
            clf;
            pmat(proj.covSource);
            box off;
            title('Source Covariance');
            hold off;
        end
        
        function plotBasisMixtures(proj, varargin)
            p = inputParser;
            p.addParameter('basisIdx', [1:min([10 proj.nBasesProj])], @(x) isvector(x) && ...
                all(inRange(x, [1 proj.nBasesSource])));
            p.addParameter('mixtureColors', [], @ismatrix);
            p.parse(varargin{:});
            basisIdx = p.Results.basisIdx;
            if islogical(basisIdx)
                basisIdx = find(basisIdx);
            end
            
            clf;
            
            cumBasisMix = cumsum(proj.latentMarginalized(basisIdx,:), 2);
            nCov = size(proj.latentMarginalized, 2);
            %nBases = proj.nBasesProj;
            nBases = length(basisIdx);
            rowHeight = 0.8;
            xMin = 0;
            xMax = max(cumBasisMix(:));
            hPatch = nan(nCov, 1);
            
            if isempty(p.Results.mixtureColors)
                cmap = [0.5 0.5 0.5; distinguishable_colors(nCov, {'w', [0.5 0.5 0.5]})];
            else
                cmap = p.Results.mixtureColors;                
                assert(size(cmap, 1) >= nCov && size(cmap, 2) == 3, 'mixtureColors must be nCov x 3 matrix');
            end
            
            for iIdxB = 1:nBases
                iB = basisIdx(iIdxB);
                for iCov = 1:nCov
                    % SW, NW, NE, SE order for rectangle
                    if iCov == 1
                        patchX = [0 0 cumBasisMix(iB, 1) cumBasisMix(iB, 1)];
                    else 
                        patchX = [cumBasisMix(iB, iCov-1) cumBasisMix(iB, iCov-1) ...
                            cumBasisMix(iB, iCov) cumBasisMix(iB, iCov)];
                    end
                    patchY = [nBases-iB nBases-iB+rowHeight nBases-iB+rowHeight nBases-iB];
                    
                    hPatch(iCov) = patch(patchX, patchY, cmap(iCov, :), 'EdgeColor', 'none');
                    hold on
                end
                
%                 h = text(0, nBases-iB+rowHeight/2, proj.basisNamesProj{iB}, ...
%                     'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Right');
%                 extent = get(h, 'Extent');
%                 %xMin = min(extent(1), xMin);
                
                text(cumBasisMix(iB, end), nBases-iB+rowHeight/2, sprintf('  %.2f%%', proj.explained(iB)*100), ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Left'); 
               % extent = get(h, 'Extent');
                %xMax = max(extent(1)+extent(3), xMax);
            end
            
            xlim([xMin xMax*1.1]);
            ylim([0 nBases]);
            axis off;
            title('Basis Mixtures');
            
            au = AutoAxis(gca);
            au.addTitle();
            au.addTicklessLabels('y', 'tick', nBases-0.5:-1:0.5, 'tickLabel', proj.basisNamesProj);
            
            au.addColoredLabels(proj.covMarginalizedNames, cmap, ...
                'posY', AutoAxis.PositionType.Bottom, 'posX', AutoAxis.PositionType.Right);

            au.update();
            au.installCallbacks();
            
%             legend(hPatch, nCov, proj.covMarginalizedNames, 'Location', 'NorthEastOutside');
%             legend boxoff;
            hold off;
        end
    end
end
