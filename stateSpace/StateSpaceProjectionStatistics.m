classdef StateSpaceProjectionStatistics

    properties
        nBasesSource
        nBasesProj
        
        % nBases x nBases covariance matrix
        corrSource
        covSource
        latent
        explained

        covMarginalized
        covMarginalizedNames
        latentMarginalized
        
        basisMixtures % normalized version of above, nBases x nMarginalizations
        basisMixtureNames
        basisMixtureColors
    end

    methods
        function plotExplained(proj, varargin)
            p = inputParser();
            p.addParamValue('threshold', 90, @issclalar);
            p.parse(varargin{:});
            threshold = p.Results.threshold;
            
            cumExplained = cumsum(proj.explained) * 100;
            cla;
            plot(1:length(cumExplained), cumExplained, 'x--', ...
                'Color', [0.7 0.7 0.7], 'MarkerEdgeColor', 'r', ...
                'LineWidth', 2);
            box off;
            xlabel('Basis');
            ylabel('Cumulative % variance explained')

            if(cumExplained(end) >= threshold)
                hold on
                xl = get(gca, 'XLim');
                plot(xl, [threshold threshold], '--', 'Color', 0.8*ones(3,1)); 

                % find the first pc that explains threshold variance
                indCross = find(cumExplained > threshold, 1, 'first');
                if ~isempty(indCross)
                    title(sprintf('%d bases explain %.1f%%  of variance', indCross, threshold));
                end

                yl = get(gca, 'YLim');
                yl(2) = 100;
                ylim(yl);
            end
        end
        
        function plotCovSource(proj, varargin)
            clf;
            pmat(proj.covSource);
            box off;
            title('Source Covariance');
        end
        
        function plotBasisMixtures(proj, varargin)
            p = inputParser;
            p.addParamValue('basisIdx', [1:min([10 proj.nBasesProj])], @(x) isvector(x) && ...
                all(inRange(x, [1 proj.nBasesSource])));
            p.parse(varargin{:});
            basisIdx = p.Results.basisIdx;
            if islogical(basisIdx)
                basisIdx = find(basisIdx);
            end
            
            proj.assertCommonValidTime();
            cla;
%             p = panel();
%             p.pack(1,1);
%             p(1,1).select();
            
            cumBasisMix = cumsum(proj.basisMixtures(basisIdx,:), 2);
            nCov = size(proj.basisMixtures, 2);
            %nBases = proj.nBasesProj;
            nBases = length(basisIdx);
            rowHeight = 0.8;
            xMin = 0;
            xMax = 1;
            hPatch = nan(nCov, 1);
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
                    
                    hPatch(iCov) = patch(patchX, patchY, proj.basisMixtureColors(iCov, :)); %, 'EdgeColor', 'none');
                    hold on
                end
                
                h = text(-0.03, nBases-iB+rowHeight/2, proj.basisNames{iB}, ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Right'); 
                extent = get(h, 'Extent');
                xMin = min(extent(1), xMin);
                
                h = text(1.03, nBases-iB+rowHeight/2, sprintf('%.2f%%', proj.explained(iB)*100), ...
                    'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Left'); 
                extent = get(h, 'Extent');
                xMax = max(extent(1)+extent(3), xMax);
            end
            
            xlim([xMin xMax]);
            ylim([0 nBases]);
            axis off;
            title('Basis Mixtures');
            legend(hPatch, nCov, proj.basisMixtureNames, 'Location', 'NorthEastOutside');
            legend boxoff;
            
        end
    end
end
