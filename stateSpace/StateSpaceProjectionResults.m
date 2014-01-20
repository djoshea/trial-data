classdef StateSpaceProjectionResults

    properties
        % nBases x nBases covariance matrix
        corrSource
        covSource
        latent
        explained

        % only populated if useCommonValidTimeWindow is true
        covMarginalized
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
    end
end
