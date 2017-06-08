classdef TiledAxisTool < handle
    % The goal of this class is to handle computations needed to draw
    % multiple tiled plots on a shared axis. Tiles are arranged top to
    % bottom and left to right in tiles. The offsets of each tile are
    % computed automatically. You must specify the scaling of each tile --
    % in both x and y axes, tiles can be scaled in 2 possible ways:
    %   - fixed: in the same scale as the data scale of the axes. 1:1
    %       scaling implies veridical scaling of data.
    %   - normalized: scaled relative to the axes to achieve a relative
    %       size with respect to the other tiles
    %
    % Normalized scaling allows the data of a set of tiles to be equalized
    % (or in relative scale) with respect to each other's data. Fixed scaling
    % ensures that the data in the axes will be in the same scale relative
    % to each other. These can be set independently for the x and y axes
    
    properties
        axh % AutoAxis

        spacingX % cols-1 - horizontal gaps between in normalized units 
        isFixedX % logical cols x 1 vector - true if fixed, false if normalized
        scaleX % numeric cols x 1 vector:
               %   if fixed, the literal scaling with respect to the data
               %   axes. 2 implies that something with input width 1 is
               %   plotted as twice as wide in axis units
               %   if normalized, the relative scaling of this tile' size
               %     with respect to other axes (fraction of width)
        limitsX % the specified X limits of this tile, cols x 2
        
        spacingY % rows-1 - vertical gaps between in normalized units
        isFixedY % rows x 1
        scaleY % rows x 1
        limitsY % rows x 1
    end
    
    properties(SetAccess=protected)
        rows
        cols
        
        % these are computed in actual data units
        cNormToDataX
        cLimitsX
        cScaleX
        cOffsetX
        cAxisLimitsX
        
        cNormToDataY
        cLimitsY 
        cScaleY
        cOffsetY
        cAxisLimitsY
    end
    
    properties(Dependent)
        cFractionX % fraction of the X axis actually taken up by each column
        cFractionY % fraction of the Y axis actually taken up by each row
    end
    
    methods
        function g = TiledAxisTool(axh, rows, cols) %#ok<*INUSL>
            g.axh = axh;
            g.initialize(rows, cols);
        end
        
        function initialize(g, rows, cols)
            % reset the grid to use Rows x Cols, preserving existing values
            % if needed
            
            g.rows = rows;
            g.cols = cols;
            
            % default to fixed scale on X (typically time)
            g.spacingX = TiledAxisTool.expandOrTruncateToSize(g.spacingX, cols-1, 1, 0.1);
            g.isFixedX = TiledAxisTool.expandOrTruncateToSize(g.isFixedX, cols, 1, true);
            g.scaleX = TiledAxisTool.expandOrTruncateToSize(g.scaleX, cols, 1, 1);
            g.limitsX = TiledAxisTool.expandOrTruncateToSize(g.limitsX, cols, 2, NaN);
            
            % default to fixed scale on Y
            g.spacingY = TiledAxisTool.expandOrTruncateToSize(g.spacingY, rows-1, 1, 0.1);
            g.isFixedY = TiledAxisTool.expandOrTruncateToSize(g.isFixedY, rows, 1, true);
            g.scaleY = TiledAxisTool.expandOrTruncateToSize(g.scaleY, rows, 1, 1);
            g.limitsY = TiledAxisTool.expandOrTruncateToSize(g.limitsY, rows, 2, NaN);
        end
            
        function update(g)
            % update all internal computations to reflect changes made to
            % setable properties
            
            % find the normalized to data unit conversion
            % if there are any fixed tiles, then the one with the smallest
            % limit span will be considered 1 normalized unit
            limitSpanX = g.limitsX(:, 2) - g.limitsX(:, 1);
            if any(g.isFixedX)
                g.cNormToDataX = min(limitSpanX(g.isFixedX));
            else
                g.cNormToDataX = 1; % nothing is fixed, so this is arbitrary
            end
            
            g.cLimitsX = nan(g.cols, 2);
            g.cScaleX = nan(g.cols, 1);
            g.cOffsetX = nan(g.cols, 1);
            offset = 0;
            for c = 1:g.cols
                if g.isFixedX(c)
                    g.cScaleX(c) = g.scaleX(c);
                    width = limitSpanX(c) * g.scaleX(c);
                else
                    width = g.scaleX(c) * g.cNormToDataX;
                    g.cScaleX(c) = limitSpanX(c) / width;
                end
                g.cLimitsX(c, :) = [offset offset+width];
                
                if c < g.cols
                    offset = offset + width + g.cNormToDataX * g.spacingX(c);
                else
                     offset = offset + width;
                end
            end
            g.cAxisLimitsX = [0 offset];
            
            % do the same for y
            limitSpanY = g.limitsY(:, 2) - g.limitsY(:, 1);
            if any(g.isFixedY)
                g.cNormToDataY = min(limitSpanY(g.isFixedY));
            else
                g.cNormToDataY = 1; % nothing is fixed, so this is arbitrary
            end
           
            g.cLimitsY = nan(g.rows, 2);
            g.cScaleY = nan(g.rows, 1);
            g.cOffsetY = nan(g.rows, 1);
            offset = 0;
            for r = 1:g.rows
                if g.isFixedY(r)
                    g.cScaleY(r) = g.scaleY(r);
                    height = limitSpanY(r) * g.scaleY(r);
                else
                    height = g.scaleY(r) * g.cNormToDataY;
                    g.cScaleY(r) = limitSpanY(r) / height;
                end
                g.cLimitsY(r, :) = [offset offset+height];
                
                if r < g.rows
                    offset = offset + height + g.cNormToDataY * g.spacingY(r);
                else
                    offset = offset + height;
                end
            end
            g.cAxisLimitsY = [0 offset];
            % reverse order since we go from top to bottom
            g.cLimitsY = offset - g.cLimitsY(:, [2 1]); 
        end
    end
    
    methods % Specifying the setup
        function colXLimits(g, col, lims)
            assert(lims(2) > lims(1));
            g.limitsX(col, :) = lims;
        end
        
        function colXLimitsNice(g, col, lims)
            lims = AutoAxisUtilities.expandLimitsToNiceNumber(lims);
            g.colXLimits(col, lims);
        end
        
        function rowYLimits(g, row, lims)
            assert(lims(2) > lims(1));
            g.limitsY(row, :) = lims;
        end
        
        function rowYLimitsNice(g, row, lims)
            lims = AutoAxisUtilities.expandLimitsToNiceNumber(lims);
            g.rowYLimits(row, lims);
        end
    end
    
    methods % Plotting on axes and plotting data
        function setAxisLimits(g)
            xlim(g.axh, g.cAxisLimitsX);
            ylim(g.axh, g.cAxisLimitsY);
        end
        
        function tx = colTransformX(g, col, x)
            tx = (x - g.limitsX(col, 1) + g.cLimitsX(col, 1)) * g.cScaleX(col);
        end
        
        function ty = rowTransformY(g, row, y)
            ty = (y - g.limitsY(row, 1) + g.cLimitsY(row, 1)) * g.cScaleY(row);
        end
        
        function addYTickBridges(g, row, varargin)
%             ceil(3*g.rows*g.cFractionY(row)))
            p = inputParser();
            p.addParameter('numTicks', 2, @isscalar);
            p.addParameter('units', '', @ischar);
            p.parse(varargin{:});
            
            ax = AutoAxis(g.axh); 
            
            ticksOrig = AutoAxisUtilities.pickNiceTickValues(g.limitsY(row, :), p.Results.numTicks);
            ticksTransform = g.rowTransformY(row, ticksOrig);
            labels = sprintfc('%g', ticksOrig);
            
            if ~isempty(p.Results.units)
                labels{end} = sprintf('%s %s', labels{end}, p.Results.units);
            end
            ax.addTickBridge('y', 'tick', ticksTransform, 'tickLabel', labels, 'otherSide', true , 'alignOuterLabelsInwards', true);
            ax.axisMarginRight = ax.axisMarginBottom;
        end
        
        function addXTickBridges(g, col, varargin)
            %  max(2, ceil(4*g.cols*g.cFractionX(col)))
            p = inputParser();
            p.addParameter('numTicks', 2, @isscalar);
            p.addParameter('units', '', @ischar);
            p.parse(varargin{:});
            
            ax = AutoAxis(g.axh); 
            ticksOrig = AutoAxisUtilities.pickNiceTickValues(g.limitsX(col, :), p.Results.numTicks);
            ticksTransform = g.colTransformX(col, ticksOrig);
            labels = sprintfc('%g', ticksOrig);
            
            if ~isempty(p.Results.units)
                labels{end} = sprintf('%s %s', labels{end}, p.Results.units);
            end
            ax.addTickBridge('x', 'tick', ticksTransform, 'tickLabel', labels, 'otherSide', false, 'alignOuterLabelsInwards', true);
        end
        
        function debug(g)
            hold(g.axh, 'on');
            g.setAxisLimits();
            % draw rectangles for each tile
            for r = 1:g.rows
                for c = 1:g.cols
                    pos = [g.cLimitsX(c, 1), g.cLimitsY(r, 1), g.cLimitsX(c, 2) - g.cLimitsX(c, 1), g.cLimitsY(r, 2) - g.cLimitsY(r, 1)];
                    rectangle('Position', pos, 'FaceColor', [1 0 0 0.2], 'EdgeColor', 'none');
                end
            end
        end
    end
    
    methods
        % Dependent properties
        function v = get.cFractionX(g)
            span = g.cLimitsX(:, 2) - g.cLimitsX(:, 1);
            v =  span / sum(span);
        end
            
        function v = get.cFractionY(g)
            span = g.cLimitsY(:, 2) - g.cLimitsY(:, 1);
            v =  span / sum(span);
        end
    end
    
    methods(Static)
        function g = demo()
            rng(1);
            clf;
            R = 5;
            C = 3;
            g = TiledAxisTool(gca, R, C);
            amp = nanvec(R);
            for r = 1:R
                amp(r) = abs(randn(1) * 5* r);
                g.rowYLimitsNice(r, [-amp(r) amp(r)]);
            end
            for c = 1:C
                g.colXLimitsNice(c, [0 10*c]);
            end
            g.update();
            
            for r = 1:R
                for c = 1:C
                    x = linspace(0, 10*c, 50);
                    y = amp(r)*sin(x);
                    tx = g.colTransformX(c, x);
                    ty = g.rowTransformY(r, y);
                    plot(tx, ty, 'k-');
                    hold on;
                end
            end
            
            for r = 1:R
                g.addYTickBridges(r);
            end
            
            for c = 1:C
                g.addXTickBridges(c);
            end
            
            ax = AutoAxis(g.axh);
            ax.update();
            
            g.debug();
        end
    end
    
    methods(Static, Hidden)
        function out = expandOrTruncateToSize(in, R, C, with)
            out = repmat(with, R, C);
            Rkeep = min(size(in, 1), R);
            Ckeep = min(size(in, 2), C); 
            if ~isempty(in)
                out(1:Rkeep, 1:Ckeep) = in(1:Rkeep, 1:Ckeep);
            end
        end
    end
end