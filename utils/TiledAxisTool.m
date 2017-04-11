classdef TiledAxisTool < handle
    % The goal of this class is to handle computations needed to draw
    % multiple tiled plots on a shared axis. Tiles are arranged top to
    % bottom and left to right in tiles. The offsets of each tile are
    % computed automatically. You must specify the scaling of each tile --
    % in both axes, tiles can be scaled in 2 possible ways:
    %   -- fixed: in the same scale as the axes. 1:1 scaling.
    %   -- normalized: scaled relative to the axes to achieve a relative
    %   size with respect to the other tiles
    %
    % Normalized scaling allows the sizes of a set of tiles to be equalized
    % (or in relative scale) with respect to each other. Fixed scaling
    % ensures that the data in the axes will be in the same scale relative
    % to each other. These can be set independently for the x and y axes
    
    properties
        axh % AutoAxis

        spacingX % cols+1 - horizontal gaps at left, between, and at right in fraction 
        isFixedX % logical cols x 1 vector - true if fixed, false if normalized
        scaleX % numeric cols x 1 vector:
               %   if fixed, the literal scaling with respect to the data axes
               %   if normalized, the relative scaling of this tile' size
               %     with respect to other axes (fraction of width)
        limitsX % the specified X limits of this tile, cols x 2
        
        spacingY % rows+1 - vertical gaps at top, between, and at bottom in cm 
        isFixedY % rows x 1
        scaleY % rows x 1
        limitsY % rows x 1
    end
    
    properties(SetAccess=protected)
        rows
        cols
        
        % these are computed in actual data units
        cLimitsX
        cScaleX
        cOffsetX
        
        cLimitsY 
        cScaleY
        cOffsetY
    end
    
    methods
        function g = TiledAxisTool(g, rows, cols) %#ok<*INUSL>
            g.initialize(rows, cols);
        end
        
        function initialize(g, rows, cols)
            % reset the grid to use Rows x Cols, preserving existing values
            % if needed
            
            % default to fixed scale on X (typically time)
            g.spacingX = GriddedAxisTool.expandOrTruncateToSize(g.spacingX, cols, 1, 0);
            g.isFixedX = GriddedAxisTool.expandOrTruncateToSize(g.isFixedX, cols, 1, true);
            g.scaleX = GriddedAxisTool.expandOrTruncateToSize(g.scaleX, cols, 1, 1);
            g.limitsX = GriddedAxisTool.expandOrTruncateToSize(g.limits, cols, 2, NaN);
            
            % default to fixed scale on Y
            g.spacingY = GriddedAxisTool.expandOrTruncateToSize(g.spacingY, cols, 1, 0);
            g.isFixedY = GriddedAxisTool.expandOrTruncateToSize(g.isFixedY, cols, 1, true);
            g.scaleY = GriddedAxisTool.expandOrTruncateToSize(g.scaleY, cols, 1, 1);
            g.limitsY = GriddedAxisTool.expandOrTruncateToSize(g.limitsY, cols, 2, NaN);
        end
            
        function update(g)
             % compute total span in X
             
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
endf