classdef AppearanceSpec
    properties(Dependent)
        Color
        LineWidth
        LineStyle
        TextOffsetX
        TextOffsetY
        HorizontalAlignment
        VerticalAlignment
    end
    
    properties(Hidden)
        mColor
        mLineWidth
        mLineStyle
        mTextOffsetX
        mTextOffsetY
        mHorizontalAlignment
        mVerticalAlignment
    end
    
    methods
        function v = get.Color(app)
            if isempty(app.mColor)
                v = app.defaultColor();
            else
                v = AppearanceSpec.convertColor(app.mColor);
            end
        end

        function app = set.Color(app, v)
            app.mColor = v;
        end

        function c = defaultColor(~)
           c = [0 0 0];
        end
       
        function v = get.LineWidth(app)
            if isempty(app.mLineWidth)
                v = app.defaultLineWidth();
            else
                v = app.mLineWidth;
            end
        end

        function app = set.LineWidth(app, v)
            app.mLineWidth = v;
        end
        
        function c = defaultLineWidth(~)
            c = 1;
        end
        
        function v = get.LineStyle(app)
            if isempty(app.mLineStyle)
                v = app.defaultLineStyle();
            else
                v = app.mLineStyle;
            end
        end

        function app = set.LineStyle(app, v)
            app.mLineStyle = v;
        end
        
        function s = defaultLineStyle(app)
            s = '-';
        end
        
        function v = get.HorizontalAlignment(app)
            if isempty(app.mHorizontalAlignment)
                v = app.defaultHorizontalAlignment();
            else
                v = app.mHorizontalAlignment;
            end
        end

        function app = set.HorizontalAlignment(app, v)
            app.mHorizontalAlignment = v;
        end
        
        function c = defaultHorizontalAlignment(~)
            c = 'center';
        end
        
         function v = get.VerticalAlignment(app)
            if isempty(app.mVerticalAlignment)
                v = app.defaultVerticalAlignment();
            else
                v = app.mVerticalAlignment;
            end
        end

        function app = set.VerticalAlignment(app, v)
            app.mVerticalAlignment = v;
        end
        
        function c = defaultVerticalAlignment(~)
            c = 'middle';
        end
        
        function v = get.TextOffsetX(app)
            if isempty(app.mTextOffsetX)
                v = app.defaultTextOffsetX();
            else
                v = app.mTextOffsetX;
            end
        end

        function app = set.TextOffsetX(app, v)
            app.mTextOffsetX = v;
        end
        
        function c = defaultTextOffsetX(~)
            c = 0;
        end
        
        function v = get.TextOffsetY(app)
            if isempty(app.mTextOffsetY)
                v = app.defaultTextOffsetY();
            else
                v = app.mTextOffsetY;
            end
        end

        function app = set.TextOffsetY(app, v)
            app.mTextOffsetY = v;
        end
        
        function c = defaultTextOffsetY(~)
            c = 0;
        end
    end
        
    methods
        function obj = AppearanceSpec(varargin)
            % instantiate with property value pair args
            nArgs = floor(numel(varargin) / 2);
            for iA = 1:nArgs
                obj.(varargin{2*iA-1}) = varargin{2*iA};
            end
        end
        
        function args = getNonEmptyArgsByName(app, argNames)
            args = {};
            for i = 1:numel(argNames)
                argName = argNames{i};
                mArgName = ['m' argName];
                
                if ~isempty(app.(mArgName))
                    args = [args, {argName, app.(mArgName)}];  %#ok<AGROW>
                end
            end
        end
       
       function args = getPlotArgs(app)
           args = app.getNonEmptyArgsByName({'Color', 'LineWidth', 'LineStyle'});
       end
       
       function args = getMarkerPlotArgs(app, showEdges)
           if nargin < 2
               showEdges = false;
           end
           if showEdges
               edge = AppearanceSpec.darkenColor(app.Color, 0.5);
           else
               edge = 'none';
           end
           args = {'MarkerFaceColor', app.Color, 'MarkerEdgeColor', edge};
       end
       
    end
       
    methods(Static)
        function cvec = convertColor(c)
            if ~ischar(c)
                cvec = c;
            else
                switch c
                    case 'b'
                        cvec = [0 0 1];
                    case 'g'
                        cvec = [0 1 0];
                    case 'r'
                        cvec = [1 0 0];
                    case 'c'
                        cvec = [0 1 1];
                    case 'm'
                        cvec = [1 0 1];
                    case 'y'
                        cvec = [1 1 0];
                    case 'k'
                        cvec = [0 0 0];
                    case 'w'
                        cvec = [1 1 1];
                    otherwise
                        error('Unknown color string');
                end
            end
        end
        
        function cvec = desaturateColor(c, gamma)
            c = AppearanceSpec.convertColor(c);
            cvec = 1 - ((1-c) .* (1-gamma));
        end
        
        function cvec = brightenColor(c, alpha)
            % brighten color by shifting fractionally towards white
            c = AppearanceSpec.convertColor(c);
            cvec = 1 - ((1-c) .* (alpha));
        end
        
        function cvec = darkenColor(c, alpha)
            % brighten color by shifting fractionally towards white
            c = AppearanceSpec.convertColor(c);
            cvec = c .* alpha;
        end
    end
end
