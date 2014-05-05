classdef AppearanceSpec
    properties(Dependent)
        Color
        LineWidth
    end
    
    properties(Hidden)
        mColor
        mLineWidth
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
           args = app.getNonEmptyArgsByName({'Color', 'LineWidth'});
       end
       
       function args = getMarkerPlotArgs(app)
           args = {'MarkerFaceColor', app.Color, 'MarkerEdgeColor', 'none'};
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
    end
end
