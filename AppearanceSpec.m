classdef AppearanceSpec
    properties
        Color
        LineWidth
        Marker
        MarkerSize
        MarkerFaceColor
        MarkerEdgeColor       
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
                
                if ~isempty(app.(argName))
                    args = [args, {argName, app.(argName)}]; 
                end
            end
        end
       
       function args = getPlotArgs(app)
           args = app.getNonEmptyArgsByName({'Color', 'LineWidth', 'Marker', ...
               'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'});
       end
       
       function c = get.Color(app)
           if isempty(app.Color)
               c = [0 0 0];
           else
               c = AppearanceSpec.convertColor(app.Color);
           end
       end
       
       function c = getMarkerFaceColor(app)
           if isempty(app.MarkerFaceColor)
               c = AppearanceSpec.convertColor(app.Color);
           else
               c = AppearanceSpec.convertColor(app.MarkerFaceColor);
           end
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
    end
end
