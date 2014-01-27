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
    end
       
       
end
