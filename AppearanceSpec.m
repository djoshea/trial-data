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
