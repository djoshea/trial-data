if ~exist('td', 'var')
    key = 'testTrialDataOlafWithHandKinematicsAndUnits';
    td = cacheLoad(key);
end

tdca = TrialDataConditionAlign(td);
tdca = tdca.groupBy('target', 'delay', 'isStim');
tdca = tdca.align('GoCue:MoveEnd').round(1);
tdca = tdca.mark('Move', 'appear', AppearanceSpec('Color', 'r', 'Marker', 'o', 'MarkerSize', 5));
tdca = tdca.colorByAttributes('delay');

%%

tdca.plotAnalogGroupMeans('handX', 'minTrials', 5, 'alpha', 0.5);
