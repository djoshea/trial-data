if ~exist('td', 'var')
    key = 'testTrialDataOlafWithHandKinematicsAndUnits';
    td = cacheLoad(key);
end

tdca = TrialDataConditionAlign(td);
tdca = tdca.groupBy('target', 'delay', 'isStim');
tdca = tdca.align('Move:MoveEnd').round(1);
tdca = tdca.colorByAttributes('target');

%%

tdca.plotAnalogGroupMeans('handX', 'minTrials', 5, 'alpha', 0.5);