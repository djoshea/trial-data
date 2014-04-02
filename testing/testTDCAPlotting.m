if ~exist('tdca', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialDataWithUnits');

    debug('Constructing TrialDataConditionAlign\n');
    tdca = TrialDataConditionAlign(td);
end

specs = getOptoReachEventAppearanceSpec();

tdca = tdca.groupBy('isStim', 'target', 'delayNominal');
tdca = tdca.align('TargetOnset-100:GoCue+100');
tdca = tdca.mark('GoCue', 'appear', specs.goCue);
tdca = tdca.mark('TargetOnset', 'appear', specs.targetOnset);
tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', specs.stim);
tdca = tdca.addAlign('Move-100:Move+400');
tdca = tdca.mark('Move', 'appear', specs.move);
tdca = tdca.round(1);

%tdca = tdca.colorByAxes('isStim',{'k', 'g'});
tdca = tdca.colorByAttributes('target');
tdca = tdca.reshapeAxes('delayNominal', [], {{'target'}, {'isStim'}});

% tdca = tdca.align('TargetOnset-100:TargetOnset+200, truncateAfter GoCue');
% tdca = tdca.addAlign('GoCue-100:GoCue+100, truncateBefore TargetOnset');
% tdca = tdca.addAlign('Move-100:Move+400');


%%
%tdca = tdca.matchSelectConditionsAlongAxis('target', {'NE', 'NW'});

tdca.plotAnalogGroupMeansByConditionAxes('handSpeed', 'minTrials', 5);