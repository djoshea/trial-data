if ~exist('tdca', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialDataWithUnits');

    debug('Constructing TrialDataConditionAlign\n');
    tdca = TrialDataConditionAlign(td);
end

tdca = tdca.groupBy('isStim', 'target', 'delayNominal');
tdca = tdca.align('TargetOnset-100:GoCue+100');
tdca = tdca.addAlign('Move-100:Move+400');
% tdca = tdca.align('TargetOnset-100:TargetOnset+200, truncateAfter GoCue');
% tdca = tdca.addAlign('GoCue-100:GoCue+100, truncateBefore TargetOnset');
% tdca = tdca.addAlign('Move-100:Move+400');
tdca = tdca.mark('GoCue', 'appear', AppearanceSpec('Color', 'g'));
tdca = tdca.mark('TargetOnset');
tdca = tdca.interval('Stim', 'StimEnd');
tdca = tdca.round(1);

tdca = tdca.colorByAxes('isStim',{'k', 'g'});
tdca = tdca.reshapeAxes('delayNominal', 'target', 'isStim');

%%
%tdca = tdca.matchSelectConditionsAlongAxis('target', {'NE', 'NW'});

tdca.plotAnalogGroupMeansByConditionAxes('handX');