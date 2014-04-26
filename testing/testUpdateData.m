if ~exist('td', 'var')
    key = 'testTrialDataWithUnits';
    td = cacheLoad(key);
end

tdca = TrialDataConditionAlign(td);

%%

[handX, time] = tdca.getAnalog('handX');
tdca = tdca.addAnalog('handXCopy', handX, 'handX', 'mm');

% select only some trials via condiiton info
tdca = tdca.setAttributeValueList('delayNominal', 300);
[handX, time] = tdca.getAnalog('handXCopy');
handX = cellfun(@(x) 3*abs(x), handX, 'UniformOutput', false);
% modify only those trials
tdca = tdca.setAnalog('handXCopy', handX, 'updateValidOnly', true);
tdca = tdca.setParam('target', 'new', 'updateValidOnly', true);

tdca = tdca.align('@TargetOnset');
tdca = tdca.resetConditionInfo();
tdca = tdca.groupBy('target').colorByAttributes('target');
tdca.plotAnalogGroupedEachTrial('handXCopy');
