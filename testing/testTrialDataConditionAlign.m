if ~exist('td', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialDataWithUnits');
end

debug('Constructing TrialDataConditionAlign\n');
tdca = TrialDataConditionAlign(td);

%tdca = tdca.groupBy({'targetX', 'targetY'});
tdca = tdca.align('@TrialStart');

tdca = tdca.selectTrials(tdca.getParam('duration') < 8000);
%tdca = tdca.round(10);

tdca = tdca.start('GoCueCommand', -200).stop('GoCueCommand', 1800).zero('GoCueCommand');

testEvent = arrayfun(@(i) round(20*rand(randi(5), 1)), 1:tdca.nTrials, ...
    'UniformOutput', false);
tdca = tdca.addEvent('TestEvent', testEvent); 

tdca = tdca.mark('GoCueCommand', +130);
tdca = tdca.mark('GoCueCommand', +170);
tdca = tdca.mark('TestEvent');
tdca = tdca.mark('TargetHeld');
tdca = tdca.mark('Success');

%tdca = tdca.interval('GoCueCommand', 'MoveOnsetOnline');

tdca = tdca.interval('GoCueCommand', 'GoCueCommand', ...
    'offsetStart', -100, 'offsetStop', 150);

tdca = tdca.interval('GoCueCommand', 'GoCueCommand', ...
    'offsetStart', -70, 'offsetStop', 75);

tdca = tdca.groupBy('targetDirection', 'delay');