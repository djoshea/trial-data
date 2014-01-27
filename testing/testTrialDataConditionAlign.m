if ~exist('td', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialData');
end

debug('Constructing TrialDataConditionAlign\n');
tdca = TrialDataConditionAlign(td);

%tdca = tdca.groupBy({'targetX', 'targetY'});
tdca = tdca.align('@TrialStart');

tdca = tdca.selectTrials(tdca.getParam('success'));
%tdca = tdca.selectTrials(tdca.getParam('duration') < 8000);
%tdca = tdca.round(10);

tdca = tdca.start('GoCue', -200).stop('GoCue', 1800).zero('GoCue');

testEvent = arrayfun(@(i) round(20*rand(randi(5), 1)), 1:tdca.nTrials, ...
    'UniformOutput', false);
tdca = tdca.addEvent('TestEvent', testEvent); 

tdca = tdca.mark('GoCue', +130);
tdca = tdca.mark('GoCue', +170);
tdca = tdca.mark('TestEvent');
tdca = tdca.mark('Move');

%tdca = tdca.interval('GoCueCommand', 'MoveOnsetOnline');

tdca = tdca.interval('Stim', 'StimEnd', ...
    'offsetStart', -100, 'offsetStop', 150);

tdca = tdca.interval('StimPulseOn', 'StimPulseOff', ...
    'offsetStart', -70, 'offsetStop', 75);

tdca = tdca.groupBy('targetDirection', 'delay');