if ~exist('td', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialData');
end

debug('Constructing TrialDataConditionAlign\n');
tdca = TrialDataConditionAlign(td);

%tdca = tdca.groupBy({'targetX', 'targetY'});
tdca = tdca.align('@TrialStart');
tdca = tdca.round(10);
