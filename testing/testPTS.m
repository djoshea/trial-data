if ~exist('td', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialDataWithUnits');
end

%pset = PopulationTrajectorySetBuilder.fromAllUnitsInTrialData(td);
pset = PopulationTrajectorySetBuilder.fromAnalogChannelsInTrialData(td);

cd = ConditionDescriptor();
cd = cd.addAttribute('targetDirection');
%cd = cd.addAttribute('delayNominal');
cd = cd.groupByAll();

cd = cd.fixValueListsByApplyingToTrialData(td);
pset = pset.setConditionDescriptor(cd);

adPlan = AlignDescriptor('TargetOnset-100:GoCue+100');
adMove = AlignDescriptor('Move-100:Move+500');


pset = pset.setAlignDescriptorSet({adPlan, adMove});

%% 

