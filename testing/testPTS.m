if ~exist('td', 'var')
    debug('Loading test TrialData\n');
    td = cacheLoad('testTrialDataWithUnits');
end

pset = PopulationTrajectorySetBuilder.fromAllUnitsInTrialData(td);

cd = ConditionDescriptor();
cd = cd.addAttribute('targetDirection');
cd = cd.addAttribute('delayNominal');
cd = cd.groupByAll();

cd = cd.fixValueListsByApplyingToTrialData(td);

adPlan = AlignDescriptor('TargetOnset-100:GoCue+100');
adMove = AlignDescriptor('Move-100:Move+500');

pset = pset.setConditionDescriptor(cd);
pset = pset.setAlignDescriptorSet({adPlan, adMove});
