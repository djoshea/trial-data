if ~exist('pset', 'var')
    if ~cacheExists('trialData_testPset')
        debug('Loading test TrialData\n');
        td = cacheLoad('testTrialDataWithUnits');

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
        
        debug('cacheSave pset to %s\n', 'trialData_testPset');
        cacheSave('trialData_testPset', pset);
    else
        debug('Loading test PopulationTrajectorySet\n');
        pset = cacheLoad('trialData_testPset');
        
    end
end

%%
proj = ProjPCA();
proj.K = 4;
[proj, stats] = proj.buildFromPopulationTrajectorySet(pset);

[psetProj, stats2] = proj.projectPopulationTrajectorySet(pset);
