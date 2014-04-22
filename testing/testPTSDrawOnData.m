if ~exist('pset', 'var')
    key = 'psetOlaf_PMdP';
    pset = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();

    cd = ConditionDescriptor();
    cd = cd.addAttribute('target');
    %cd = cd.addAttribute('delayNominal');
    cd = cd.addAttribute('isStim');
    cd = cd.groupBy('target', 'isStim');
    cd = cd.addAttribute('discretizedStimTrialType', 'valueList', ...
        {'delay300', 'delay300optoRelTarg320'});
    
    cd = cd.fixValueListsByApplyingToTrialData(pset.dataSources{1});
    cd = cd.colorByAttributes('target');
    
    ad = AlignDescriptor('TargetOnset-100:GoCue+100 @ GoCue');
    ad = ad.truncateAfter('GoCue+100');
    ad = ad.round(1);
    ad = ad.mark('TargetOnset', 'appear', eventAppear.targetOnset);
    ad = ad.mark('GoCue', 'appear', eventAppear.goCue);
    ad = ad.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, ...
        'showOnAxis', true);
    ad = ad.interval('TargetOnset', 'GoCue', 'as', 'Plan', 'appear', ...
        AppearanceSpec('Color', [0.9 0.7 0.7]), 'showOnAxis', false);
    
    adM = AlignDescriptor('Move-100:Move+500');
    ad = ad.round(1);
    adM = adM.mark('Move', 'appear', eventAppear.move, 'showOnAxis', false);
    adM = adM.mark('MoveEnd', 'appear', eventAppear.moveEnd, 'showOnAxis', false);
    adM = adM.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, ...
        'showOnAxis', true);
    
    %adM = adM.interval('Move', 'Move+100', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);
    
    pset = pset.setConditionDescriptor(cd);
    pset = pset.setAlignDescriptorSet({ad, adM});
    
    done = true;
end

%%
if ~exist('donePCA', 'var') || ~donePCA
    proj = ProjPCA();
    [pca, proj, statsBuild, statsProject] = proj.buildFromAndProjectPopulationTrajectorySet(pset);
    donePCA = true;
end

%%

figure(1), clf, set(gca, 'Color', 'w');
pca.plotStateSpace('markAlpha', 0.5, 'alpha', 0.6, ...
    'conditionIdx', [1 5], 'alignIdx', [1 2]);

%%
pca.plotStateSpace('markAlpha', 0.5, 'alpha', 1);


% figure(1), clf, set(gcf, 'Color', 'w');
% pset.plotAnalogGroupMeans('handX', 'minTrials', 5, 'alpha', 0.5, ...
%     'timeAxisStyle', 'marker', 'markAlpha', 0.8);
% 
% return;
% 
% %%
% 
% saveFigureSvg('testReachMeans.pdf');

% %%
% 
% figure(1), clf, set(gcf, 'Color', 'w');
% %pset.plotAnalogGroupedEachTrial('handX');
% pset.plotAnalogGroupedEachTrial('handX', 'alpha', 0.5, 'timeAxisStyle', 'marker', ...
%     'markAlpha', 0.5);
% 
% %%
% figure(2), clf, set(gcf, 'Color', 'w');
% pset.plotAnalogGroupedEachTrial2D('handX', 'handY', 'alpha', 0.5, 'markAlpha', 0.5);

%%
% figure(3), clf, set(gcf, 'Color', 'w');
% pset.plotAnalogGroupedEachTrial3D('handX', 'handY', 'handZ', 'alpha', 0.5, 'markAlpha', 0.5);

