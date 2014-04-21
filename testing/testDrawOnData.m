if ~exist('td', 'var')
    key = 'testTrialDataOlafWithHandKinematicsAndUnits';
    td = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();

    tdca = TrialDataConditionAlign(td);
    tdca = tdca.groupBy('target', 'delayNominal', 'isStim');
    tdca = tdca.setAttributeValueList('discretizedStimTrialType', {'delay300', 'delay300optoRelTarg320'});
    tdca = tdca.align('TargetOnset-100:GoCue+100 @ GoCue').round(1);
    tdca = tdca.truncateBefore('GoCue');
    tdca = tdca.mark('GoCue', 'appear', eventAppear.goCue, 'showOnAxis', true);
    tdca = tdca.mark('TargetOnset', 'appear', eventAppear.targetOnset);
    tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, 'showOnAxis', true);
    
    tdca = tdca.addAlign('Move-100:Move+500').round(1);
    tdca = tdca.mark('Move', 'appear', eventAppear.move, 'showOnAxis', true);
    tdca = tdca.mark('MoveEnd', 'appear', eventAppear.moveEnd, 'showOnAxis', true);
    tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, 'showOnAxis', true);
    %tdca = tdca.interval('Move', 'Move+100', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);
    tdca = tdca.colorByAttributes('target');
%     
    done = true;
end

%%
% 
% figure(1), clf, set(gcf, 'Color', 'w');
% tdca.plotAnalogGroupMeans('handX', 'minTrials', 5, 'alpha', 0.5, ...
%     'timeAxisStyle', 'marker', 'markAlpha', 0.8);
% 
% %%
% 
% saveFigureSvg('testReachMeans.pdf');

%%

figure(5), clf, set(gcf, 'Color', 'w');
tdca.plotAnalogGroupMeans('handX', 'alpha', 0.5, 'markAlpha', 0.7);

return;

figure(1), clf, set(gcf, 'Color', 'w');
%tdca.plotAnalogGroupedEachTrial('handX');
tdca.plotAnalogGroupedEachTrial('handX', 'alpha', 0.5, 'timeAxisStyle', 'marker', ...
    'markAlpha', 0.5);

% %%
% figure(2), clf, set(gcf, 'Color', 'w');
% tdca.plotAnalogGroupedEachTrial2D('handX', 'handY', 'alpha', 0.5, 'markAlpha', 0.5);
% 
% %%
% figure(3), clf, set(gcf, 'Color', 'w');
% tdca.plotAnalogGroupedEachTrial3D('handX', 'handY', 'handZ', 'alpha', 0.5, 'markAlpha', 0.5);
% 
