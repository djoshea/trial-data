if ~exist('td', 'var')
    key = 'testTrialDataOlafWithHandKinematicsAndUnits';
    td = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();

    tdca = TrialDataConditionAlign(td);
    tdca = tdca.groupBy('target', 'delayNominal', 'isStim');
    tdca = tdca.setAttributeValueList('delayNominal', 300);
    tdca = tdca.align('TargetOnset-100:MoveEnd+100').round(1);
    tdca = tdca.mark('Move', 'appear', eventAppear.move, 'showOnAxis', false);
    tdca = tdca.mark('TargetOnset', 'appear', eventAppear.targetOnset);
    tdca = tdca.mark('MoveEnd', 'appear', eventAppear.moveEnd, 'showOnAxis', false);
    tdca = tdca.mark('GoCue', 'appear', eventAppear.goCue, 'showOnAxis', false);
    tdca = tdca.interval('Stim', 'StimEnd', 'appear', eventAppear.stim, 'showOnAxis', false);
    tdca = tdca.interval('Move', 'Move+100', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);
    tdca = tdca.colorByAttributes('target');
    
    done = true;
end

%%

% figure(1), clf, set(gcf, 'Color', 'w');
% tdca.plotAnalogGroupMeans('handX', 'minTrials', 5, 'alpha', 0.5, ...
%     'timeAxisStyle', 'marker', 'markAlpha', 0.8);
% 
% return;
% 
% %%
% 
% saveFigureSvg('testReachMeans.pdf');

%%

figure(1), clf, set(gcf, 'Color', 'w');
%tdca.plotAnalogGroupedEachTrial('handX');
tdca.plotAnalogGroupedEachTrial('handX', 'alpha', 0.5, 'timeAxisStyle', 'marker', ...
    'markAlpha', 0.5);

%%
figure(2), clf, set(gcf, 'Color', 'w');
tdca.plotAnalogGroupedEachTrial2D('handX', 'handY');

