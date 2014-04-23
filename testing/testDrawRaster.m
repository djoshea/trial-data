if ~exist('td', 'var')
    key = 'testTrialDataWithUnits';
    td = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();

    tdca = TrialDataConditionAlign(td);
    %tdca = tdca.setAttributeValueList('discretizedStimTrialType', {'delay300', 'delay300optoRelMove50'});
    
    tdca = tdca.align('TargetOnset-100:GoCue+100 @ GoCue').round(1);
    tdca = tdca.mark('GoCue', 'appear', eventAppear.goCue, 'showOnAxis', true);
    tdca = tdca.mark('TargetOnset', 'appear', eventAppear.targetOnset);
    tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, 'showOnAxis', true);
    
    tdca = tdca.addAlign('Move-100:Move+500').round(1);
    tdca = tdca.mark('Move', 'appear', eventAppear.move, 'showOnAxis', true);
    tdca = tdca.mark('MoveEnd', 'appear', eventAppear.moveEnd, 'showOnAxis', true);
    tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, 'showOnAxis', true);
    %tdca = tdca.interval('Move', 'Move+100', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);
    
    tdca = tdca.setConditionAppearanceFn(@appearFn_colorByTarget);
    
    moveDuration = tdca.useAlign(2).getEventFirst('MoveEnd');
    tdca = tdca.addParam('moveDuration', moveDuration);
    tdca = tdca.sortWithinConditionsBy('moveDuration');
    
    tdca = tdca.setAttributeValueList('discretizedStimTrialType', {'delay300', 'delay300optoRelTarg320'});
    tdca = tdca.groupBy('target', 'isStim');
    
    done = true;
end

%%

figure(6), clf; set(gcf, 'Color', 'w');
tdca.plotRaster('unit97_1', 'colorSpikes', false);

