if ~exist('td', 'var')
    key = 'testTrialDataWithUnits';
    td = cacheLoad(key);
    td = td.selectTrials(~isnan(td.getEventFirst('Move')));
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();
    tdca = TrialDataConditionAlign(td);
    
    tdca = tdca.align('GoCue-300:MoveEnd+200 @ GoCue').round(1);
    tdca = tdca.mark('GoCue', 'appear', eventAppear.goCue, 'showOnAxis', true);
    tdca = tdca.mark('TargetOnset', 'appear', eventAppear.targetOnset, 'showOnAxis', false);
    tdca = tdca.mark('Move', 'appear', eventAppear.move, 'showOnAxis', true);
    tdca = tdca.mark('MoveEnd', 'appear', eventAppear.moveEnd, 'showOnAxis', true);
    tdca = tdca.interval('TargetOnset', 'GoCue', 'as', 'Plan', 'appear', AppearanceSpec('Color', 'r'), 'showOnAxis', false);	
    
    tdca = tdca.interval('GoCue', 'Move', 'as', 'RT', 'appear', AppearanceSpec('Color', [235 223 126] / 255), 'showOnAxis', false);	
    tdca = tdca.interval('Move', 'MoveEnd', 'as', 'Move', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);	
    
    %tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', 'appear', eventAppear.stim, 'showOnAxis', true);
    %tdca = tdca.interval('Move', 'Move+100', 'appear', AppearanceSpec('Color', 'm'), 'showOnAxis', false);
    
    tdca = tdca.sortWithinConditionsBy({'-rt'});
    
    tdca = tdca.setAttributeValueList('isStim', 0);
    tdca = tdca.groupBy('target');
    
    done = true;
end

%%

figure(6), clf; set(gcf, 'Color', 'w');
tdca.plotRaster('unit97_1', 'colorSpikes', false);

