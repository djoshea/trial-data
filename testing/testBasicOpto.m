if ~exist('td', 'var')
    key = 'tdBasicOptoExample';
    td = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();
    tdca = TrialDataConditionAlign(td);
    
    tdca = tdca.align('Stim-200:StimEnd+200 @ Stim').round(1);
    tdca = tdca.mark('Stim', 'as', 'Stim Start', 'appear', eventAppear.stim, ...
        'showOnData', false);
    tdca = tdca.mark('StimPulse', 'appear', eventAppear.stim, 'showOnAxis', false);
%     tdca = tdca.interval('Stim', 'StimEnd', 'as', 'Stim', ...
%         'showOnAxis', true, 'showOnData', false, 'appear', eventAppear.stim);
    tdca = tdca.groupBy('stimPulseFreq', 'stimDuration');
    %tdca = tdca.setAttributeValueList('stimPulseFreq', 80);
    
    done = true;
end

%%

figure(1), clf; set(gcf, 'Color', 'w');
tdca.plotRaster('unit97_1', 'colorSpikes', false);
au = AutoAxis(gca);
au.axisInset(1) = 5;
au.update();

%%

figure(2), clf; set(gcf, 'Color', 'w');
tdca.plotPSTH('unit97_1', 'markAlpha', 0.5);
