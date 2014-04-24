if ~exist('td', 'var')
    key = 'tdBasicOptoExample';
    td = cacheLoad(key);
end

if ~exist('done', 'var') || ~done
    eventAppear = getOptoReachEventAppearanceSpec();
    tdca = TrialDataConditionAlign(td);
    
    tdca = tdca.align('Stim-200:StimEnd+200 @ Stim').round(1);
    tdca = tdca.interval('StimPulse', 'StimPulseEnd+3', 'appear', eventAppear.stim, ...
        'showOnData', true, 'showOnAxis', false);
    tdca = tdca.groupBy('stimDuration', 'stimPulseFreq');
    done = true;
end

%%

figure(1), clf; set(gcf, 'Color', 'w');
tdca.plotRaster('unit97_1', 'colorSpikes', false, 'intervalAlpha', 0.8);
au = AutoAxis(gca);
au.axisInset(1) = 6;
au.update();

%%

% figure(2), clf; set(gcf, 'Color', 'w');
% tdca.plotPSTH('unit97_1', 'markAlpha', 0.5);
