% should look like 3 sinusoids, blue offset by -10 ms, red by +10, black
% with green x on top

tx = (-600:800)';
x = sin(tx / 100);

newDelta = 20;

% method = 'filter';
method = 'interp';

clf;
plot(tx, x, 'k-');
hold on;

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Causal, 'resampleMethod', method);
plot(ty, y, 'rx-');

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Acausal, 'resampleMethod', method);
plot(ty, y, 'bx-');

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Centered, 'resampleMethod', method);
plot(ty, y, 'gx');

% if causal, in is 5 and out is 4, say 0:5:100 vs. 0:4:100, then we want to
% interp the centers of the input to have centers 
% then we want 


%% % should look like 3 sinusoids superimposed, with blue shifted -1 ms, red shifted + 1 ms


tx = (-600:20:800)';
x = sin(tx / 100);

newDelta = 1;

%method = 'filter';
method = 'interp';

figure(2);
clf;
plot(tx, x, 'k-');
hold on;

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx+10, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Causal, 'resampleMethod', method);
plot(ty, y, 'rx-');

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx-10, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Acausal, 'resampleMethod', method);
plot(ty, y, 'bx-');

[y, ty] = TrialDataUtilities.Data.resampleTensorInTime(x, 1, tx, 'timeDelta', newDelta, ...
    'binAlignmentMode', BinAlignmentMode.Centered, 'resampleMethod', method);
plot(ty, y, 'gx');

% if causal, in is 5 and out is 4, say 0:5:100 vs. 0:4:100, then we want to
% interp the centers of the input to have centers 
% then we want 
