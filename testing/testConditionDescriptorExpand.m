
cdT = ConditionDescriptor;
cdT = cdT.addAttribute('success', 'valueList', true);
cdT = cdT.addAttribute('perturbDirection', 'valueList', {'None'});
cdT = cdT.addAttribute('targetDirection', 'valueList', {'Up', 'Down'});
cdT = cdT.groupBy('targetDirection');
cdT = cdT.fixAllAxisValueLists();

cdP = ConditionDescriptor;
cdP = cdP.addAttribute('success', 'valueList', true);
cdP = cdP.addAttribute('targetDirection', 'valueList', {'Up', 'Down'});
cdP = cdP.addAttribute('perturbDirection', 'valueList', {'Left', 'Right'});
cdP = cdP.groupBy('perturbDirection', 'targetDirection');
cdP = cdP.fixAllAxisValueLists();

%%
cdE = ConditionDescriptor.expandAxesToMatch(cdT, cdP);

%%
% approach - figure out where each attribute is going to live on each axis and which attributes will be retained as filters

% start with the axes attributes on first cd
% for input i = 2:end
%   take axis a attributes. if they are a subset of another axis, then we'll add this to the axis

