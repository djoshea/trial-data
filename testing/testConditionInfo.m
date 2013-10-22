%clear classes;

if ~exist('S', 'var')
    load ~/Downloads/tempS.mat
end

ci = ConditionInfo();
ci = ci.addAttributes(fieldnames(S));
ci = ci.applyToTrialData(S);
%ci

%ci = ci.addAttribute('delay');

ci = ci.addAxis('isOptical');

ci = ci.binAttributeQuantiles('rt', 5);
ci = ci.binAttributeUniform('delay', 5);