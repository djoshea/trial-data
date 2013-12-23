%clear classes;

if ~exist('S', 'var')
    load ~/Downloads/tempS.mat
end

ci = ConditionInfo();
ci = ci.addAttributes(fieldnames(S));
ci = ci.applyToTrialData(S);

ci = ci.binAttributeQuantiles('rt', 5);
ci = ci.binAttributeUniform('delay', 5);

ci = ci.groupByAll();
ci = ci.axisShuffle('target');
