%clear classes;

%ci = ConditionDescriptor();
ci = ConditionInfo();

ci = ci.addAttribute('a');
ci = ci.addAttribute('b');

ci = ci.binAttributeQuantiles('b', 5);

ci = ci.addAxis('a');
ci = ci.addAxis('b');

ci

ci = ci.addAxis({'a', 'b'});