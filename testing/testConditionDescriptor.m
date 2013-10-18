%clear classes;

ci = ConditionDescriptor();

ci = ci.addAttribute('a');
ci = ci.addAttribute('b');

ci = ci.binAttributeQuantiles('b', 5);

ci = ci.addAxis('a');