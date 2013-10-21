classdef ConditionDescriptorOnDemandCache < handle & matlab.mixin.Copyable
% This is a handle class used by condition descriptor for caching
% compute-on-demand property values. Value classes can compute
% properties on the fly in a property get method, but they can't store them
% for future quick recall. This class handles that, and making copies of it
% before modification alleviates any dependency between separate value
% instances.

    properties(Transient)
        conditions
        conditionsRelevantAttributesOnly
        
        names
        appearances
        attributeValueList
        attributeValueListAsStrings
        
        axisValueLists
    end
    
end