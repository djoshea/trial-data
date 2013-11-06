classdef ConditionInfoInternalHandleStore < handle
% this class holds onto data in ConditionInfo that would otherwise need to be recomputed on the fly

            ci.conditions = [];
            ci.appearances = {};
            ci.names = {};
            ci.conditionIdx = [];
            ci.conditionSubs = [];
            ci.conditionSubsIncludingManualInvalid = [];
            ci.listByCondition = [];
end
