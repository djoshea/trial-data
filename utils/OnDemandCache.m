classdef OnDemandCache < handle & matlab.mixin.Copyable
% this is a handle class designed to be held in a property by a value class,
% which enables that class to compute its properties on demand and memoize their values

    properties
        data = struct();
    end
end
