classdef TrialDataInterface < handle & matlab.mixin.Copyable
% TrialDataInterface is an AbstractImplementor class which acts as a wrapper
% around a set of trials. Subclasses describe how to access the data within
%
    properties(SetAccess=protected)
        data % nTrials x 1 struct array of trial data
        
        channelDescriptors % struct with ChannelDescriptor for each channel 
        
        channels % cell array of channel names
        
        groups % cell array of channel groups
    end

    properties
        meta % whatever you'd like
        
        name % string describing entire collection of trials dataset
    end

    methods(Abstract)
        % Describe the various logical groups of related channels
        % groups: scalar struct. fields are group names, values are cellstr of channel names
        groups = getChannelGroups(tdi, varargin);
        
        % Describe the channels present in the dataset 
        % channelDescriptors: scalar struct. fields are channel names, values are ChannelDescriptor 
        channelDescriptors = getChannelDescriptors(tdi, varargin);
        
        channelData = getDataForChannel(tdi, channelName, varargin);
    end
    
    % Memoized (cached) values that are obtained upon request in their get method
    properties(SetAccess=protecte, Access=public)
        nTrials 
    end

    methods % Auto-request-and-cache for above properties
        function nTrials = get.nTrials(tdi)
            if isempty(tdi.nTrials)
                tdi.nTrials = tdi.getTrialCount();
            end
            nTrials = tdi.nTrials;
        end
    end

    methods
        function tdi = TrialDataInterface(data)
            if nargin > 0
                tdi.data = data;
            end
            tdi.groups = tdi.getChannelGroups(); 
            tdi.channelDescriptors = tdi.getChannelDescriptors();
            tdi.channels = fieldnames(tdi.channelDescriptors);
        end
    end

end
