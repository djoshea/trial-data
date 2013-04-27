classdef ChannelDescriptor < handle & matlab.mixin.Copyable
    properties
        type % ChannelType value
        
        name = '';
        description = '';
        units = '';
        
        Fs = NaN; % sampling frequency in Hz 
        
        meta % anything you'd like
    end
        
    methods
        function d = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.addOptional('type', ChannelType.Unknown, @(x) isa(x, 'ChannelType'));
            p.parse(varargin{:});
            
            d.name = p.Results.name;
            d.type = p.Results.type;
        end
    end
end
