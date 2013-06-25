classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    properties
        name = ''; % short name, must be valid field name

        groupName = ''; % name of group to which this channel belongs 

        description = ''; % extended description 

        units = ''; % string describing units
        
        Fs = NaN; % sampling frequency in Hz 
        
        meta % anything you'd like

        storageDataClass = 'double'; % original class name, for storage purposes

        dataClass = 'double'; % classname, e.g. 'double' or 'char'

        defaultValue = NaN;

        special = false; % whether this channel is a "special" identifier channel used by TrialData
    end

    methods(Abstract)
        % return a string with this channels type
        type = getType(cdesc);

        str = describe(cdesc);

        dataFields = getExtraDataFields(cdesc);
    end

    methods
        function cd = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.parse(varargin{:});
            
            cd.name = p.Results.name;
        end
    end
end
