classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    properties
        name = ''; % short name, must be valid field name

        groupName = ''; % name of group to which this channel belongs 

        description = ''; % extended description 

        units = ''; % string describing units
        
        meta % anything you'd like

        Fs = NaN; % sampling frequency in Hz 

        storageDataClass = 'double'; % original class name, for storage purposes

        datenum = false; % if true, treat as datenum value 

        scalar = true; % if true, this will be treated as a scalar quantity

        defaultValue = [];

        special = false; % whether this channel is a "special" identifier channel used by TrialData
    end

    properties(Dependent)
        dataClass
        numeric
        string
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

        function tf = get.numeric(cd)
            tf = ismember(cd.storageDataClass, {'double', 'single', 'logical', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'});
        end

        function get.string(cd)
            tf = strcmp(cd.storageDataClass, 'char');
        end

        % convert all numeric types to double
        function cls = get.dataClass(cd)
            if ismember(cd.storageDataClass, {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'});
                cls = 'double';
            else 
                cls = cd.storageDataClass;
            end
        end
    end
end
