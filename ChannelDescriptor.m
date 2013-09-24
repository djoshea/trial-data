classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    % Manually specified meta data
    properties
        name = ''; % short name, must be valid field name

        groupName = ''; % name of group to which this channel belongs 

        description = ''; % extended description 

        units = ''; % string describing units

        meta % anything you'd like

        special = false; % whether this channel is a "special" identifier channel used by TrialData

        dfd % DataFieldDescriptor instance
    end

    % Inferred from data by inferAttributesFromData below
    properties(SetAccess=protected)
        storageDataClass = 'double'; % original class name, for storage purposes
    end
    
    properties(Dependent)
        collectAsCell
        
        missingValue
        
        dataClass
    end

    methods(Abstract)
        % return a string with this channels type
        type = getType(cdesc);

        str = describe(cdesc);

        dataFields = getDataFields(cdesc);
        
        cd = inferAttributesFromData(cd, dataCell);
    end

    methods
        function cd = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.addOptional('dfd', [], @(x) isa(x, 'DataFieldDescriptor'));
            p.parse(varargin{:});

            cd.name = p.Results.name;
            cd.dfd = p.Results.dfd;
        end

        function cd = set.dfd(cd, dfd)
            assert(isempty(dfd) || isa(dfd, 'DataFieldDescriptor'), 'dfd must be a DataFieldDescriptor');
            cd.dfd = dfd;
        end
        
        function val = get.missingValue(cd)
            if isempty(cd.dfd)
                val = [];
            else
                val =  cd.dfd.getEmptyValueElement();
            end
        end
        
        function tf = get.collectAsCell(cd)
            if isempty(cd.dfd)
                tf = true;
            else
                tf = ~cd.dfd.matrix;
            end
        end
        
        function dataClass = get.dataClass(cd)
            switch cd.storageDataClass
                case {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'}
                    dataClass = 'double';
                otherwise
                    dataClass = cd.storageDataClass;
            end
        end

        function str = getAxisLabel(cd)
            if isempty(cd.units)
                str = sprintf('%s', cd.name);
            else
                str = sprintf('%s (%s)', cd.name, cd.units);
            end
        end

        % look at the data struct and check that the correct data fields
        % are present in the R struct
        function [ok, msg] = checkData(cd, R)
            fields = cd.getDataFields();
            if ~all(isfield(R, fields))
                ok = false;
                msg = sprintf('Missing fields %s', strjoin(fields(~isfield(R, fields)), ', '));
            else
                ok = true;
                msg = '';
            end
        end

        % do any replacement of missing values, etc.
        function R = repairData(cd, R)
            % replace empty values with default in all data fields
            fields = cd.getDataFields();
            
            % for now, only repair field 1
            fld = fields{1};
            missingValue = cd.missingValue;
            if ~isempty(missingValue)
                for iR = 1:numel(R)
                    if isempty(R(iR).(fld)) || isequaln(R(iR).(fld), NaN)
                        R(iR).(fld) = missingValue;
                    end
                end
            end
        end

        function data = convertData(cd, data)
            fields = cd.getDataFields();
            data = structConvertFieldValues(data, 'double', fields);
        end 
    end
end
