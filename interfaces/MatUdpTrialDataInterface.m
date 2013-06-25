classdef MatUdpTrialDataInterface < TrialDataInterface

    properties(SetAccess=protected)
        R % original struct array over trials

        nTrials
        fields % fields of R that are not directly in support of another field, e.g _time or _tag field
        fieldsAll % fieldnames(R)

        % each is struct with fields .type .groupName, .name
        fieldInfo
        paramInfo
        eventInfo
        analogInfo

        tUnits 
    end

    methods
        % Constructor : bind to MatUdp R struct and parse channel info
        function td = MatUdpTrialDataInterface(R)
            assert(isvector(R) && ~isempty(R) && isstruct(R), 'Trial data must be a struct vector');
            td.R = makecol(R);
            td.nTrials = numel(R);

            % parse field info
            td.parseChannelInfo();             

            r = R(1);
            td.tUnits = r.tUnits;
        end

    end
    
    % TrialDataInterface implementation
    methods
        % return a string describing the data set wrapped by this TDI
        function datasetName = getDatasetName(tdi, varargin)
            datasetName = '';
        end

        % return a scalar struct containing any arbitrary metadata
        function datasetMeta = getDatasetMeta(tdi, varargin)
            datasetMeta = [];
        end

        % return the number of trials wrapped by this interface
        function nTrials = getTrialCount(tdi, varargin)
            nTrials = numel(tdi.R);
        end

        % return the name of the time unit used by this interface
        function timeUnitName = getTimeUnitName(tdi, varargin)
            timeUnitName = td.tUnits;
        end

        % return the time conversion factor, i.e. number of time units in 1 second
        function N = getTimeUnitsPerSecond(tdi, varargin)
            timeUnitName = tdi.getTimeUnitName();

            switch(timeUnitName)
                case 'ms'
                    N = 1000;
                case 's'
                    N = 1;
                otherwise
                    error('Unrecognized timeUnits %s', timeUnitName);
            end
        end

        % Describe the channels present in the dataset 
        % channelDescriptors: scalar struct. fields are channel names, values are ChannelDescriptor 
        function channelDescriptors = getChannelDescriptors(tdi, varargin)
            nFields = numel(tdi.fieldInfo);
            for iF = 1:nFields
                info = tdi.fieldInfo(iF);
                switch(info.type)
                    case 'analog'
                        cd = AnalogChannelDescriptor(info.name);
                    case 'event'
                        cd = EventChannelDescriptor(info.name);
                    case 'param'
                        cd = ParamChannelDescriptor(info.name);
                    otherwise
                        error('Unknown field type %s for channel %s', info.type, info.name);
                end

                cd.groupName = info.groupName;

                % store original field name to ease lookup in getDataForChannel()
                cd.meta.originalField = info.originalField;

                channelDescriptors(iF) = cd;
            end
        end
        
        % return a nTrials x 1 struct with all the data for each specified channel 
        % in channelNames cellstr array
        %
        % the fields of this struct are determined as:
        % .channelName for the primary data associated with that channel
        % .channelName_extraData for extra data fields associated with that channel
        %   e.g. for an AnalogChannel, .channelName_time
        %   e.g. for an EventChannel, .channelName_tags
        %
        % channelNames may include any of the channels returned by getChannelDescriptors()
        %   as well as any of the following "special" channels used by TrialData:
        %
        %   trialId : a unique numeric identifier for this trial
        %   subject : string, subject from whom the data were collected
        %   protocol : string, protocol in which the data were collected
        %   protocolVersion: numeric version identifier for that protocol
        %   saveTag : a numeric identifier for the containing block of trials
        %   duration : time length of each trial in tUnits
        %   tStartWallclock : wallclock datenum indicating when this trial began
        %   tStopWallclock : wallclock datenum indicating when this trial ended
        %
        function channelData = getChannelData(tdi, channelDescriptors, varargin)
            channelData = [];
            for iC = 1:length(channelDescriptors)
                cd = channelDescriptors(iC);

                origField = cd.meta.originalField;
                newField = cd.name;
                
                if strcmp(newField, 'handX')
                   a = 1; 
                end

                % check for original field
                assert(isfield(tdi.R, origField), 'R is missing field %s', origField);

                % handle extra fields associated with this channel
                extraFields = makecol(cd.getExtraDataFields());
                if ~isempty(extraFields)
                    origExtraFields = cellfun(@(extra) [origField '_' extra], extraFields, ...
                        'UniformOutput', false);
                    newExtraFields = cellfun(@(extra) [newField '_' extra], extraFields, ...
                        'UniformOutput', false);

                    for i = 1:length(extraFields)

                        % if original time field doesn't exist, use trial .time field
                        if strcmp(extraFields{i}, 'time')
                            if ~isfield(tdi.R, origExtraFields{i})
                                % TODO may need to check this on a per-trial basis
                                origExtraFields{i} = 'time';
                            end

                        % if this extra field doesn't exist, create it and leave it blank
                        elseif ~isfield(tdi.R, origExtraFields{i})
                            tdi.R(1).(origExtraFields{i}) = [];
                        end
                    end
                else
                    origExtraFields = {};
                    newExtraFields = {};
                end
                origFields = [{origField}; origExtraFields];
                newFields = [{newField}; newExtraFields];

                overwriting = isfield(channelData, newFields);
                if any(overwriting)
                    fprintf('WARNING: overwriting field(s) %s\n', strjoin(newFields(overwriting), ', '));
                end

                % copy over and translate field names
                channelData = copyStructField(tdi.R, channelData, origFields, newFields);
            end
        end
    end

    methods

        function parseChannelInfo(td)
            fields = fieldnames(td.R);
            td.fieldsAll = fields;
            fieldInfo = cell2mat(cellfun(@(x) td.parseFieldName(x, true), fields, 'UniformOutput', false));
            
            % filter out special derivative fields
            mask = strcmp({fieldInfo.special}, '');
            td.fields = fields(mask);
            td.fieldInfo = fieldInfo(mask);

            td.paramInfo = td.fieldInfo(strcmp({td.fieldInfo.type}, 'param'));
            td.eventInfo = td.fieldInfo(strcmp({td.fieldInfo.type}, 'event'));
            td.analogInfo = td.fieldInfo(strcmp({td.fieldInfo.type}, 'analog'));
        end

        function info = parseFieldName(td, field, withType)
            % returns the type, groupName, name for a given field
            if ~isempty(strfind(field, '__'))
                separator = '__';
                validChars = '\w_';
            else
                separator = '_';
                validChars = '^_';
            end

            if isempty(strfind(field, separator))
                % no separators, it's a meta field, so pretend its a parameter
                info.type = 'param'; 
                info.groupName = 'meta';
                info.name = field;
                info.special = '';
            else
                pattern = ['(?<groupName>[' validChars ']+)' separator '(?<name>[' validChars ']+)'];
                if withType
                    pattern = ['(?<type>[' validChars ']+)' separator pattern];
                end
                pattern = [pattern '(?<special>' separator 'time)?']; 

                info = regexp(field, pattern, 'names');

                assert(length(info) == 1, 'Error parsing field name %s', field);
                if ~isempty(info.special)
                    % remove included leading _
                    info.special = info.special(length(separator)+1:end);
                end

            end

            info.originalField = field;
        end

        function field = getFieldName(td, type, name)
            % first try finding it directly
            if ismember(td.fields, name)
                field = name;
                return;
            end

            % otherwise search all fields
            matches = strcmp({td.fieldInfo.type}, type) & ...
                      strcmp({td.fieldInfo.special}, '') & ...
                      strcmp({td.fieldInfo.name}, name);

            if any(matches)
                if nnz(matches) > 1
                    error('Multiple channels with type %s, name %s found. Specify groupName also', type, name);
                end

            else
                % try parsing it into group name and name
                info = td.parseFieldName(name, false);
                matches = strcmp({td.fieldInfo.type}, type) & ...
                          strcmp({td.fieldInfo.groupName}, info.groupName) & ...
                          strcmp({td.fieldInfo.name}, info.name);

                assert(nnz(matches) < 2, 'Multiple channels with type %s, groupName %s, name %s found', type, info.groupName, info.name);
            end

            if ~any(matches)
                error('No channel with name %s found', name);
            end

            field = td.fields{matches};
        end

        function field = getTimeFieldName(td, type, name)
            field = [td.getFieldName(type, name) '_time'];
            if ~ismember(field, td.fieldsAll)
                % no time field specified, use trial default
                field = 'time';
            end
        end
        
        function v = convertToDouble(td, v)
            if iscell(v)
                v = cellfun(@double, v, 'UniformOutput', false);
            else
                v = double(v);
            end
        end
    end
end
