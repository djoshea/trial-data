classdef MatUdpTrialDataInterfaceV3 < TrialDataInterface

    properties(SetAccess=protected)
        R % original struct array over trials
        meta % meta data merged over all trials

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
        function td = MatUdpTrialDataInterfaceV3(R, meta)
            assert(isvector(R) && ~isempty(R) && isstruct(R), 'Trial data must be a struct vector');
            assert(isstruct(meta) && ~isempty(meta) && isvector(meta), 'Meta data must be a vector struct');
            td.R = makecol(R);
            td.meta = meta;
            td.nTrials = numel(R);

            td.tUnits = R(1).tUnits;
        end
    end
    
    % TrialDataInterface implementation
    methods
        % return a string describing the data set wrapped by this TDI
        function datasetName = getDatasetName(tdi, varargin) %#ok<INUSD>
            datasetName = '';
        end

        % return a scalar struct containing any arbitrary metadata
        function datasetMeta = getDatasetMeta(tdi, varargin) %#ok<INUSD>
            datasetMeta = [];
        end

        % return the number of trials wrapped by this interface
        function nTrials = getTrialCount(tdi, varargin)
            nTrials = numel(tdi.R);
        end

        % return the name of the time unit used by this interface
        function timeUnitName = getTimeUnitName(tdi, varargin)
            timeUnitName = tdi.tUnits;
        end

        % Describe the channels present in the dataset 
        % channelDescriptors: scalar struct. fields are channel names, values are ChannelDescriptor 
        function channelDescriptors = getChannelDescriptors(tdi, varargin)
            % remove special fields
%             maskSpecial = ismember({fieldInfo.name}, ...
%                 {'subject', 'protocol', 'protocolVersion', 'trialId', 'duration', ...
%                  'saveTag', 'tsStartWallclock', 'tsStopWallclock', ...
%                  'tUnits', 'version', 'time'});
            iChannel = 1;
            groups = tdi.meta(1).groups;
            groupNames = fieldnames(groups);
            nGroups = numel(groupNames);
            fprintf('Inferring channel data characteristics...\n');
            for iG = 1:nGroups
                group = groups.(groupNames{iG});

                if strcmp(group.type, 'event')
                    % for event groups, the meta signalNames field is unreliable unless we take
                    % the union of all signal names
                    names = arrayfun(@(m) m.groups.(groupNames{iG}).signalNames, ...
                        tdi.meta, 'UniformOutput', false, 'ErrorHandler', @(varargin) {});
                    group.signalNames = unique(cat(1, names{:}));
                end

                nSignals = numel(group.signalNames);
                for iS = 1:nSignals
                    name = group.signalNames{iS};
                    if strncmp(name, [group.name '_'], length(group.name)+1)
                        dataFieldMain = name;
                        name = name(length(group.name)+2:end);
                    else
                        dataFieldMain = sprintf('%s_%s', group.name, name);
                    end
                    
                    dataCell = {tdi.R.(dataFieldMain)};
                    
                    switch(group.type)
                        case 'analog'
                            timeField = sprintf('%s_time', group.name);
                            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, '', tdi.tUnits);
                        case 'event'
                            cd = EventChannelDescriptor.buildMultipleEvent(name, tdi.tUnits);
                        case 'param'
                            cd = ParamChannelDescriptor(name);
                        otherwise
                            error('Unknown field type %s for channel %s', group.type, name);
                    end

                    cd.groupName = group.name;

                    
                    if strcmp(group.type, 'analog')
                        timeCell = {tdi.R.(timeField)};
                        cd = cd.inferAttributesFromData(dataCell, timeCell);
                    else
                        cd = cd.inferAttributesFromData(dataCell);
                    end

                    % store original field name to ease lookup in getDataForChannel()
                    cd.meta.originalField = dataFieldMain;

                    channelDescriptors(iChannel) = cd; %#ok<AGROW>
                    iChannel = iChannel + 1;
                end
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
        %   trialIdStr : a unique string for describing this trial
        %   subject : string, subject from whom the data were collected
        %   protocol : string, protocol in which the data were collected
        %   protocolVersion: numeric version identifier for that protocol
        %   saveTag : a numeric identifier for the containing block of trials
        %   duration : time length of each trial in tUnits
        %   timeStartWallclock : wallclock datenum indicating when this trial began
        %
        function channelData = getChannelData(tdi, channelDescriptors, varargin)
            channelData = tdi.R;
            
            % rename special channels
            channelData = mvfield(channelData, 'wallclockStart', 'timeStartWallclock');
            
            channelData = assignIntoStructArray(channelData, 'TrialStart', 0);
            channelData = copyStructField(channelData, channelData, 'duration', 'TrialEnd');
            
            for iC = 1:numel(channelDescriptors)
                cd = channelDescriptors(iC);
                if cd.special
                    continue;
                end
                
                % rename field to remove group name
                if isfield(channelData, cd.name)
                    warning('Duplicate channel name %s', cd.name);
                    continue;
                end
                channelData = mvfield(channelData, cd.meta.originalField, cd.name);
            end
        end
    end
end
