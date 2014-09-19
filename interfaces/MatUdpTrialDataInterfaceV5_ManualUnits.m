classdef MatUdpTrialDataInterfaceV5_ManualUnits < TrialDataInterface
% like V3, but uses millisecond timestamps as doubles
    properties(SetAccess=protected)
        R % original struct array over trials
        meta % meta data merged over all trials
        nTrials
        timeUnits 
    end

    methods
        % Constructor : bind to MatUdp R struct and parse channel info
        function td = MatUdpTrialDataInterfaceV5_ManualUnits(R, meta)
            assert(isvector(R) && ~isempty(R) && isstruct(R), 'Trial data must be a struct vector');
            assert(isstruct(meta) && ~isempty(meta) && isvector(meta), 'Meta data must be a vector struct');
            td.R = makecol(R);
            td.meta = meta;
            td.nTrials = numel(R);

            td.timeUnits = R(1).timeUnits;
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
            timeUnitName = tdi.timeUnits;
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
            signals = tdi.meta(1).signals;
            groupNames = fieldnames(groups);
            nGroups = numel(groupNames);
            
            nChannels = sum(cellfun(@(name) numel(groups.(name).signalNames), groupNames));
            
            prog = ProgressBar(nChannels, 'Inferring channel data characteristics');
            iSignal = 0;
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
                    iSignal = iSignal + 1;
                    prog.update(iSignal);
                    name = group.signalNames{iS};

                    if ~strcmp(group.type, 'event')
                        % event fields won't have an entry in signals table
                        signalInfo = signals.(name);
                    end
                    dataFieldMain = name;
                    
                    dataCell = {tdi.R.(dataFieldMain)};
                    
                    switch(group.type)
                        case 'analog'
                            timeField = signalInfo.timeFieldName;
                            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, signalInfo.units, tdi.timeUnits);
                            timeCell = {tdi.R.(timeField)};
                            cd = cd.inferAttributesFromData(dataCell, timeCell);
                        case 'event'
                            cd = EventChannelDescriptor.buildMultipleEvent(name, tdi.timeUnits);
                            cd = cd.inferAttributesFromData(dataCell);
                        case 'param'
                            cd = ParamChannelDescriptor.buildFromValues(name, dataCell, signalInfo.units);
                        otherwise
                            error('Unknown field type %s for channel %s', group.type, name);
                    end

                    cd.groupName = group.name;

                    % store original field name to ease lookup in getDataForChannel()
                    cd.meta.originalField = dataFieldMain;

                    channelDescriptors(iChannel) = cd; %#ok<AGROW>
                    iChannel = iChannel + 1;
                end
            end
            prog.finish();
            
            % now detect units
            fieldList = fieldnames(tdi.R);
            mask = ~cellfun(@isempty, regexp(fieldList, '^unit[\d_]+$', 'start'));
            unitFields = fieldList(mask);
            
            prog = ProgressBar(numel(unitFields), 'Adding spike units and waveforms');
            for iU = 1:numel(unitFields)
                prog.update(iU);
                cd = SpikeChannelDescriptor(unitFields{iU});
                wavefield = sprintf('%s_waveforms', unitFields{iU});
                if isfield(tdi.R, wavefield)
                    cd.waveformsField = wavefield;
                    cd.waveformsTvec = (-10:21)' / 30;
                end
                
                channelDescriptors(iChannel) = cd; %#ok<AGROW>
                iChannel = iChannel + 1;
            end
            prog.finish();
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
            debug('Converting / repairing channel data...\n');
            channelData = tdi.R;            
            
            % rename special channels
            for i = 1:numel(channelData)
                channelData(i).timeStartWallclock = channelData(i).wallclockStart;
                channelData(i).TrialStart = 0;
                channelData(i).TrialEnd = channelData(i).duration;
            end
            
            debug('Fixing waveform data\n');
            % transpose waveforms data
            for iC = 1:numel(channelDescriptors)
                cd = channelDescriptors(iC);
                if isa(cd, 'SpikeChannelDescriptor')
                    if ~isempty(cd.waveformsField)
                        for iR = 1:numel(channelData)
                            channelData(iR).(cd.waveformsField) = channelData(iR).(cd.waveformsField)';
                        end
                    end
                end
            end
        end
    end
end
