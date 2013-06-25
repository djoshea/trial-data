classdef TrialData
% TrialData represents a collection of trials, whose data is accessed via
% a TrialDataInterface. 
%
% TrialData is not a handle class, meaning that methods which modify this instance
% will return the new TrialData. Any changes made to the underlying TrialDataStore
% will be done via copy-on-write, so that changes will not propagate to another
% TrialData instance.
    
    % Obtained from TrialDataInterface
    properties% (SetAccess=protected)
        initialized = false; % has initialize() been called yet?

        trialDataInterfaceClass = '';

        data = struct();  % standardized format nTrials x 1 struct with all trial data  

        datasetName = ''; % string describing entire collection of trials dataset

        datasetMeta = struct();
        
        channelDescriptors % struct with ChannelDescriptor for each channel, by name
        
        channelNames % cell array of channel names
    end

    properties(Dependent) 
        nTrials
        nChannels
    end
    
    methods % get. accessors for above properties which simply refer to tdi.?
        function nTrials = get.nTrials(td)
            nTrials = numel(td.data);
        end

        function nChannels = get.nChannels(td)
            nChannels = numel(td.channelDescriptors);
        end
    end
    
    methods
        function td = TrialData(varargin)
            if ~isempty(varargin)
                td.initialize(varargin{:});
            end

        % copy everything over from the TrialDataInterface
        function initialize(varargin)
            p = inputParser();
            p.addOptional('trialDataInterface', @(tdi) isa(tdi, 'TrialDataInterface'));
            p.parse(varargin);
            tdi = p.Results.trialDataInterface;
            
            % copy over basic details from the TrialDataInterface
            td.trialDataInterfaceClass = class(tdi);
            td.datasetName = tdi.getDatasetName();
            td.datasetMeta = tdi.getDatasetMeta();
            td.timeUnitName = tdi.getTimeUnitName();
            td.timeUnitsPerSecond = tdi.getTimeUnitsPerSecond();
            
            % request channel descriptors for both special params and regular channels
            specialParams = makecol(tdi.getSpecialParamChannelDescriptors());
            specialNames = {specialParams.name};
            regularChannels = makecol(tdi.getChannelDescriptors());
            regularNames = {regularChannels.name};

            % check for reserved channel names
            overlap = intersect(specialNames, regularNames);
            if ~isempty(overlap)
                error('getChannelDescriptors() returned the following reserved channel names: %s', strjoin(overlap, ', '));
            end

            % combine all channelDescriptors together
            td.channelDescriptors = [specialParams; regularChannels];

            nTrials = tdi.getTrialCount();
            nChannels = numel(td.channelDescriptors);

            % request all channel data at once
            td.channelNames = {td.channelDescriptors.name};
            data = tdi.getChannelData(td.channelNames); 

            % loop over channels and verify
            for iChannel = 1:nChannels
                chd = td.channelDescriptors(iChannel); 
                name = chd.name;

                % check for main field
                assert(isfield(data, name), 'getChannelData missing field %s', name);
                
                % these are the data subfields for this channel
                chExtraFields = chd.getExtraDataFields();
                for iEx = 1:numel(chExtraFields)
                    subfield = sprintf('%s_%s', name, chExtraFields{iEx});
                    assert(isfield(data, subfield), 'getChannelData missing field %s', subfield);
                end
            end

            % TODO:
            % Replace empty values with defaults?
            % Convert from storage type to data type?

            td.data = data;

            td.initialized = true;
        end
    end
    
end

