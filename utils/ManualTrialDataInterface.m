classdef ManualTrialDataInterface < TrialDataInterface

    properties
        datasetName = '';
        datasetMeta = [];
        timeUnitName
        
        % nTrials x 1
        TrialStart = []; 
        TrialEnd = [];
    end
        
    methods
        function tdi = ManualTrialDataInterface(varargin)
            p = inputParser();
            p.addParameter('timeUnitName', 'ms', @ischar);
            p.addParameter('TrialStart', [], @isvector);
            p.addParameter('TrialEnd', [], @isvector);
            p.parse(varargin{:});
            
            tdi.timeUnitName = p.Results.timeUnitName;
            tdi.TrialStart = p.Results.TrialStart;
            tdi.TrialEnd = p.Results.TrialEnd;
        end
        
        % return a string describing the data set wrapped by this TDI
        function datasetName = getDatasetName(tdi, varargin)
            datasetName = tdi.datasetName;
        end

        % return a scalar struct containing any arbitrary metadata
        function datasetMeta = getDatasetMeta(tdi, varargin)
            datasetMeta = tdi.datasetMeta;
        end

        % return the number of trials wrapped by this interface
        function nTrials = getTrialCount(tdi, varargin)
            nTrials = numel(tdi.TrialStart);
        end

        function timeUnitName = getTimeUnitName(tdi, varargin)
            timeUnitName = tdi.timeUnitName;
        end
        
        function channelDescriptors = getChannelDescriptors(tdi, varargin)
            channelDescriptors = [];
        end
        
        function channelData = getChannelData(tdi, channelNames, varargin)
            channelData = struct('TrialStart', num2cell(makecol(tdi.TrialStart)), ...
                'TrialEnd', num2cell(makecol(tdi.TrialEnd)));
        end
    end
end
