function [Q, analogLookup] = loadNev_Generic(filename, varargin)

p = inputParser;
p.addParameter('trialSegmentationInfo', [], @(x) true);
p.addParameter('analogInfo', [], @(x) true);
p.addParameter('eventInfo', [], @(x) true);
p.addParameter('trialIdFn', [], @(x) true);
p.addParameter('type', [], @(x) true); % gets assigned as .type to each element of Q
p.addParameter('version', [], @(x) true); % gets assigned as .version to each element of Q

p.addParameter('spikeWindowPre', 0, @isscalar); % time before trial start to grab spikes from 
p.addParameter('spikeWindowPost', 0, @isscalar); % time after trial end to grab spikes from

p.addParameter('nsxExts',  {'.ns1', '.ns2', '.ns3', '.ns4', '.ns5', '.ns6'}, @iscellstr);

p.addParameter('keepExtraAnalogChannels', true, @islogical); % if false, channels not mentioned in analogInfo will be discarded
p.parse(varargin{:});

trialSegmentationInfo = p.Results.trialSegmentationInfo;
analogInfo = p.Results.analogInfo;
eventInfo = p.Results.eventInfo;
trialIdFn = p.Results.trialIdFn;
type = p.Results.type;
version = p.Results.version;
spikeWindowPre = p.Results.spikeWindowPre;
spikeWindowPost = p.Results.spikeWindowPost;
nsxExts = p.Results.nsxExts;

if isempty(trialSegmentationInfo)
    error('No trialSegmentationInfo specified');
end

[~, name, ext] = fileparts(filename);

% append .nev extension if necessary
if isempty(ext) || ~strcmp(ext, '.nev')
    filenameNev = [filename '.nev'];
else
    filenameNev = filename;
end

% load spiking data from the nev file
debug('Loading nev file %s\n', filenameNev);
data = openNEV(filenameNev, 'nosave'); % , 'nomat');

spikeData.timestamp = double(data.Data.Spikes.TimeStamp) / 30.0;
spikeData.electrode = data.Data.Spikes.Electrode;
spikeData.unit = data.Data.Spikes.Unit;
spikeData.waveform = data.Data.Spikes.Waveform;
spikeData.waveformScaleBy = 0.25;

eventData.timestamp = double(data.Data.SerialDigitalIO.TimeStamp) / 30.0;
eventData.code = data.Data.SerialDigitalIO.UnparsedData;

if p.Results.keepExtraAnalogChannels
    channelIdList = [];
else
    % figure out which analog channel ids to load to save memory
    channelIdList = LoadNev.buildAnalogChannelIdList(analogInfo);
end

% load analog data from nsx files
nsxData = LoadNev.nevExtractAnalog(filenameNev, 'channelIds', channelIdList, 'nsxExts', nsxExts);

% Do trial segmentation
trialInfo = LoadNev.getTrialSegmentation(trialSegmentationInfo, spikeData, eventData, nsxData);

if isempty(trialInfo)
    warning('Could not segment trials');
    Q = [];
    return;
end

% Use this trialInfo to build a Q with .CerebusInfo
Q = LoadNev.segmentTrials(trialInfo);

% Add the nev file name to CerebusInfo
for iQ = 1:length(Q)
    Q(iQ).CerebusInfo.nevFile = filenameNev;
end

% segment spikes and waveforms
Q = LoadNev.addSegmentedSpikes(Q, spikeData, 'spikeWindowPre', spikeWindowPre, 'spikeWindowPost', spikeWindowPost);
clear spikeData;

% add the analog channels in nicely formatted channel groups
% with time vectors and lookup tables
[Q, analogLookup] = LoadNev.addSegmentedAnalog(Q, analogInfo, nsxData);
clear analogInfo;
clear nsxData;
    
% grab events within time period
Q = LoadNev.addSegmentedEvents(Q, eventInfo, eventData); 
clear eventInfo;
clear eventData;

if ~isempty(Q)
    % Add trial ids according to callback function provided
    if ~isempty(trialIdFn) 
        for iq = 1:length(Q)
            Q(iq).trialId = trialIdFn(Q(iq));
        end
    else
        % empty by default
        [Q.trialId] = deal([]);
    end

    % Add .type and .version fields as specified
    [Q.type] = deal(type);
    [Q.version] = deal(version);

    % include short nev file name in CerebusInfo field 
    for iq = 1:length(Q)
        Q(iq).CerebusInfo.nevNameShort = [name ext];
    end
end

