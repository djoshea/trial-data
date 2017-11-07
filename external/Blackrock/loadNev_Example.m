function [Q, info] = loadNev_Example(fname, varargin)
% load nev for Eric's arm dynamics task

import(getPackageImportString);
info.version = 20171107;

% special calling mode to get current version
if nargin == 0
	Q = [];
	return;
end

% specify names of Cerebus recorded analog channels
analogInfo.analog.photobox = 129;
analogInfo.lfp = 1:128;

% trial segmentation via serial text
trialSeg.mode = 'StartEndPatterns';
trialSeg.startPattern = double('TrialStart;');
trialSeg.endPattern = double('');

trialSeg.maxTrialLength = 100e3; % ms

% no event codes directly in serial
eventInfo.lowOrderByte = false;
eventInfo.skipPulsesStart = numel(trialSeg.startPattern); 
eventInfo.skipPulsesEnd = 0;

% don't grab extra data outside of trial boundaries (since trials are
% contiguous)
spikeWindowPre = 0;
spikeWindowPost = 0;

fname = GetFullPath(fname);

Q = LoadNev.loadNev_Generic(fname, 'trialSegmentationInfo', trialSeg, ... 
	'analogInfo', analogInfo, 'eventInfo', eventInfo, ...
	'trialIdFn', [], 'type', 'ChannelReach', 'version', info.version, ...
	'spikeWindowPre', spikeWindowPre, 'spikeWindowPost', spikeWindowPost, varargin{:});

if ~isempty(Q)
    fprintf('Parsing nev serial data\n');
    % regexp to parse the contents of the serial data after 'TrialStart;'
    % above
    pat = 'TrialId=(?<trialId>\d+)[\.;]*Timestamp=(?<timestamp>\d+)[\.;]*ConditionId=(?<conditionId>\d+)';
    serial = arrayfun(@(e) char(e.codes'), [Q.evc]', 'UniformOutput', false);
    tokens = regexp(serial, pat, 'names');

    prog = ProgressBar(numel(Q), 'Parsing nev serial data');
    for iq = 1:length(Q)
        prog.update(iq);
        % copy length field
        Q(iq).cerebusLength = Q(iq).CerebusInfo.length;
        
        token = tokens{iq};
        
        if isempty(token)
           warning('Could not parse serial data for trial %d', iq);
           Q(iq).trialId = NaN;
           Q(iq).cerebusTimestamp = NaN;
           Q(iq).conditionId = NaN;
        else
           Q(iq).trialId = str2double(token(1).trialId);
           Q(iq).cerebusTimestamp = str2double(token(1).timestamp);
           Q(iq).conditionId = str2double(token(1).conditionId);
        end
    end
else
    warning('Empty Q struct returned for %s', fname);
end
prog.finish();

