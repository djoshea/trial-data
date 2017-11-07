function Q = segmentTrials(trialInfo)
% Q = segmentTrials(Q)
% Builds a Q struct from trial info that has the appropriate start and end times
% 	in the CerebusInfo field

Q = struct('CerebusInfo', cell(length(trialInfo),1));

for iq = 1:length(trialInfo)
    % get trial times
    startTime = trialInfo(iq).startTime;
    endTime = trialInfo(iq).endTime;

    % add CerebusInfo field for FOMASH
    Q(iq).CerebusInfo.startTime = startTime;
    Q(iq).CerebusInfo.endTime = endTime;
	Q(iq).CerebusInfo.length = endTime - startTime;
end

