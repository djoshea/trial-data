function threshEst = estimateThresholdFromSpikeWaveforms(tdca, chName, varargin)
% threshEst = estimateThresholdFromSpikeWaveforms(tdca, chName)
% chName is the name of a continuous neural channel or spike channel
% options:
%   includeChannelsOnSameArrayElectrode [ true ] : if chName refers to a
%     SpikeChannelDescriptor, look at all channels on the same array,
%     electrode as this channel. this is necessarily the case for when chName
%     is a continuous neural channel
%
% Does not touch tdca.valid (does not call reset), and so only uses spikes
% on the trials currently valid.

p = inputParser();
p.addParameter('ignoreZeroUnit', false, @islogical);
p.parse(varargin{:});
            
spikeChList = tdca.listSpikeChannelsOnSameArrayElectrodeAs(chName, 'ignoreZeroUnit', p.Results.ignoreZeroUnit);
valAtThreshCell = cellvec(numel(spikeChList));

for iC = 1:numel(spikeChList)
    [waveMat, waveTvec] = tdca.getSpikeWaveformMatrix(spikeChList{iC});

    % we take the timepoint one after zero, right after the threshold is
    % crossed
    [~, idxThreshCross] = min(abs(waveTvec));
    valAtThreshCell{iC} = waveMat(:, idxThreshCross+1);
end

valAtThresh = cat(1, valAtThreshCell{:});
if nanmean(valAtThresh) < 0
    % negative threshold
    threshEst = nanmax(valAtThresh);
else
    % positive threshold
    threshEst = nanmin(valAtThresh);
end

end
