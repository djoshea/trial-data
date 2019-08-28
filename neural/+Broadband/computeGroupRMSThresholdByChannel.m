function rmsThreshold = computeGroupRMSThresholdByChannel(td, group, rmsMultiplier, varargin)
% extract spiking data from continuous signal by using specific threshold
% or rms value. does not sort the resulting waveforms

if nargin < 3
    rmsMultiplier = 1;
end

p = inputParser();
p.addParameter('perTrial', false, @islogical);
p.addParameter('smoothOverTrials', 0, @isscalar);
p.parse(varargin{:});
    
% fetch the data to be thresholded
contData = td.getAnalogChannelGroup(group, 'applyScaling', false);

emptyMask = cellfun(@isempty, contData);
mask = ~emptyMask & td.valid;

ssqByTrial = cellfun(@(x) sum(x-mean(x, 1, 'omitnan', 'native'), 1, 'omitnan').^2, contData(mask), 'UniformOutput', false);
ssqByTrial = double(cat(1, ssqByTrial{:}));
countByTrial = cellfun(@(x) sum(~isnan(x), 1), contData(mask), 'UniformOutput', false);
countByTrial = double(cat(1, countByTrial{:}));

if p.Results.perTrial
    rms = sqrt(ssqByTrial ./ countByTrial);
    
    if p.Results.smoothOverTrials > 1
        for iC = 1:size(rms, 2)
            rms(:, iC) = smooth(rms(:, iC), p.Results.smoothOverTrials);
        end
    end
    
    % inflate back to full nTrials size
    rms = TensorUtils.inflateMaskedTensor(rms, 1, mask);
else
    rms = sqrt(sum(ssqByTrial, 1, 'omitnan') ./ sum(countByTrial, 1, 'omitnan'));
end

% apply scaling
cd = td.getChannelDescriptor(group);
impl = cd.getImpl();
rms = impl.convertDataSingleOnAccess(1, rms);

rmsThreshold = rmsMultiplier .* rms;