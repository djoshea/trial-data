function rmsThreshold = computeRMSThreshold(tdca, chNameAR, rmsMultiplier, varargin)
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
contData = tdca.getAnalog(chNameAR);

emptyMask = cellfun(@isempty, contData);

mask = ~emptyMask & tdca.valid;

ssqByTrial = cellfun(@(x) nansum((x-nanmean(x)).^2), contData(mask));
countByTrial = cellfun(@(x) nnz(~isnan(x)), contData(mask));

if p.Results.perTrial
    rms = sqrt(ssqByTrial ./ countByTrial);
    
    if p.Results.smoothOverTrials > 1
        rms = smooth(rms, p.Results.smoothOverTrials);
    end
    
    % inflate back to full nTrials size
    rms = TensorUtils.inflateMaskedTensor(rms, 1, mask);
else
    rms = sqrt(nansum(ssqByTrial) ./ nansum(countByTrial));
end

rmsThreshold = rmsMultiplier .* rms;