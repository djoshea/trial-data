function [waves, maskKeep] = cleanWaveforms(waves, varargin)
% cleanWaves = cleanWaveforms(waves)
% Parameters:
%   meanSDThresh 
%   maxSDThresh
%
% Returns waveforms that meet certain criteria for being "clean". First
% takes z-scores at each time point. Each wave must have: 
% (1) a maximum absolute z-score less than maxSDThresh (default 2), and 
% (2) a mean absolute z-score less than meanSDThresh (default 0.4).
%
% Used in the TTP tool.
%
% This file is part of the spike sorting software package MKsort, licensed
% under GPL version 2.
% 
% Copyright (C) 2010, 2011, 2012, 2013 The Board of Trustees of The Leland
% Stanford Junior University
% 
% Written by Matt Kaufman, modified by @djoshea within trial-data
% 
% Please see mksort.m for full license and contact details.

p = inputParser();
p.addParameter('meanSDThresh', 0.4, @isscalar);
p.addParameter('maxSDThresh', 2, @isscalar);
p.parse(varargin{:});

meanSDThresh = p.Results.meanSDThresh;
maxSDThresh = p.Results.maxSDThresh;

% Find how far we can go before we lose 10% of our waveforms (from
% alignment), that's our waveform length
waveLen = 0;
for t = 1:size(waves, 2)
  if (sum(waves(:, t) == 0) < 0.1 * size(waves, 1))
    waveLen = waveLen + 1;
  end
end

waves = waves(:, 1:waveLen);

% Find waves where the mean deviation from the mean wave is < meanSDThresh
% SDs, and where the farthest outlier for each wave is < maxSDThresh SDs
% from the mean wave at that point
meanWave = nanmean(waves, 1);
stdWave = nanstd(waves, 0, 1);
stdWave(stdWave == 0) = 1;

maskKeep = false(size(waves, 1), 1);
for i = 1:size(waves, 1)
  zscores = (waves(i, :) - meanWave)./stdWave;
  if (abs(nanmean(zscores)) < meanSDThresh && nanmax(zscores) < maxSDThresh && ...
      nanmin(zscores) > -maxSDThresh)
    maskKeep(i) = true;
  end
end

waves = waves(maskKeep, :);
