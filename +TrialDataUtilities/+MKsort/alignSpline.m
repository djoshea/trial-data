function [awaves, newTvec] = alignSpline(waves, waveTvec, varargin)
% awaves = alignSpline(waves, waveTvec, varargin)
%
% Parameters:
%   interpTroughThresh
%   alignMethod
%
% This function is not particularly fast, and could probably be optimized.
%
% waves is a matrix where each row is a waveform
%
% trough is the approximate point at which the trough occurs, in samples.
%
% alignMethod is:
% 1 for downward slope
% 2 for trough
% 3 for upward slope
% 4 first aligns by trough. Aligns by upward slope unless downward slope is
%   > 1.5 x upward slope. If so, uses downward.
% 5 for center of mass
%
% This file is part of the spike sorting software package MKsort, licensed
% under GPL version 2.
%
% Copyright (C) 2010, 2011, 2012, 2013 The Board of Trustees of The Leland
% Stanford Junior University
%
% Written by Matt Kaufman, Modified by djoshea to work within trial-data
%
% Please see mksort.m for full license and contact details.

alignMethods = {'downward', 'trough', 'upward', 'upwardDownward', 'com'};

p = inputParser();
p.addParameter('interpTroughThresh', [], @(x) isempty(x) || isscalar(x));
p.addParameter('interpFactor', 10, @isscalar);
p.addParameter('alignMethod', 'upwardDownward', @(x) ischar(x) && ismember(x, alignMethods));
p.addParameter('fillNaN', true, @islogical); 
p.parse(varargin{:});

interpTroughThresh = p.Results.interpTroughThresh;
alignMethod = p.Results.alignMethod;
fillNaN = p.Results.fillNaN;

troughPoint = find(waveTvec >= 0, 1, 'first');

waveLen = size(waves,2);
assert(waveLen == numel(waveTvec), 'Waveform time vector does not match size(waves, 2)');
nwaves = size(waves,1);
trough = troughPoint * p.Results.interpFactor;
alignpt = trough / 2;
upFactor = 1.5;

xs = waveTvec;
xis = linspace(min(waveTvec), max(waveTvec), numel(waveTvec) * p.Results.interpFactor)';
newTvec = xis;

% If interpTroughThresh is specified, clip off all points for each wave
% with values below interpTroughThresh, and spline interpolate so that we
% get a decent estimate of the trough location. If interpTroughThresh is
% >0, do the same for the post-peak.
if ~isempty(interpTroughThresh)
    for i = 1:nwaves
        if interpTroughThresh < 0
            good = waves(i, :) > interpTroughThresh;
        else
            good = waves(i, :) < interpTroughThresh;
        end
        waves(i, :) = interp1(find(good), waves(i, good), 1:waveLen, 'spline');
    end
end


% Interpolate everybody
iwaves = zeros(nwaves, length(xis));
bases = zeros(1, nwaves);
for w = 1:nwaves
    mask = ~isnan(waves(w, :));
    bases(w) = mean(waves(w, find(mask, 3)));
    if ~all(isnan(waves(w, mask)))
        iwaves(w, :) = interp1(xs(mask),waves(w,mask),xis,'spline', NaN);
    end
end

% Align. Wackiness for method 4.
if ~strcmp(alignMethod, 'upwardDownward')
    awaves = alignWaves(iwaves, trough, bases, alignpt, alignMethod, fillNaN);
    
else
    % Align on trough
    awaves = alignWaves(iwaves, trough, bases, alignpt, 'trough', fillNaN);
    cwaves = TrialDataUtilities.MKsort.cleanWaveforms(awaves);
    meanWave = nanmean(cwaves, 1);
    
    [wmin, mini] = nanmin(meanWave);
    base = mean(bases);
    % Downward slope
    down = find(meanWave(1:mini) > 0.5*(wmin-base)+base, 1, 'last');
    try
        downSlope = meanWave(down-1) - meanWave(down+1);
    catch
        downSlope = -Inf;
    end
    % Upward slope
    up = mini - 1 + find(meanWave(mini:end) > 0.5*(max(meanWave)-wmin) + wmin, 1);
    upSlope = meanWave(up+1) - meanWave(up-1);
    
    if downSlope > upFactor*upSlope
        awaves = alignWaves(iwaves, trough, bases, alignpt, 'downward', fillNaN);
    else
        awaves = alignWaves(iwaves, trough, bases, alignpt, 'upward', fillNaN);
    end
end
end

function awaves = alignWaves(waves, trough, bases, alignpt, alignMethod, fillNaN)

nwaves = size(waves,1);
waveLen = size(waves,2);
if fillNaN
    awaves = nan(nwaves, waveLen);
else
    awaves = zeros(nwaves, waveLen);
end

if strcmp(alignMethod, 'trough')
    alignpt = round(alignpt * 1.5);
elseif strcmp(alignMethod, 'upward')
    alignpt = alignpt * 3;
elseif strcmp(alignMethod, 'com')
    alignpt = alignpt * 2;
end

for w = 1:nwaves
    [wmin, mini] = min(waves(w, trough-50:trough+50));
    mini = mini + trough - 50 - 1;
    
    if strcmp(alignMethod, 'downward')  % Down
        pre = find(waves(w, 1:mini) > 0.5*(wmin-bases(w))+bases(w), 1, 'last');
    elseif strcmp(alignMethod, 'trough')  % trough
        pre = mini;
    elseif strcmp(alignMethod, 'upward')  % Up
        pre = mini - 1 + find(waves(w, mini:end) > 0.5*(nanmax(waves(w,:))-wmin) + wmin, 1);
    elseif strcmp(alignMethod, 'com')
        waves(w, isnan(waves(w, :))) = 0;
        threshold = 0.3 * wmin;  % COM
        pre = 1 + find(waves(w, 1:trough) > threshold, 1, 'last');
        post = trough - 2 + find(waves(w, trough:end) > threshold, 1);
        if isempty(post) || isempty(pre) || pre > post
            pre = mini;
        else
            cums = cumsum(waves(w, pre:post) .* (pre:post));
            pre = pre - 1 + find(cums > cums(end)/2, 1, 'last');
        end
    end
    
    if pre < alignpt
%         awaves(w, 1:(alignpt-pre)) = 0;
        awaves(w, alignpt-pre+1:end) = waves(w, 1:end-(alignpt-pre));
    else
        awaves(w, 1:(end-pre+alignpt)) = waves(w, (1+pre-alignpt):end);
    end
end
end
