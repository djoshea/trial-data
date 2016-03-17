function data = removeDelaysMultiAnalogTensor(data, delays, varargin)
% delays is as given by findDelaysMultiAnalogTensor(data, delays, varargin)
% data - R (trials) x T (time) x C (channels)
% delays - R (trials) x 1
%
% parameters:
%   fillMode: can be a scalar value like NaN or 0 to fill the edges of the
%     signals as they shift left or right. or can be 'hold' to hold the
%     edge values of the signals out to the end

p = inputParser();
p.addParameter('fillMode', NaN, @(x) isscalar(x) || ischar(x));
p.parse(varargin{:});

R = size(data, 1);
T = size(data, 2);
C = size(data, 3);

hold = false;
if ischar(p.Results.fillMode)
    if strcmp(p.Results.fillMode, 'hold')
        hold = true;
        fill = NaN;
    else
        error('Unknown fillMode');
    end
else
    fill = p.Results.fillMode;
end

for r = 1:R
    idxGrab = 1:T;
    idxInsert = (1:T) - delays(r);
    mask = idxInsert > 1 & idxInsert < T;
    idxInsert = idxInsert(mask);
    idxGrab = idxGrab(mask);
    slice = data(r, :, :);
    data(r, :, :) = fill;
    data(r, idxInsert, :) = slice(1, idxGrab, :);
    
    if hold
        data(r, 1:idxInsert(1)-1, :) = repmat(data(r, idxInsert(1), :), [1 numel(1:idxInsert(1)-1) 1]);
        data(r, idxInsert(end)+1:T, :) = repmat(data(r, idxInsert(end), :), [1 numel(idxInsert(end)+1:T) 1]);
    end
end

end