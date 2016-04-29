function [hMean, hError] = stairsError(x,y, varargin)
% like stairs, but with error regions.
% stairsError(x, y, yerror, ...) - plot +/- yerror 
% stairsError(x, y, ylo, yhi, ...) - plot confidence intervals
%
% params:
%   color - width of stairs line, defaults to next color in ColorOrder
%   line width 
%   errorColor - defaults to lighter version of color
%   errorAlpha - defaults to 0.8

    p = inputParser();
    p.addRequired('err1', [], @ismatrix);
    p.addOptional('err2', [], @ismatrix);
    p.addParameter('color', [], @(x) true);
    p.addParameter('errorColor', [], @(x) true);
    p.addParameter('errorAlpha', 0.8, @isscalar);
    p.addParamValue('lineWidth', 1, @isscalar);
    p.KeepUnmatched = true;
    p.CaseSensitive = false;
    p.parse(varargin{:});

    if isempty(p.Results.err2)
        yhi = p.Results.err1;
        ylo = p.Results.err1;
    else
        ylo = p.Results.err1;
        yhi = p.Results.err2; 
    end
    
    if isempty(p.Results.color)
        colorArg = {};
    else
        colorArg = {'Color', p.Results.color};
    end

    hMean = stairs(x, y, colorArg{:}, 'LineWidth', p.Results.lineWidth);
    
    % build polygons for fill 
    TrialDataUtilities.Plotting.hideInLegend(hError);
end
