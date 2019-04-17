function [fn, szOut] = buildAnalogAffineTransformFn(WinByOut, bOut, varargin)
    % fn = buildAnalogAffineTransformFn(WinByOut, bOut, varargin)
    p = inputParser();
    p.addParameter('outputFn', [], @(fn) isempty(fn) || isa(fn, 'function_handle'));
    p.parse(varargin{:});
    outFn = p.Results.outputFn;
    
    szW = size(WinByOut);
    szOut = szW(2:end);
    bOut = makecol(bOut);
    if isscalar(bOut)
        bOut = repmat(bOut, TensorUtils.expandScalarSize(szOut));
    end
    assert(TensorUtils.compareSizeVectors(szOut, size(bOut)), 'Sizes of W and b are not compatible');
    
    fn = @(dataCell, timeCell, varargin) affineTransformFn(WinByOut, bOut, outFn, dataCell, timeCell, varargin{:});
end


function  [dataCell, timeCell] = affineTransformFn(W, b, outFn, dataCell, timeCell, varargin)
    p = inputParser();
    p.addParameter('scalingApplied', true, @islogical);
    p.addParameter('slice', {}, @iscell);
    p.parse(varargin{:});
    
    % out = outFn(in * WinByOut + b)
    slice = p.Results.slice;
    if ~isempty(slice)
        W = W(:, slice{:});
        b = b(slice{:});
    end
    szW = size(W);
    
    for iT = 1:numel(dataCell)
        if isempty(dataCell{iT}), continue; end
        nT = size(dataCell{iT}, 1);
        sz = [nT szW(2:end)];
        if isempty(outFn)
            dataCell{iT} = reshape(dataCell{iT}(:, :) * W, sz) + shiftdim(b, 1); 
        else
            dataCell{iT} = outFn(reshape(dataCell{iT}(:, :) * W, sz) + shiftdim(b, 1)); 
        end
    end
    
end