function [fn, szOut] = buildAnalogSliceTransformFn(sliceArgs, varargin)
    % [fn, szOut] = buildAnalogSliceTransformFn(WinByOut, bOut, varargin)
    % build analog channel group that simply slices into another channel
    % group
    
    if ~iscell(sliceArgs)
        sliceArgs = {sliceArgs};
    end
    szOut = cellfun(@(a) numel(a), sliceArgs);
    
    fn = @(dataCell, timeCell, varargin) sliceTransformFn(sliceArgs, dataCell, timeCell, varargin{:});
end


function  [dataCell, timeCell] = sliceTransformFn(sliceArgs, dataCell, timeCell, varargin)
    p = inputParser();
    p.addParameter('scalingApplied', true, @islogical);
    p.addParameter('slice', {}, @iscell);
    p.parse(varargin{:});
    
    % slice into the original sliceArgs
    slice = p.Results.slice;
    if ~isempty(slice)
        for iD = 1:numel(slice)
            sliceArgs{iD} = sliceArgs{iD}(slice{iD});
        end
    end
    
    for iT = 1:numel(dataCell)
        if isempty(dataCell{iT}), continue; end
        dataCell{iT} = dataCell{iT}(:, sliceArgs{:});
    end
    
end