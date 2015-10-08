function startIdx = findFirstConsecutiveRun(v, nConsecutive, dim, nToReturn)
% finds the start locations of consecutive runs of nConsecutive true values
% in a row

    if nargin < 3
        dim = find(size(v) ~= 1, 1);
        if isempty(dim), dim = 1; end
    end
    if nargin < 4
        nToReturn = 1;
    end
    
    % pad on either side with zeros
    sz = size(v);
    szPad = sz;
    szPad(dim) = 1;
    zeroPad = cast(zeros(szPad), 'like', v);
    vAux = cat(dim, zeroPad, v, zeroPad);
    dvAux = diff(vAux, 1, dim);
   
    
    szOut = sz;
    szOut(dim) = nToReturn;

    startIdxCell = TensorUtils.mapSlices(@findSingleLongRuns, dim, dvAux);
    startIdx = cell2mat(startIdxCell);
    
    function startIdxSelected = findSingleLongRuns(daux)
        runStartIdx = find(daux == 1);
        if isempty(runStartIdx)
            startIdxSelected = nanvec(nToReturn);
        else
            runLengths = find(daux == -1) - runStartIdx;
            idx = find(runLengths >= nConsecutive, nToReturn, 'first');
            startIdxSelected = cat(1, makecol(runStartIdx(idx)), nanvec(nToReturn-numel(idx)));        
        end
    end
end
