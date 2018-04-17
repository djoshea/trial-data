function [dataCell, tvec] = equalizeTimeVectorsForTimeseries(dataCell, timeCell, timeDim, varargin)
% given an N x 1 data cell and time cell, where time cell contains vectors
% that describe time along dataCell{:} dim timeDim, equalize the contents
% time vectors

N = numel(dataCell);
assert(numel(dataCell) == N, 'Sizes of dataCell and timeCell must match');

if timeDim ~= 1
    ndimsOrig = nanvec(N);
    for i = 1:numel(dataCell)
        [dataCell{i}, ndimsOrig(i)] = TensorUtils.shiftdimToFirstDim(dataCell{i}, timeDim);
        assert(isvector(timeCell{i}) && size(dataCell{i}, 1) == numel(timeCell{i}), 'timeCell must be time vector whose sizes match size(dimCell, timeDim)');
    end
end

[dataCat, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(dataCell(:), timeCell(:), varargin{:});
dataCat = TensorUtils.splitAlongDimension(dataCat, 1, [], true);

if timeDim ~= 1
    for i = 1:numel(dataCell)
        dataCell{i} = TensorUtils.unshiftdimToFirstDim(dataCat{i}, timeDim, ndimsOrig(i));
    end
end

