function [dataCellProjected dataCellRemaining percentVarianceIn] = projectDataCellIntoSubspace(dataCell, vectors)

% given an nBases x nConditions x nAlign dataCell where all vectors have the same length,
% project onto subspace spanned by vector (column vectors)

nBases = size(dataCell, 1);
nConditions = size(dataCell, 2);
nAlign = size(dataCell, 3);

P = vectors*(vectors'*vectors)^(-1)*vectors';

[dataCellProjected, dataCellRemaining] = deal(cell(size(dataCell)));

for iAlign = 1:nAlign
    for iCondition = 1:nConditions
        dataMat = cell2mat(dataCell(:, iCondition, iAlign)')';
        dataProj = P*dataMat;
        dataRem = dataMat - dataProj;
        % dataProj/Rem are nUnits x nTimePoints
        
        dataCellProjected(:, iCondition, iAlign) = mat2cell(dataProj', size(dataProj, 2), ones(nBases, 1));
        dataCellRemaining(:, iCondition, iAlign) = mat2cell(dataRem', size(dataRem, 2), ones(nBases, 1));
    end
end

