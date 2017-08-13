function dataCat = simpleCatJoin(dataPre, dataPost, joinIdxInPre, nextIdxInPost)
% time along dim 2, bases along dim 1, joinIdxInPre and nextIdxInPost must have size of dims 3 and
% beyond, and are usually returned by computeBestJointPoint

    nPre = size(dataPre, 2);
    nPost = size(dataPost, 2);
    nOverlap = nPre - joinIdxInPre + nextIdxInPost - 1;
    
    sz = size(dataPre);
    szCat = sz;
    szCat(2) = nPre + nPost - min(nOverlap(:));
    nTraj = prod(sz(3:end));
    
    dataCat = nan(szCat);
        
    for c = 1:nTraj
        insert = cat(2, dataPre(:, 1:joinIdxInPre(c), c), dataPost(:, nextIdxInPost(c):end, c));
        dataCat(:, 1:size(insert, 2), c) = insert;
    end
    
end
