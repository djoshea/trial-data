function [Qcell, infoCell] = loadNevMulti(fnameList, loadNevFn, varargin)

% check for existence of all file
for ifile = 1:length(fnameList)
    if ~exist(fnameList{ifile}, 'file')
        error('Could not find file %s', fnameList{ifile});
    end
end

% load each nev file
[Qcell, infoCell] = cellvec(numel(fnameList));
for ifile = 1:length(fnameList)
    [Qcell{ifile}, infoCell{ifile}] = loadNevFn(fnameList{ifile}, varargin{:});
end

end
