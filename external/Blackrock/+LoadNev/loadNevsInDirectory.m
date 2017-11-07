function [Qcell, infoCell, fileList] = loadNevsInDirectory(directory, loadNevFn, varargin)

p = inputParser();
p.addParameter('excludeFilenameContains', '', @ischar);
p.KeepUnmatched = true;
p.parse(varargin{:});

fileList = LoadNev.listNevsInDirectory(directory, 'excludeFilenameContains', p.Results.excludeFilenameContains);
if isempty(fileList)
    fprintf('\tWarning: no .nev files found in %s\n', directory);
    Qcell = {};
    infoCell = {};
    fileList = {};
    return;
end

[Qcell, infoCell] = LoadNev.loadNevMulti(fileList, loadNevFn, p.Unmatched);

end
