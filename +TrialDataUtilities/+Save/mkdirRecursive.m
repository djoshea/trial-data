function mkdirRecursive(dirPath)
% like mkdir -p : creates intermediate directories as required

s = warning('off', 'MATLAB:MKDIR:DirectoryExists');

if exist(dirPath, 'dir')
    return;
else
    parent = fileparts(dirPath);
    if ~isempty(parent)
        TrialDataUtilities.Save.mkdirRecursive(parent);
    end

    mkdir(dirPath);
end

warning(s);

