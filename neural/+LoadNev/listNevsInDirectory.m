function fileList = listNevsInDirectory(directory, varargin)
% fileList = listNevsInDirectory(directory, varargin)
% 
% Returns a list of fully-qualified nev file names in directory, ordered according to their numerical suffix
% If a nev of equivalent name exists within sortedNevDirectory (inside the dataContext), 
%   the path to this sorted nev file

p = inputParser();
p.addParameter('excludeFilenameContains', '', @ischar);
p.parse(varargin{:});

% order the nevs according to their numerical suffix before loading
fileList = makecol(sort(listFilesInDirectory(directory, 'nev')));

if ~isempty(p.Results.excludeFilenameContains)
    where = strfind(fileList, p.Results.excludeFilenameContains);
    mask = cellfun(@isempty, where);
    fileList = fileList(mask);
end
