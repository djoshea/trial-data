function target = resolveSymlink(file)

    if ismac
        target = file;
        return;
    else
        target = getResolved(file);
    end
    
    return;

    if ~exist(file, 'file')
        % try recursively on its parent
        [parent leaf ext] = fileparts(file);
        parent = TrialDataUtilities.Save.resolveSymlink(parent);
        target = fullfile(parent, [leaf ext]);
    else
        target = getResolved(file);
    end

end

function path = escapePathForShell(path)
% path = escapePathForShell(path)
% Escape a path to a file or directory for embedding within a shell command
% passed to cmd or unix.

path = strrep(path, ' ', '\ ');

end

function result = getResolved(file)
    cmd = sprintf('readlink -m %s', escapePathForShell(GetFullPath(file)));
    [status, result] = system(cmd);
    if status || isempty(result)
        fprintf(result);
        error('Error resolving symbolic link');
    end 

    if result(end) == newline
        result = result(1:end-1);
    end
end
