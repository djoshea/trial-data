function success = symlink(src, link)
% symlink(src, linkDest)

    src = TrialDataUtilities.Save.resolveSymlink(GetFullPath(src));
    assert(exist(src, 'file') == 2, 'Source is a dangling symlink');
    
    link = GetFullPath(link);
    mkdirRecursive(fileparts(link));
    if exist(link, 'file')
        delete(link);
    end
    cmd = sprintf('ln -s "%s" "%s"', src, link);
    [status, output] = unix(cmd);
    
    if status
        fprintf('Error creating symlink: \n');
        fprintf('%s\n', output);
    end

    success = ~status;
end
