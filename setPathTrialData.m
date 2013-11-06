function setPathTrialData()
    codeRoot = pathToThisFile();
    fprintf('Path: Adding trial-data code at %s\n', codeRoot);
    addPathRecursive(codeRoot);
end
