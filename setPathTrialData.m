function setPathTrialData()
    codeRoot = pathToThisFile();
    fprintf('Path: Adding trial-data at %s\n', codeRoot);
    addPathRecursive(codeRoot);
end
