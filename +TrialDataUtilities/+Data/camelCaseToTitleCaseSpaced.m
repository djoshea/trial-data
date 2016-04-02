function newName = camelCaseToTitleCaseSpaced(name)
% here we assume the name is TitleCased, camelCased, or
% snake_cased and convert to Spaced Words
    pattern = '([A-Z]*[a-z]+)';
    words = regexp(name, pattern, 'match');
    
    if iscell(name)
        newName = cellfun(@(word) strjoin(upperFirst(word), ' '), words, 'UniformOutput', false);
    else
        newName = strjoin(upperFirst(words), ' ');
    end
end