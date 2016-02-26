function map = expandWrapColormap(map, n)
    if isa(map, 'function_handle')
        map = map(n);
    elseif size(map, 1) < n
        % make at least big enough
        map = repmat(map, ceil(n / size(map, 1)), 1);
%         % cut it down to size - actually no need
%         map = map(1:n, :);
    end
end