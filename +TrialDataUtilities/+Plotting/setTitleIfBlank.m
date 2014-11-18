function h = setTitleIfBlank(axh, fmat, varargin)
    if isempty(varargin)
        str = fmat;
    else
        str = sprintf(fmat, varargin{:});
    end

    h = get(axh, 'Title');
    if isempty(h) || isempty(get(h, 'String'))
        h = title(str);
        set(h, 'Interpreter', 'none');
    end
end