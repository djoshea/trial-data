function cmat = toRGB(c)
    if isnumeric(c)
        cmat = c;
        assert(size(c, 2) == 3 || size(c, 2) == 4, 'Color matrix is not the correct size');
    elseif iscell(c)
        cmat = cell2mat(cellfun(@convertSingle, makecol(c), 'UniformOutput', false));
    elseif ischar(c)
        cmat = convertSingle(c);
    end
    
    function cvec = convertSingle(c)
        if isnumeric(c)
            cvec = c;
            return;
        end
        switch c
            case 'b'
                cvec = [0 0 1];
            case 'g'
                cvec = [0 1 0];
            case 'r'
                cvec = [1 0 0];
            case 'c'
                cvec = [0 1 1];
            case 'm'
                cvec = [1 0 1];
            case 'y'
                cvec = [1 1 0];
            case 'k'
                cvec = [0 0 0];
            case 'w'
                cvec = [1 1 1];
            otherwise
                error('Unknown color string');
        end
    end
end
