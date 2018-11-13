function prefix = longestCommonPrefix(entries)
    entries = string(entries);
    nCommon = nan(numel(entries)-1, 1);
    for iS = 1:numel(entries)-1
        nCommon(iS) = getNCommon(entries(iS), entries(iS+1));
    end
    
    nCommon = min(nCommon);
    prefix = extractBefore(entries(1), nCommon+1);
end

function n = getNCommon(a, b)
    len = min(strlength(a), strlength(b));
    for i = 1:len
        if ~strncmp(a, b, i)
            n = i-1;
            return;
        end
    end
    
    n = len;
end
        