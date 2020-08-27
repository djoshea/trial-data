function tf = isstringlike(v)
    tf = ischar(v) || isstring(v) || iscellstr(v);
end