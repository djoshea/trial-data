function hvec = allocateGraphicsHandleVector(num)
    if verLessThan('matlab','8.4.0')
        hvec = nan(num, 1);
    else
        hvec = gobjects(num, 1);
    end
end
        
        