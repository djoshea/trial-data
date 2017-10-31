function tf = isintegertol(mat)
    tf = abs(mat - round(mat)) < 1e-6;
end