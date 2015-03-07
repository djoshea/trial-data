function hn = getNullGraphicsHandle()
    if verLessThan('matlab','8.4.0')
        hn = NaN;
    else
        hn = matlab.graphics.GraphicsPlaceholder();
    end
end
