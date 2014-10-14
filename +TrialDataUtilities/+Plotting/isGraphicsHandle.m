function mask = isGraphicHandle(h)
    if isobject(h)
        mask = isgraphics(h);
    else
        mask = ~isnan(h);
    end
end