function setMarkerLayerFront(s)
    for i = 1:length(s)
        if ~verLessThan('matlab', '8.4')
            mh = s.MarkerHandle;
            if ~isa(mh, 'matlab.graphics.GraphicsPlaceholder')
                mh.Layer = 'front';
            else
                % defer until ready
                addlistener(s(i),'MarkedClean',...
                    @(h, EventData) setLayer(h));
            end
        end
    end
end

function setLayer(src)  
    src.MarkerHandle.Layer = 'front';
end

