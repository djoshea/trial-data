function setArrowOpacity(s, faceAlpha, edgeAlpha)
% stores information in UserData struct to cause saveFigure to render
% marker points as translucent when exporting to svg

    if nargin < 3
        edgeAlpha = faceAlpha;
    end
    
    for i = 1:length(s)
        if ~verLessThan('matlab', '8.4')
            for iH = 1:numel(s.NodeChildren)
                h = s.NodeChildren(iH);
                switch(class(h))
                    case 'matlab.graphics.primitive.world.Marker'
                        if ~isempty(h.EdgeColorData)
                            h.EdgeColorType = 'truecoloralpha';
                            h.EdgeColorData(4) = uint8(edgeAlpha*255);
                        end
                        if ~isempty(h.FaceColorData)
                            h.FaceColorType = 'truecoloralpha';
                            h.FaceColorData(4) = uint8(faceAlpha*255);
                        end
% 
%                         % keep transparent
%                         addlistener(h, 'MarkedClean', @(ObjH, EventData) keepMarkerAlpha(ObjH, EventData, alpha));

                    case 'matlab.graphics.shape.internal.arrow.ArrowHead'
                        h.FaceAlpha = faceAlpha;
                        h.EdgeHandle.ColorType = 'truecoloralpha';
                        h.EdgeHandle.ColorData(4) = uint8(edgeAlpha*255);
                        % keep transparent
%                         addlistener(h, 'MarkedClean', @(ObjH, EventData) keepArrowHeadAlpha(ObjH, EventData, alpha));

                    case 'matlab.graphics.primitive.world.LineStrip'
                        h.ColorType = 'truecoloralpha';
                        h.ColorData(4) = uint8(edgeAlpha*255);
                end
            end
        end
    end
end

function keepMarkerAlpha(src, ~, faceAlpha, edgeAlpha)  
    mh = src;
    if ~isempty(mh.EdgeColorData)
        mh.EdgeColorType = 'truecoloralpha';
        mh.EdgeColorData(4) = uint8(edgeAlpha*255);
    end
    if ~isempty(mh.FaceColorData)
        mh.FaceColorType = 'truecoloralpha';
        mh.FaceColorData(4) = uint8(faceAlpha*255);
    end
end

function keepArrowHeadAlpha(src, ~, faceAlpha, edgeAlpha)  
    src.FaceAlpha = faceAlpha;
    if ~isempty(src.EdgeHandle)
        src.EdgeHandle.ColorType = 'truecoloralpha';
        src.EdgeHandle.ColorData(4) = uint8(edgeAlpha*255);
    end
end

