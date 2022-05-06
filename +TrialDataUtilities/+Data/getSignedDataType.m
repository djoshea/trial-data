function outType = getSignedDataType(inType, signed)
    arguments
        inType (1, 1) string
        signed (1, 1) logical = true;
    end
    
    switch inType
        case {"uint8", "int8"}
            if signed
                outType = "int8";
            else
                outType = "uint8";
            end
        
        case {"uint16", "int16"}
            if signed
                outType = "int16";
            else
                outType = "uint16";
            end

        case {"uint32", "int32"}
            if signed
                outType = "int32";
            else
                outType = "uint32";
            end

        case {"uint64", "int64"}
            if signed
                outType = "int64";
            else
                outType = "uint64";
            end

        case {"single", "double"}
            outType = inType;

        otherwise
            error("Invalid numeric type %s", inType);
    end

end