function desc = get.attributeDescriptions(ci)
    desc = cellvec(ci.nAttriubutes);
    modes = ci.getAttribute
    for i = 1:ci.nAttributes
        name = ci.attributeNames{i};  
        nValues = ci.nValuesByAttribute(i);
        nAutoBins = ci.attributeValueBinsAutoCount(i);

        switch modes{i}
            case AttributeValueListManual
                suffix = sprintf('(%d)', nValues);
            case AttributeValueListAuto
                suffix = '(?)';
            case AttributeValueBinsManual
                suffix = sprintf('(%d bins)', nValues);
            case AttributeValueBinsAutoUniform
                suffix = sprintf('(%d unif bins)', nAutoBins);
            case AttributeValueBinsAutoQuantile
                suffix = sprintf('(%d quantiles)', nAutoBins);
        end

        desc{i} = sprintf('%s %s', name, suffix);
    end
end

function desc = get.axisDescriptions(ci)
    desc = cellvec(ci.nAxes);

    for i = 1:ci.nAxes
        attr = ci.axisAttributes{i}; 
        idx = ci.getAttributeIdx(attr);   
        attrDesc = ci.attributeDescriptions(idx);
        nValuesAttr = ci.nValuesByAttribute(idx);

        desc{i} = strjoin(attrDesc, ' x '); 
    end
end
