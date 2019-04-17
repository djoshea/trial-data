classdef ChannelData
    
properties
    cd
    impl
    conditionInfo
    alignInfoSet

    fieldData % nFields x 1 cell of data, each of which should be nTrials x ...
    fieldIds
end

methods
    function chdata = ChannelData()
        
    end
end

end