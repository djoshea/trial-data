classdef ChannelFieldSpec 
    % Utility class for identifying the type and structure of data within a single field belonging to a channel
    % Instances are built by ChannelDescriptors to describe their own data, and 
    % methods of this class are used by TrialData to validate, convert, and clean data corresponding to that field
    
    properties
        elementClass
        
         % nFields x 1 arrays or cells
        collectAsCellByField
        missingValueByField
        isNumericScalarByField % true for logical as well
        isBooleanByField
        isStringByField
        isVectorByField
        isVectorizableByField % elements may be concatenated into vector
        isCategoricalByField

        % indicates whether a given field is shareable between multiple
        % channels. If a field is marked as shareable, it will be copied
        % when this channel's data is updated, so that the other channel's
        % are not affected. See getIsShareableByField for implementation
        isShareableByField

        % data class of each field as accessed
        accessClassByField % cell array specifying data class to convert each

        % data class of each field as stored in .data
        memoryClassByField

        % persistent storage data class of each field (or empty if unknown / mixed)
        storageClassByField

        % these are shortcuts to the first element of the properties above
        % since the first data field is considered the primary data field
        unitsPrimary
        dataFieldPrimary
