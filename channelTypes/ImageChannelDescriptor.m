classdef ImageChannelDescriptor < AnalogChannelGroupDescriptor
    properties
        fileName = '';
        fileInfo = struct();
    end
    
    methods
        function cd = ImageChannelDescriptor(name, timeField)
            cd = cd@AnalogChannelGroupDescriptor(name, timeField);
        end
        
        function impl = getImpl(cd)
            impl = ImageChannelImpl(cd);
        end
    end
    
    % these methods simply wrap the AnalogChannelDescriptor ones to get the
    % class name right
    methods(Static)
        function cd = buildImage(name, timeField, units, timeUnits, varargin)
            cd = ImageChannelDescriptor(name, timeField);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroup(name, timeField, units, timeUnits, ...
                'channelDescriptor', cd, varargin{:});
        end
        
        function cd = buildImageFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            cd = ImageChannelDescriptor(name, timeField);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroupFromValues(name, timeField, ...
                units, timeUnits, dataCell, timeCell, ...
                'channelDescriptor', cd, varargin{:});
        end
        
        function cls = getSubChannelClass()
            cls = '';
        end
    end
end