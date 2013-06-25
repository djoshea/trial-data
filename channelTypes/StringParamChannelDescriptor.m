classdef StringParamChannelDescriptor < ChannelDescriptor

    methods
        function type = getType(cdesc)
            type = 'stringParam';
        end

        function str = describe(cdesc)
            str = sprintf('StringParam');  
        end

        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});

            cd.dataClass = 'char';
            cd.storageDataClass = 'char';
            cd.defaultValue = '';
        end
    end

end
