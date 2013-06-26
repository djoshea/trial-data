classdef StringParamChannelDescriptor < ParamChannelDescriptor

    methods
        function type = getType(cdesc)
            type = 'stringParam';
        end

        function str = describe(cdesc)
            str = sprintf('StringParam');  
        end

        function cd = StringParamChannelDescriptor(varargin)
            cd = cd@ParamChannelDescriptor(varargin{:});

            cd.storageDataClass = 'char';
            cd.defaultValue = '';
        end
    end

end
