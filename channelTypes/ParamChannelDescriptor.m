classdef ParamChannelDescriptor < ChannelDescriptor

    methods
        function type = getType(cdesc)
            type = 'param';
        end

        function str = describe(cdesc)
            str = sprintf('Param (%s)', cdesc.name, cdesc.units);  
        end

        function dataFields = getExtraDataFields(cdesc)
            dataFields = {};
        end

        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.defaultValue = NaN;
            cd.scalar = true; % by default, change this if not true
        end
    end

    methods(Static) % infer channel descriptor from values
        function [cd cleanedValues] = inferFromValues(values)
            cd = ParamChannelDescriptor();
            assert(isvector(values), 'Values must be a vector');
            
            if ~iscell(values)
                cd.scalar = true;
            else
                % TODO deal with numeric vector type
                [cd.scalar mat] = isScalarCell(values);
                if cd.scalar
                    values = mat;
                end
            end

            cd.storageDataClass = class(values);
        end

    end
end
