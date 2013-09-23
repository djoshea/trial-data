classdef AnalogChannelDescriptor < ChannelDescriptor
    properties
        timeField
        
        isVector = true; % is this a vector analog field?
    end

    methods
        % default to name_time if not manually specified
        function field = get.timeField(cdesc)
            if isempty(cdesc.timeField)
                if isempty(cdesc.name)
                    field = '';
                else
                    field = sprintf('%s_time', cdesc.name);
                end
            else
                field = cdesc.timeField;
            end
        end

        function type = getType(cdesc)
            type = 'analog';
        end

        function fields = getDataFields(cdesc)
            fields = {cdesc.name, cdesc.timeField};
        end

        function str = describe(cdesc)
            str = sprintf('Analog (%s)', cdesc.name, cdesc.units);  
        end

        function cd = AnalogChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
        
        function cd = inferAttributesFromData(cd, dataCell)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            cd.isVector = true;
            for i = 1:numel(dataCell)
                if ~isvector(dataCell{i})
                    cd.isVector(false)
                    break;
                end
            end
                   
            if cd.isVector
                cd.dfd = NumericVectorField();
            else
                cd.dfd = NumericField();
            end
            
            if isempty(dataCell)
                cd.storageDataClass = 'double';
            else
                cd.storageDataClass = class(dataCell{1});
            end
        end
        
        function tf = getIsVector(cd)
            tf = cd.isVector;
        end
    end
    
     methods(Static)
        function cd = buildVectorAnalog(name, timeField, units)
            cd = AnalogChannelDescriptor(name);
            cd.timeField = timeField;
            cd.dfd = NumericVectorField();
            cd.isVector = true;
            if nargin > 2
                cd.units = units;
            end
        end
    end

end
