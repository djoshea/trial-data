classdef ParamChannelDescriptor < ChannelDescriptor
    
    methods
        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.elementTypeByField = cd.UNKNOWN;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};
        end
        
        function type = getType(~)
            type = 'param';
        end

        function str = describe(cd)
            str = sprintf('Param %s (%s)', cd.name, cd.unitsPrimary);  
        end

        function cd = inferAttributesFromData(cd, varargin)
            % THIS ASSUMES THAT DATACELL IS HOMOGENOUS
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');
            
            if isempty(varargin)
                error('Must provide at least 1 data cell');
            end
           
            assert(numel(varargin) == 1, 'ParamChannels take only 1 data cell');
            dataCell = varargin{1};
            
            cd.dataFields = {cd.name};
            if isempty(cd.unitsByField)
                cd.unitsByField = {''};
            end
            
            if ~iscell(dataCell)
                cd.originalDataClassByField = {class(dataCell)};               
                if islogical(dataCell)
                    cd.elementTypeByField = cd.BOOLEAN;
                else
                    cd.elementTypeByField = cd.SCALAR;
                end

            else
                scalar = all(cellfun(@isscalar, dataCell));
                vector = all(cellfun(@isvector, dataCell));
                numeric = all(cellfun(@isnumeric, dataCell));
                cls = ChannelDescriptor.getCellElementClass(dataCell);
                
                cd.originalDataClassByField = cls;

                if scalar
                    if strcmp(cls, 'logical')
                        cd.elementTypeByField = cd.BOOLEAN;
                    else
                        cd.elementTypeByField = cd.SCALAR;
                    end

                elseif vector
                    cd.elementTypeByField = cd.VECTOR;

                elseif numeric
                    cd.elementTypeByField = cd.NUMERIC;

                elseif strcmp(cls, 'char')
                    cd.elementTypeByField = cd.STRING;

                else
                    warning('Inconsistent data types encountered');
                    cd.elementTypeByField = cd.UNKNOWN;
                    cd.originalDataClassByField = '';
                end 
            end
        end
        
    end
    
    methods(Static)
        function cd = buildStringParam(name)
            cd = ParamChannelDescriptor(name);
            cd.originalDataClassByField = {'char'};
            cd.elementTypeByField = cd.STRING;
        end
        
        function cd = buildScalarParam(name, units)
            cd = ParamChannelDescriptor(name);
            if nargin > 1
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            end
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.SCALAR;
        end 
        
        function cd = buildDatenumParam(name)
            cd = ParamChannelDescriptor(name);
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.DATENUM;
        end 
        
        function cd = buildBooleanParam(name)
            cd = ParamChannelDescriptor(name);
            cd.originalDataClassByField = {'logical'};
            cd.elementTypeByField = cd.BOOLEAN;
        end 
        
        function cd = buildFromValues(name, values)
            cd = ParamChannelDescriptor(name);
            cd = cd.inferAttributesFromData(values);
        end
    end
end
