classdef ParamChannelDescriptor < ChannelDescriptor
    methods(Access=protected)
        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.elementTypeByField = cd.UNKNOWN;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};
            cd.fieldIds = {'data'};
        end
    end
    
    methods
        function cd = initialize(cd)
            % check and repair internal consistency
            cd.warnIfNoArgOut(nargout);
            if isempty(cd.fieldIds)
                cd.fieldIds = {'data'};
            end
            
            cd = initialize@ChannelDescriptor(cd);
        end
        
        function impl = getImpl(cd)
            impl = ParamChannelImpl(cd);
        end
        
        function type = getType(~)
            type = 'param';
        end

        function str = describe(cd)
            if isempty(cd.unitsPrimary)
                str = sprintf('Param %s', cd.name);  
            else
                str = sprintf('Param %s (%s)', cd.name, cd.unitsPrimary);  
            end
        end

        function cd = inferAttributesFromData(cd, varargin)
            % THIS ASSUMES THAT DATACELL IS HOMOGENOUS
            cd.warnIfNoArgOut(nargout);
            
            if isempty(varargin)
                error('Must provide at least 1 data cell');
            end
           
            assert(numel(varargin) == 1, 'ParamChannels take only 1 data cell');
            dataCell = varargin{1};
            
            cd.dataFields = {cd.name};
            if isempty(cd.unitsByField)
                cd.unitsByField = {''};
            end
            
            if strcmp(cd.name, 'protocolVersion')
                a = 1;
            end
            
            if ~iscell(dataCell)
                data = dataCell;
                cd.originalDataClassByField = {class(data)};               
                if islogical(data)
                    cd.elementTypeByField = cd.BOOLEAN;
                elseif isstring(data)
                    cd.elementTypeByField = cd.STRING;
                elseif iscategorical(data)
                    cd.elementTypeByField = cd.CATEGORICAL;
                else
                    cd.elementTypeByField = cd.SCALAR;
                end

            else
                scalar = all(cellfun(@(x) isempty(x) || isscalar(x), dataCell));
                vector = all(cellfun(@(x) isempty(x) || isvector(x), dataCell));
                numeric = all(cellfun(@isnumeric, dataCell));
                cls = ChannelImpl.getCellElementClass(dataCell);
                
                cd.originalDataClassByField = {cls};

                if strcmp(cls, 'cell')
                    cd.elementTypeByField = cd.CELL;
                    
                elseif strcmp(cls, 'char')
                    cd.elementTypeByField = cd.STRING;
                
                elseif scalar
                    catData = cat(1, dataCell{:});
                    if strcmp(cls, 'logical') || (strcmp(cls, 'uint8') && all(ismember(catData, uint8([0 1]))))
                        cd.elementTypeByField = cd.BOOLEAN;
                    else
                        cd.elementTypeByField = cd.SCALAR;
                    end
                    
                elseif vector
                    cd.elementTypeByField = cd.VECTOR;

                elseif numeric
                    cd.elementTypeByField = cd.NUMERIC;

                else
                    warning('Inconsistent data types encountered');
                    cd.elementTypeByField = cd.UNKNOWN;
                    cd.originalDataClassByField = {''};
                end 
            end
        end
        
    end
    
    methods(Static)
        function cd = buildLike(cdTemplate, name)
            cd = cdTemplate;
            cd.name = name;
            cd.dataFields{1} = name;
        end
        
        function cd = buildStringParam(name, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            cd.originalDataClassByField = {'char'};
            cd.elementTypeByField = cd.STRING;
        end
        
        function cd = buildScalarParam(name, dataClass, units, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            if nargin < 2
                dataClass = 'double';
            end
            if nargin > 2
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {dataClass};
            cd.elementTypeByField = cd.SCALAR;
        end
        
        function cd = buildVectorParamAccessAsMatrix(varargin)
            cd = ParamChannelDescriptor.buildVectorParam(varargin{:});
            cd.catAlongFirstDimByField = true;
        end 
        
        function cd = buildVectorParam(name, dataClass, units, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            if nargin < 2
                dataClass = 'double';
            end
                
            if nargin > 2
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {dataClass};
            cd.elementTypeByField = cd.VECTOR;
        end 
        
        function cd = buildNumericParam(name, dataClass, units, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            if nargin < 2
                dataClass = 'double';
            end
            if nargin > 2
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {dataClass};
            cd.elementTypeByField = cd.NUMERIC;
        end 
        
        function cd = buildDatenumParam(name, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.DATENUM;
        end 
        
        function cd = buildBooleanParam(name, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            cd.originalDataClassByField = {'logical'};
            cd.elementTypeByField = cd.BOOLEAN;
        end 
        
        function cd = buildFromValues(name, values, units, varargin)
            cd = ParamChannelDescriptor(name, varargin{:});
            cd = cd.inferAttributesFromData(values);
            if nargin < 3
                units = '';
            end
            cd.unitsByField = {units};
        end
        
        function cls = getSubChannelClass()
            cls = 'ChannelDescriptor';
        end
    end
end
