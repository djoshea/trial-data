classdef ParamChannelDescriptor < ChannelDescriptor
    
    methods(Access=protected)
        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
            cd.dataFields = {cd.name};
            cd.elementTypeByField = cd.UNKNOWN;
            cd.originalDataClassByField = {''};
            cd.unitsByField = {''};
        end
    end
    
    methods
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
                cd.originalDataClassByField = {class(dataCell)};               
                if islogical(dataCell)
                    cd.elementTypeByField = cd.BOOLEAN;
                else
                    cd.elementTypeByField = cd.SCALAR;
                end

            else
                scalar = all(cellfun(@(x) isempty(x) || isscalar(x), dataCell));
                vector = all(cellfun(@(x) isempty(x) || isvector(x), dataCell));
                numeric = all(cellfun(@isnumeric, dataCell));
                cls = ChannelDescriptor.getCellElementClass(dataCell);
                
                cd.originalDataClassByField = {cls};

                if strcmp(cls, 'cell')
                    cd.elementTypeByField = cd.CELL;
                    
                elseif scalar
                    catData = cat(1, dataCell{:});
                    if strcmp(cls, 'logical') || (strcmp(cls, 'uint8') && all(ismember(catData, uint8([0 1]))))
                        cd.elementTypeByField = cd.BOOLEAN;
                    else
                        cd.elementTypeByField = cd.SCALAR;
                    end

                elseif strcmp(cls, 'char')
                    cd.elementTypeByField = cd.STRING;
                    
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
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.SCALAR;
        end
        
        function cd = buildVectorParamAccessAsMatrix(varargin)
            cd = ParamChannelDescriptor.buildVectorParam(varargin{:});
            cd.catAlongFirstDimByField = true;
        end 
        
        function cd = buildVectorParam(name, units)
            cd = ParamChannelDescriptor(name);
            if nargin > 1
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.VECTOR;
        end 
        
        function cd = buildNumericParam(name, units)
            cd = ParamChannelDescriptor(name);
            if nargin > 1
                assert(ischar(units), 'Units must be string');
                cd.unitsByField = {units};
            else
                cd.unitsByField = {''};
            end
            cd.originalDataClassByField = {'double'};
            cd.elementTypeByField = cd.NUMERIC;
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
        
        function cd = buildFromValues(name, values, units)
            cd = ParamChannelDescriptor(name);
            cd = cd.inferAttributesFromData(values);
            if nargin < 3
                units = '';
            end
            cd.unitsByField = {units};
        end
    end
end
