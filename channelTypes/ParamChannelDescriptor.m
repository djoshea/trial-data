classdef ParamChannelDescriptor < ChannelDescriptor
    
    properties
        isString = false;
        isScalar = false;
        isBoolean = false;
    end
    
    methods
        function type = getType(cdesc)
            type = 'param';
        end

        function fields = getDataFields(cdesc)
            fields = {cdesc.name};
        end

        function str = describe(cdesc)
            str = sprintf('Param (%s)', cdesc.name, cdesc.units);  
        end

        function cd = ParamChannelDescriptor(varargin)
            cd = cd@ChannelDescriptor(varargin{:});
        end
        
        function cd = inferAttributesFromData(cd, dataCell)
            % THIS ASSUMES THAT DATACELL IS HOMOGENOUS
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');
            
            v = []; i = 1;
            while isempty(v)
                v = dataCell{i};
                i = i+1;
            end
            
            cd.isString = false;
            cd.isScalar = false;
            if ischar(v)
                cd.dfd = StringField();
                cd.isString = true;
                
            elseif isscalar(v)
                if islogical(v)
                    cd.dfd = BooleanField();
                    cd.isBoolean = true;
                else
                    cd.dfd = ScalarField();
                end
                cd.isScalar = true;
                
            elseif isvector(v)
                cd.dfd = NumericVectorField();
              
            elseif isnumeric(v) || islogical(v)
                cd.dfd = NumericField();
            else
                error('ParameterChannel attributes could not be inferred from data');
            end
            
            cd.storageDataClass = class(v);
        end
        
    end
    
    methods(Static)
        function cd = buildStringParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd = StringField();
            cd.isString = true;
            cd.storageDataClass = 'char';
        end
        
        function cd = buildScalarParam(name, units)
            cd = ParamChannelDescriptor(name);
            if nargin > 1
                cd.units = units;
            end
            cd.dfd = ScalarField();
            cd.isScalar = true;
            cd.storageDataClass = 'double';
        end 
        
        function cd = buildDatenumParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd  = DateTimeField();
            cd.isScalar = true;
            cd.storageDataClass = 'double';
        end 
        
        function cd = buildBooleanParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd  = BooleanField();
            cd.isScalar = true;
            cd.isBoolean = true;
            cd.storageDataClass = 'logical';
        end 
    end
 
%     methods(Static) % infer channel descriptor from values
%         function [cd cleanedValues] = inferFromValues(values)
%             cd = ParamChannelDescriptor();
%             assert(isvector(values), 'Values must be a vector');
%             
%             if ~iscell(values)
%                 cd.scalar = true;
%             else
%                 % TODO deal with numeric vector type
%                 [cd.scalar mat] = isScalarCell(values);
%                 if cd.scalar
%                     values = mat;
%                 end
%             end
% 
%             cd.storageDataClass = class(values);
%         end
% 
%     end
end
