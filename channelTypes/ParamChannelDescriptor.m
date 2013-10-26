classdef ParamChannelDescriptor < ChannelDescriptor
    
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
            
            if ~iscell(dataCell)
                cd.dfd = ScalarField();
                cd.storageDataClass = class(dataCell);
                
            else
                v = []; i = 1;
                while isempty(v)
                    v = dataCell{i};
                    i = i+1;
                end
                
                if ischar(v)
                    cd.dfd = StringField();

                elseif isscalar(v)
                    if islogical(v)
                        cd.dfd = BooleanField();
                    else
                        cd.dfd = ScalarField();
                    end

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
        
    end
    
    methods(Static)
        function cd = buildStringParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd = StringField();
            cd.storageDataClass = 'char';
        end
        
        function cd = buildScalarParam(name, units)
            cd = ParamChannelDescriptor(name);
            if nargin > 1
                cd.units = units;
            end
            cd.dfd = ScalarField();
            cd.storageDataClass = 'double';
        end 
        
        function cd = buildDatenumParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd  = DateTimeField();
            cd.storageDataClass = 'double';
        end 
        
        function cd = buildBooleanParam(name)
            cd = ParamChannelDescriptor(name);
            cd.dfd  = BooleanField();
            cd.storageDataClass = 'logical';
        end 
        
        function cd = buildFromValues(name, values)
            cd = ParamChannelDescriptor(name);
            cd = cd.inferAttributesFromData(values);
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
