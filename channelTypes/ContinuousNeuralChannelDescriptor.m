classdef ContinuousNeuralChannelDescriptor < AnalogChannelDescriptor
    properties(Dependent) 
        % these are encoded in the channel name as type_arrayElectrode,
        % e.g. lfp_vA001
        type string
        array string
        electrode
    end
    
    methods
        function cd = ContinuousNeuralChannelDescriptor(name, timeField)
            cd = cd@AnalogChannelDescriptor(name, timeField);
        end
        
        function impl = getImpl(cd)
            impl = ContinuousNeuralChannelImpl(cd);
        end
    end
    
    % these methods simply wrap the AnalogChannelDescriptor ones to get the
    % class name right
    methods(Static)
        function cd = buildVectorAnalog(name, timeField, units, timeUnits, varargin)
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
            cd = AnalogChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits, ...
                'channelDescriptor', cd, varargin{:});
        end
        
        function cd = buildSharedMatrixColumnAnalog(name, dataFieldName, dataFieldColumnIndex, ...
                timeField, units, timeUnits, varargin)
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
            cd = AnalogChannelDescriptor.buildSharedMatrixColumnAnalog(name, dataFieldName, dataFieldColumnIndex, ...
                timeField, units, timeUnits, ...
                'channelDescriptor', cd, varargin{:});
        end
        
        function cd = buildVectorAnalogFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
            cd = AnalogChannelDescriptor.buildVectorAnalogFromValues(name, timeField, ...
                units, timeUnits, dataCell, timeCell, ...
                'channelDescriptor', cd, varargin{:});
        end
    end
    
    methods
        function name = getNameWithUpdatedArray(cd, array)
            name = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(cd.type, array, cd.electrode);
        end
        
        function name = getNameWithUpdatedElectrode(cd, electrode)
            name = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(cd.type, cd.array, electrode);
        end
        
        function name = getNameWithUpdatedType(cd, type)
            name = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(type, cd.array, cd.electrode);
        end
      
        function type = get.type(cd)
            type = ContinuousNeuralChannelDescriptor.parseTypeArrayElectrode(cd.name);
        end
        
        function array = get.array(cd)
           [~, array] = ContinuousNeuralChannelDescriptor.parseTypeArrayElectrode(cd.name);
        end
        
        function elec = get.electrode(cd)
            [~, ~, elec] = ContinuousNeuralChannelDescriptor.parseTypeArrayElectrode(cd.name);
        end
        
        function cdGroup = buildGroupChannelDescriptor(cd)
            cdGroup = ContinuousNeuralChannelGroupDescriptor.buildAnalogGroup(cd.primaryDataField, cd.timeField, ...
                cd.unitsByField{1}, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2});
        end
    end
    
    methods(Static)
        function cd = buildFromTypeArrayElectrode(type, array, electrode, timeField, numElectrodes, units, timeUnits, varargin)
            if nargin < 5, numElectrodes = []; end
            name = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(type, array, electrode, numElectrodes);
            cd = ContinuousNeuralChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits, varargin{:});
        end
        
        function cd = buildFromArrayElectrode(arrayElectrodeStr, timeField, units, timeUnits, varargin)
            name = ContinuousNeuralChannelDescriptor.convertArrayElectrodeToChannelName(arrayElectrodeStr);
            cd = ContinuousNeuralChannelDescriptor.buildVectorAnalog(name, timeField, units, timeUnits, varargin{:});
        end
        
        function [chName] = convertArrayElectrodeToChannelName(arrayElectrodeStr)
            [type, array, electrode] = ContinuousNeuralChannelDescriptor.parseTypeArrayElectrode(arrayElectrodeStr);
            assert(~isnan(electrode), 'Could not parse array/electrode string %s', arrayElectrodeStr);
            chName = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(type, array, electrode, numElectrode);
        end
        
        function str = generateNameFromTypeArrayElectrode(type, array, electrode, numElectrodes)
            if nargin < 4 || isempty(numElectrodes)
                padZeros = 3;
            else
                padZeros = max(3, ceil(log10(numElectrodes)));
            end
            
            assert(isscalar(electrode));
            
           if isempty(type)
                str = sprintf("%s%0*d", array, padZeros, electrode);
            else
                str = sprintf("%s_%s%0*d", type, array, padZeros, electrode);
            end
        end

        function [type, array, electrode] = parseTypeArrayElectrode(str)
            tokens = regexp(str, '(?<type>\w+_)?(?<array>[A-Za-z_]*)(?<electrode>\d+)', 'names', 'once');
            if isempty(tokens)
                array = "";
                electrode = NaN;
                type = "";
            else
                if isempty(tokens.array) && ~isempty(tokens.type)
                    % this can happen when array ends with _
                    if contains(tokens.type, '_')
                        idx = find(tokens.type == '_', 1, 'first');
                        temp = tokens.type;
                        tokens.type = temp(1:idx);
                        tokens.array = temp(idx+1:end);
                    end
                end
                array = tokens.array;
                electrode = str2double(tokens.electrode);
                type = tokens.type(1:end-1); % strip trailing _
            end
        end
        
        function cls = getSubChannelClass()
            cls = '';
        end
    end

end