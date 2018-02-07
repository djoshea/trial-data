classdef ContinuousNeuralChannelGroupDescriptor < AnalogChannelGroupDescriptor
    properties(Dependent) 
        % these are encoded in the channel name
        array
        electrode
        type
    end
    
    methods
        function cd = ContinuousNeuralChannelGroupDescriptor(name, timeField, varargin)
            cd = cd@AnalogChannelGroupDescriptor(name, timeField, varargin{:});
        end
    end
    
    % these methods simply wrap the AnalogChannelDescriptor ones to get the
    % class name right
    methods(Static)
        function cd = buildAnalogGroup(name, timeField, units, timeUnits, varargin)
            cd = ContinuousNeuralChannelGroupDescriptor(name, timeField);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroup(name, timeField, units, timeUnits, ...
                'channelDescriptor', cd, varargin{:});
        end
        
        function cd = buildAnalogGroupFromValues(name, timeField, units, timeUnits, dataCell, timeCell, varargin)
            cd = ContinuousNeuralChannelGroupDescriptor(name, timeField);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroupFromValues(name, timeField, ...
                units, timeUnits, dataCell, timeCell, ...
                'channelDescriptor', cd, varargin{:});
        end
    end
    
    methods
        function name = getNameWithUpdatedArray(cd, array)
            name = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArrayElectrode(cd.type, array, cd.electrode);
        end
        
        function name = getNameWithUpdatedElectrode(cd, electrode)
            name = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArrayElectrode(cd.type, cd.array, electrode);
        end
        
        function name = getNameWithUpdatedType(cd, type)
            name = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArrayElectrode(type, cd.array, cd.electrode);
        end
      
        function type = get.type(cd)
            type = ContinuousNeuralChannelGroupDescriptor.parseTypeArrayElectrode(cd.name);
        end
        
        function array = get.array(cd)
           [~, array] = ContinuousNeuralChannelGroupDescriptor.parseTypeArrayElectrode(cd.name);
        end
        
        function elec = get.electrode(cd)
            [~, ~, elec] = ContinuousNeuralChannelGroupDescriptor.parseTypeArrayElectrode(cd.name);
        end
       
        function cdIndividual = buildIndividualSubChannel(cd, name, index)
            cdIndividual = ContinuousNeuralChannelDescriptor.buildSharedMatrixColumnAnalog(name, cd.name, index, cd.timeField, cd.unitsByField{1}, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2});
        end
    end
    
    methods(Static)
        function cd = buildFromTypeArray(type, array, varargin)
            name = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArray(type, array);
            timeField = sprintf('%s_time', name);
            cd = ContinuousNeuralChannelGroupDescriptor.buildAnalogGroup(name, timeField, varargin{:});
        end
        
        function cd = buildFromArray(array, varargin)
            name = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArray('', array);
            timeField = sprintf('%s_time', name);
            cd = ContinuousNeuralChannelGroupDescriptor.buildAnalogGroup(name, timeField, varargin{:});
        end
        
        function [chName] = convertTypeArrayToChannelName(array)
            [type, array] = ContinuousNeuralChannelGroupDescriptor.parseTypeArray(array);
             chName = ContinuousNeuralChannelGroupDescriptor.generateNameFromTypeArray(type, array);
        end
        
        function str = generateNameFromTypeArray(type, array)
            if isempty(type)
                str = array;
            else
                if isempty(array)
                    array = 'array';
                end
                str = sprintf('%s_%s', type, array);
            end
        end

        function [type, array] = parseTypeArray(str)
            tokens = regexp(str, '(?<type>\w+_)?(?<array>[A-Za-z_]*)', 'names', 'once');
            if isempty(tokens)
                array = '';
                type = '';
            else
                array = tokens.array;
                type = tokens.type(1:end-1); % strip trailing _
            end
        end
    end

end