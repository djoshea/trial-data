classdef ContinuousNeuralChannelGroupDescriptor < AnalogChannelGroupDescriptor
    properties(Dependent) 
        % these are encoded in the channel name
        array
        type
    end
    
    properties
        electrodes
    end
    
    methods
        function cd = ContinuousNeuralChannelGroupDescriptor(name, timeField, electrodes, varargin)
            sampleSize = numel(electrodes);
            cd = cd@AnalogChannelGroupDescriptor(name, timeField, sampleSize, varargin{:});
            cd.electrodes = makecol(electrodes);
        end
        
        function impl = getImpl(cd)
            impl = ContinuousNeuralChannelGroupImpl(cd);
        end
        
        function names = getSubChannelNames(cd, ~) 
            % allows sub classes to override
            names = arrayfun(@(elec) ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(cd.type, cd.array, elec, cd.nChannels), cd.electrodes);
        end
    end
    
    % these methods simply wrap the AnalogChannelDescriptor ones to get the
    % class name right
    methods(Static)
        function cd = buildAnalogGroup(name, timeField, electrodes, units, timeUnits, varargin)
            cd = ContinuousNeuralChannelGroupDescriptor(name, timeField, electrodes);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroup(name, timeField, units, timeUnits, ...
                'channelDescriptor', cd, 'sampleSize', numel(electrodes), varargin{:});
            
        end
        
        function cd = buildAnalogGroupFromValues(name, timeField, electrodes, units, timeUnits, dataCell, timeCell, varargin)
            cd = ContinuousNeuralChannelGroupDescriptor(name, timeField, electrodes);
            cd = AnalogChannelGroupDescriptor.buildAnalogGroupFromValues(name, timeField, ...
                units, timeUnits, dataCell, timeCell, ...
                'channelDescriptor', cd, 'sampleSize', numel(electrode), varargin{:});
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
            type = ContinuousNeuralChannelGroupDescriptor.parseTypeArray(cd.name);
        end
        
        function array = get.array(cd)
           [~, array] = ContinuousNeuralChannelGroupDescriptor.parseTypeArray(cd.name);
        end
       
        function cdIndividual = buildIndividualSubChannel(cd, index)
            units = cd.subChannelUnits{index};
            name = cd.subChannelNames{index};
            cdIndividual = ContinuousNeuralChannelDescriptor.buildSharedMatrixColumnAnalog(name, cd.name, index, cd.timeField, units, cd.unitsByField{2}, ...
                'scaleFromLims', cd.scaleFromLims, 'scaleToLims', cd.scaleToLims, ...
                'dataClass', cd.originalDataClassByField{1}, 'timeClass', cd.originalDataClassByField{2}, ...
                'transformChannelNames', cd.transformChannelNames, 'transformFn', cd.transformFn, ...
                'transformFnMode', cd.transformFnMode);
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
        
        function cls = getSubChannelClass()
            cls = 'ContinuousNeuralChannelDescriptor';
        end
    end

end