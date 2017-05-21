classdef ContinuousNeuralChannelDescriptor < AnalogChannelDescriptor
    properties(Dependent) 
        % these are encoded in the channel name
        array
        electrode
        type
    end
    
%     properties(SetAccess=protected)
%         arrayManual = '';
%         electrodeManual = [];
%         typeManual = '';
%     end
    
    methods
        function cd = ContinuousNeuralChannelDescriptor(name, timeField)
            cd = cd@AnalogChannelDescriptor(name, timeField);
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
        function cd = buildFromTypeArrayElectrode(type, array, electrode, timeField)
            name = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(type, array, electrode);
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
        end
        
        function cd = buildFromArrayElectrode(arrayElectrodeStr, timeField)
            name = ContinuousNeuralChannelDescriptor.convertArrayElectrodeToChannelName(arrayElectrodeStr);
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
        end
        
        function [chName] = convertArrayElectrodeToChannelName(arrayElectrodeStr)
            [type, array, electrode] = ContinuousNeuralChannelDescriptor.parseTypeArrayElectrode(arrayElectrodeStr);
            assert(~isnan(electrode), 'Could not parse array/electrode string %s', arrayElectrodeStr);
            chName = ContinuousNeuralChannelDescriptor.generateNameFromTypeArrayElectrode(type, array, electrode);
        end
        
        function str = generateNameFromTypeArrayElectrode(type, array, electrode)
            if isempty(type)
                str = sprintf('%s%03d', array, electrode);
            else
                str = sprintf('%s_%s%03d', type, array, electrode);
            end
        end

        function [type, array, electrode] = parseTypeArrayElectrode(str)
            tokens = regexp(str, '(?<type>\w+_)?(?<array>[A-Za-z_]*)(?<electrode>\d+)', 'names', 'once');
            if isempty(tokens)
                array = '';
                electrode = NaN;
                type = '';
            else
                array = tokens.array;
                electrode =str2double(tokens.electrode);
                type = tokens.type(1:end-1); % strip trailing _
            end
        end
    end

end