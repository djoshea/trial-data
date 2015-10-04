classdef ContinuousNeuralChannelDescriptor < AnalogChannelDescriptor
    properties(Dependent)
        array
        electrode
    end
    
    properties(SetAccess=protected)
        arrayManual = '';
        electrodeManual = [];
    end
    
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
        function cd = setArrayElectrode(cd, array, electrode)
            cd.warnIfNoArgOut(nargout);
            assert(ischar('array'));
            assert(isnumeric(electrode));
            cd.electrodeManual = electrode;
            cd.arrayManual = array;
        end
        
        function a = get.array(cd)
            if isempty(cd.arrayManual)
                [a, ~] = ContinuousNeuralChannelDescriptor.parseArrayElectrode(cd.name);
            else
                a = cd.arrayManual;
            end
        end
        
        function e = get.electrode(cd)
            if isempty(cd.electrodeManual) || isnan(cd.electrodeManual)
                [~, e] = ContinuousNeuralChannelDescriptor.parseArrayElectrode(cd.name);
            else
                e = cd.electrodeManual;
            end
        end
    end
    
    methods(Static)
        function cd = buildFromArrayElectrode(arrayElectrodeStr, timeField)
            name = ContinuousNeuralChannelDescriptor.convertArrayElectrodeToChannelName(arrayElectrodeStr);
            cd = ContinuousNeuralChannelDescriptor(name, timeField);
        end
        
        function [chName] = convertArrayElectrodeToChannelName(arrayElectrodeStr)
            [array, electrode] = ContinuousNeuralChannelDescriptor.parseArrayElectrode(arrayElectrodeStr);
            assert(~isnan(electrode), 'Could not parse array/electrode string %s', arrayElectrodeStr);
            chName = sprintf('continuous_%s%02d', array, electrode);
        end

        function [array, electrode] = parseArrayElectrode(str)
            tokens = regexp(str, '(continuous_)?(?<array>[A-Za-z_]*)(?<electrode>\d+)', 'names', 'once');
            if isempty(tokens)
                array = '';
                electrode = NaN;
            else
                array = tokens.array;
                electrode =str2double(tokens.electrode);
            end
        end
    end

end