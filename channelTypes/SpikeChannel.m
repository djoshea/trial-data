classdef SpikeChannelDescriptor < EventChannelDescriptor
    properties
        unit
        channel 
    end

    methods
        function type = getType(cdesc)
            type = 'spike';
        end

        function str = describe(cdesc)
            str = sprintf('Spike %d.%d', cdesc.channel, cdesc.unit);  
        end

        function dataFields = getExtraDataFields(cdesc)
            dataFields = {'waveforms'};
        end

        function cd = SpikeChannelDescriptor(varargin)
            cd = cd@EventChannelDescriptor(varargin{:});
        end
    end

end
