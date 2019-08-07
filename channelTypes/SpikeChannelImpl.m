classdef SpikeChannelImpl < ChannelImpl
    methods
        function impl = SpikeChannelImpl(cd)
            impl.cd = cd;
        end
        
        function data = convertDataCellOnAccess(impl, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (impl.scaleFromLims -> impl.scaleToLims)
            cd = impl.cd;
            data = convertDataCellOnAccess@ChannelImpl(impl, fieldIdx, data);
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelDescriptor.scaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
        end

        function data = convertDataSingleOnAccess(impl, fieldIdx, data)
            cd = impl.cd;
            data = convertDataSingleOnAccess@ChannelImpl(impl, fieldIdx, data);
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelImpl.scaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
        end

        function data = convertAccessDataCellToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelImpl.unscaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
            data = convertAccessDataCellToMemory@ChannelImpl(impl, fieldIdx, data);
        end

        function data = convertAccessDataSingleToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            if cd.hasWaveforms && fieldIdx == 2
                data = ChannelImpl.unscaleData(data, cd.waveformsScaleFromLims, cd.waveformsScaleToLims);
            end
            data = convertAccessDataSingleToMemory@ChannelImpl(impl, fieldIdx, data);
        end

        function waveData = scaleWaveforms(impl, waveData)
            if iscell(waveData)
                waveData = impl.convertDataCellOnAccess(2, waveData);
            else
                waveData = impl.convertDataSingleOnAccess(2, waveData);
            end
        end

        function waveData = unscaleWaveforms(impl, waveData)
            if iscell(waveData)
                waveData = impl.convertAccessDataCellToMemory(2, waveData);
            else
                waveData = impl.convertAccessDataSingleToMemory(2, waveData);
            end
        end
    end
end