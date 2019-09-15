classdef AnalogChannelImpl < ChannelImpl
    methods
        function impl = AnalogChannelImpl(cd)
            impl.cd = cd;
        end
        
        function data = convertDataCellOnAccess(impl, fieldIdx, data)
            % cast to access class, also do scaling upon request
            cd = impl.cd;
            data = convertDataCellOnAccess@ChannelImpl(impl, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
            end
        end
        
        function data = convertDataSingleOnAccess(impl, fieldIdx, data)
            cd = impl.cd;
            data = convertDataSingleOnAccess@ChannelImpl(impl, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
            end
        end
        
        function data = convertAccessDataCellToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
            end
            data = convertAccessDataCellToMemory@ChannelImpl(cd, fieldIdx, data);
        end
        
        function data = convertAccessDataSingleToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
            end
            data = convertAccessDataSingleToMemory@ChannelImpl(cd, fieldIdx, data);
        end
        
        function [dataCell, timeCell] = computeTransformDataRaw(impl, td, varargin)
            % if shared column, use the slice arg to subselect the
            % appropriate column
            cd = impl.cd;
            if cd.isColumnOfSharedMatrix
                sliceArgs = {cd.primaryDataFieldColumnIndex};
            else
                sliceArgs = {};
            end
            
            [dataCell, timeCell] = AnalogChannelGroupImpl.doComputeTransformData(impl, td, 'slice', sliceArgs, varargin{:});
        end
    end
end