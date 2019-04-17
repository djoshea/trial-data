classdef AnalogChannelGroupImpl < ChannelImpl
    methods
        function impl = AnalogChannelGroupImpl(cd)
            impl.cd = cd;
        end
        
         function data = convertDataCellOnAccess(cd, fieldIdx, data)
            % cast to access class, also do scaling upon request
            % (cd.scaleFromLims -> cd.scaleToLims)
            data = convertDataCellOnAccess@ChannelImpl(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims) && ~isempty(data)
                data = ChannelImpl.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(1);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(1);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 scaleFn = @(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
%                 if iscell(data)
%                     data = cellfun(scaleFn, data, 'UniformOutput', false);
%                 else
%                     data = scaleFn(data);
%                 end
            end
        end
        
        function data = convertDataSingleOnAccess(cd, fieldIdx, data)
            data = convertDataSingleOnAccess@ChannelImpl(cd, fieldIdx, data);
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = (data-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
            end
        end
        
        function data = convertAccessDataCellToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.unscaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = cellfun(@(d) (d-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow, data, 'UniformOutput', false);
            end
            data = convertAccessDataCellToMemory@ChannelImpl(cd, fieldIdx, data);
        end
        
        function data = convertAccessDataSingleToMemory(cd, fieldIdx, data)
            if fieldIdx == 1 && ~isempty(cd.scaleFromLims) && ~isempty(cd.scaleToLims)
                data = ChannelImpl.scaleData(data, cd.scaleFromLims, cd.scaleToLims);
%                 scaleFromLow = cd.scaleFromLims(2);
%                 scaleFromRange = cd.scaleFromLims(2) - cd.scaleFromLims(1);
%                 scaleToLow = cd.scaleToLims(2);
%                 scaleToRange = cd.scaleToLims(2) - cd.scaleToLims(1);
%                 data = (data-scaleToLow)*(scaleFromRange/scaleToRange) + scaleFromLow;
            end
            data = convertAccessDataSingleToMemory@ChannelImpl(cd, fieldIdx, data);
        end
        
        function [dataCell, timeCell] = computeTransformDataRaw(cd, td, varargin)
            [dataCell, timeCell] = AnalogChannelGroupDescriptor.doComputeTransformData(cd, td, varargin{:});
        end
             
    end
    
     methods(Static)
        function [dataCell, timeCell] = doComputeTransformData(cd, td, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical);
            p.addParameter('slice', {}, @(x) true); % this is used to index specifically into each sample
            p.addParameter('timeCell', [], @(x) true);
            p.addParameter('sort', false, @islogical);
            p.parse(varargin{:});
            
            % get raw data from trial data
            [dataCell, timeCell] = td.getAnalogChannelGroupMulti(cd.transformChannelNames, 'raw', true, ...
                'applyScaling', p.Results.applyScaling, 'sort', p.Results.sort);
            
            % use transform function
            switch cd.transformFnMode
                case 'simple'
                    % simple function, call it on one trial at a time and we'll
                    % handle the slicing
                    args = p.Results.slice;
                    if ~iscell(args)
                        args = {args};
                    end

                    prog = ProgressBar(numel(dataCell), 'Computing transform analog channel on the fly');
                    for iT = 1:numel(dataCell)
                        if ~isempty(dataCell{iT})
                            prog.update(iT);
                            dataCell{iT} = cd.transformFn(dataCell{iT});
                            if ~isempty(args)
                                dataCell{iT} = dataCell{iT}(:, args{:}); % take a slice through the data
                            end
                        end
                    end
                    prog.finish();
                    timeCell = timeCell(:, 1);
                    
                case 'default'
                    % more advanced function, can handle all data at once and slicing
                    if ~isempty(p.Results.slice)
                        args = p.Results.slice;
                        if ~iscell(args)
                            args = {args};
                        end
                    else
                        args = {};
                    end
                    
                   [dataCell, timeCell] = cd.transformFn(dataCell, timeCell, 'slice', args, 'scalingApplied', p.Results.applyScaling);
            end
        end
     end
end