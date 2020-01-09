classdef ChannelImpl 
    properties
        cd
    end
    
    methods
        function impl = ChannelImpl(cd)
            if nargin >= 1
                impl.cd = cd;
            end
        end
        
        function data = convertDataToCorrectVectorFormat(impl, fieldIdx, data)
            % if collectAsCellByField(fieldIdx) is true, convert to cell
            % if false, convert to numeric vector with appropriate missing
            % values inserted to equalize length
            
            cd = impl.cd; %#ok<*PROPLC>
            fieldIdx = cd.lookupFieldId(fieldIdx);

            if cd.collectAsCellByField(fieldIdx)
                if ~iscell(data)
                    data = num2cell(data);
                end
            else
                % convert to numeric vector, check for non-scalar values,
                % and fill missing values
                if iscell(data)
                    if ismember(cd.elementTypeByField(fieldIdx), [cd.VECTOR, cd.NUMERIC])
                        % numeric data cat along first dim
                        nonEmpty = ~cellfun(@isempty, data);
                        if cd.elementTypeByField(fieldIdx) == cd.VECTOR
                            data = cellfun(@makerow, data, 'UniformOutput', false);
                        end
                        try
                            mat = cell2mat(data(nonEmpty));
                        catch
                            error('Numeric data that is to be concatenated along first dim has uneven sizing');
                        end

                        data = TensorUtils.inflateMaskedTensor(mat, 1, nonEmpty, cd.missingValueByField{fieldIdx});
                    elseif ismember(cd.elementTypeByField(fieldIdx), cd.STRING)
                        if iscell(data)
                            % TODO address this
                            emptyMask = cellfun(@(x) isempty(x) || (isscalar(x) && ismissing(x)), data);
                            data(emptyMask) = {''};
                        end
                        assert(iscellstr(data) || isstring(data), 'String data field must be cellstr or string');
                        data = string(data);
                        assert(isempty(data) || isvector(data));

                        % replace empty values with missing?
                        emptyMask = arrayfun(@(x) x == "", data);
                        data(emptyMask) = string(missing());
                    elseif ismember(cd.elementTypeByField(fieldIdx), cd.CATEGORICAL)
                        if iscell(data)
                            emptyMask = cellfun(@(x) isempty(x) || (isscalar(x) && ismissing(x)), data);
                            data(emptyMask) = {categorical(missing)};
                        end
                        charMask = cellfun(@ischar, data);
                        data(charMask) = cellfun(@(x) categorical(string(x)), data(charMask), 'UniformOutput', false);
                        
                        data = cellfun(@categorical, data);
                        
                    else
                        missingVal = cd.missingValueByField{fieldIdx};
                        assert(isscalar(missingVal));
                        nVals = cellfun(@numel, data);
                        if any(nVals > 1)
                            throwError('Data must contain scalar values for each trial');
                        end
                        [data{nVals==0}] = deal(missingVal);
                        nel = numel(data);
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, class(missingVal));
                        data = ChannelImpl.cellCast(data, newClass);
                        data = cat(1, data{:});
                        assert(isempty(data) || isvector(data) && numel(data) == nel);
                    end
                end
            end

            function throwError(varargin)
                error(['Error in channel %s, field %d: ' varargin{1}], cd.name, fieldIdx, varargin{2:end});
            end
        end

        function [cd, data] = checkConvertDataAndUpdateMemoryClassToMakeCompatible(impl, fieldIdx, data)
            % this function does a couple of things. It takes a cell or
            % vector of data destined for a specific field of this channel
            % descriptor. It first checks whether this data is at all
            % acceptable for this field type (i.e. BOOLEAN, SCALAR, VECTOR,
            % above). It will throw an error if not.
            %
            % If the data are compatible, it will then check the class of
            % data against my .memoryClassByField{fieldIdx}. If it is
            % possible to convert data to memClass, data will be converted.
            % If not, then memClass will be updated to reflect the change.
            % For example, if memClass is single and data(i) = uint16(1),
            % data(i) will be converted to single(1). However, if memClass
            % is uint16 and data(i) = double(1.5), then memClass will be changed to
            % double. This will not change any other data set on this field
            % in the TrialData instance, but this preexisting data will be
            % cast into the access class on access anyway, so it's not
            % necessary to worry about it now.if not we
            % change the class to match the new format.

            % if meant to be collected as a vector, do that first as it
            % simplifies checking the types
            cd = impl.cd;
            fieldIdx = cd.lookupFieldId(fieldIdx);
            iF = fieldIdx;

            % convert to cell or vector depending on collectAsCellByField
            data = impl.convertDataToCorrectVectorFormat(iF, data);

            memClass = cd.memoryClassByField{iF};
            switch cd.elementTypeByField(iF)
                case cd.BOOLEAN
                    data(isnan(data)) = false;
                    convertedData = logical(data);
                    if any(convertedData ~= data)
                        throwError('Data must be logical or convertible to logical vector');
                    end
                    newClass = memClass;

                case {cd.SCALAR, cd.DATENUM}
                    newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
                    data = ChannelImpl.icast(data, newClass);

                case cd.VECTOR
                    if cd.collectAsCellByField(iF) && iscell(data)
                        okay = cellfun(@(x) isempty(x) || isvector(x), data);
                        if ~all(okay)
                            throwError('Data cell contents must be vectors or empty');
                        end
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);
                        data = ChannelImpl.cellCast(data, newClass);
                    else
                        newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
%                         if strcmp(newClass, 'logical')
%                             data(isnan(data)) = false;
%                         end
                        data = ChannelImpl.icast(data, newClass);
                    end

                    if cd.collectAsCellByField(iF) && iscell(data)
                        % we'll be cat'ing them along dim 1.
                        data = cellfun(@makecol, data, 'UniformOutput', false);
                    end

                case cd.NUMERIC
                    if cd.collectAsCellByField(iF) && iscell(data)
                        okay = cellfun(@(x) isempty(x) || (isnumeric(x) || islogical(x)), data);
                        if ~all(okay)
                            throwError('Data cell contents must be numeric');
                        end
                        newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);
                        data = ChannelImpl.cellCast(data, newClass);
                    else
                        newClass = TrialDataUtilities.Data.determineCommonClass(class(data), memClass);
%                         if strcmp(newClass, 'logical')
%                             data(isnan(data)) = false;
%                         end
                        data = ChannelImpl.icast(data, newClass);
                    end

                case cd.STRING
                    okay = cellfun(@(x) isempty(x) || (ischar(x) && isvector(x)), data);
                    if ~all(okay)
                        throwError('Data cell contents must be strings');
                    end
                    data = cellfun(@makerow, data, 'UniformOutput', false);
                    newClass = 'char';

                case cd.CELL
                    okay = cellfun(@(x) isempty(x) || iscell(x), data);
                    if ~all(okay)
                        throwError('Data cell contents must be cells');
                    end
                    % this used to be 'cell', but it wasn't correct for spike
                    % array waveforms. might need to be fixed if this isn't
                    % sufficient
                    newClass = TrialDataUtilities.Data.cellDetermineCommonClass(data, memClass);
                    
                case cd.CATEGORICAL
                    newClass = 'categorical';
                    data = ChannelImpl.icast(data, newClass);

                otherwise
                    throwError('Unknown element type')
            end

            cd.originalDataClassByField{iF} = newClass;

            function throwError(varargin)
                error(['Error in channel %s, field %d: ' varargin{1}], cd.name, fieldIdx, varargin{2:end});
            end
        end

        function data = convertDataSingleOnAccess(impl, fieldIdx, data)
            cd = impl.cd;
            fieldIdx = cd.lookupFieldId(fieldIdx);
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = ChannelImpl.icast(data, accClass);
            end
        end

        function data = convertDataCellOnAccess(impl, fieldIdx, data)
            cd = impl.cd;
            fieldIdx = cd.lookupFieldId(fieldIdx);
            
            % when data is accessed by TrialData, this is one additional
            % chance to cast or adjust the data stored in .data. This
            % default implementation casts from the memory class to the
            % access class on demand
            %memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};

            data = ChannelImpl.cellCast(data, accClass);
            % data = cellfun(@(a) cast(a, accClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            if ~cd.collectAsCellByField(fieldIdx)
                % for vector types, we can makerow the contents to ensure
                % that cell2mat works as intended
                if cd.isVectorByField(fieldIdx)
                    data = cellfun(@makerow, data, 'UniformOutput', false);
                end

                % because of the way missingValue works, empty trials will
                % be NaN for invalid trials when setting channel data
                % trials missing values will have a single Nan
                nanMask = cellfun(@(x) isscalar(x) && ismissing(x), data);
                if any(~nanMask)
                    data = cat(1, data{~nanMask});
                else
                    data = data{1}([], :, :, :, :, :, :);
%                     data = zeros(0, 1);
                end
                data = TensorUtils.inflateMaskedTensor(data, 1, ~nanMask);
            end
        end

        function data = convertAccessDataSingleToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            fieldIdx = cd.lookupFieldId(fieldIdx);
            
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = ChannelImpl.icast(data, memClass);
            end
        end

        function data = convertAccessDataCellToMemory(impl, fieldIdx, data)
            cd = impl.cd;
            fieldIdx = cd.lookupFieldId(fieldIdx);
            
            % when data is accessed by TrialData, this is one additional
            % chance to cast or adjust the data stored in .data. This
            % default implementation casts from the access class to the
            % memory class on demand
            memClass = cd.memoryClassByField{fieldIdx};
            accClass = cd.accessClassByField{fieldIdx};
            if ~strcmp(memClass, accClass)
                data = cellfun(@(a) ChannelImpl.icast(a, memClass), data, 'UniformOutput', ~cd.collectAsCellByField(fieldIdx));
            end
        end
        
        function [ok, missing] = checkData(impl, data)
            % look at the data struct and check that the correct data fields
            % are present in the R struct
        
            cd = impl.cd;
            fields = cd.dataFields;
            if ~all(isfield(data, fields))
                ok = false;
                missing = fields(~isfield(data, fields));
            else
                ok = true;
                missing = {};
            end
        end

        function data = addMissingFields(impl, data)
            % manually add the fields needed by this channel to data
            % struct, filling them with the appropriate missing value
            
            cd = impl.cd;
            for iF = 1:cd.nFields
                fld = cd.dataFields{iF};
                if ~isfield(data, fld)
                    % insert field with appropriate missing values
                    [data(:).(fld)] = deal(cd.missingValueByField{iF});
                end
            end
        end

        function data = correctMissingValueInData(impl, iF, data)
            cd = impl.cd;
            missingValue = cd.missingValueByField{iF};

            if iscell(data)
                replace = cellfun(@(x) isempty(x) || (isscalar(x) && ismissing(x)), values);
                [data(replace)] = deal(missingValue);
            else
                replace = isnan(values);
                data(replace) = missingValue;
            end
        end

    end
    
    methods
        function str = getAxisLabelPrimary(impl)
            cd = impl.cd; %#ok<*PROP>
            if isempty(cd.unitsByField{1})
                str = sprintf('%s', cd.name);
            else
                str = sprintf('%s (%s)', cd.name, cd.unitsPrimary);
            end
        end
    end
    
    methods(Static)
        function data = icast(data, newClass)
            if ~isa(data, newClass)
                if strcmp(newClass, 'logical')
                    data(isnan(data)) = false;
                    data = logical(data);
                elseif strcmp(newClass, 'categorical')
                    data = categorical(data);
                elseif strcmp(newClass, 'string')
                    data = string(data);
                else
                    data = cast(data, newClass);
                end
            end
        end

        function data = scaleData(data, scaleFromLims, scaleToLims)
            if isempty(scaleFromLims) || isempty(scaleToLims)
                return;
            end
            scaleFromLow = scaleFromLims(2);
            scaleFromRange = scaleFromLims(2) - scaleFromLims(1);
            scaleToLow = scaleToLims(2);
            scaleToRange = scaleToLims(2) - scaleToLims(1);
            if iscell(data)
                data = cellfun(@(d) (d-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow, data, 'UniformOutput', false);
            else
                data = (data-scaleFromLow)*(scaleToRange/scaleFromRange) + scaleToLow;
            end
        end

        function data = unscaleData(data, scaleFromLims, scaleToLims)
            data = ChannelImpl.scaleData(data, scaleToLims, scaleFromLims);
        end

        function cls = getCellElementClass(dataCell)
            if isempty(dataCell)
                cls = 'double';
            elseif ~iscell(dataCell)
                cls = class(dataCell);
            else
                nonEmpty = ~cellfun(@isempty, dataCell);
                if ~any(nonEmpty)
                    cls = class(dataCell{1}); % doing this instead in case its char
%                     cls = 'double'; % assume double if no values found
                else
                    first = find(nonEmpty, 1);
                    cls = class(dataCell{first});
                end
            end
        end

        function sz = getCellElementSize(dataCell)
            if isempty(dataCell)
                sz = [];
            elseif ~iscell(dataCell)
                sz = size(dataCell);
            else
                nonEmpty = ~cellfun(@isempty, dataCell);
                if ~any(nonEmpty)
                    sz = [];
                else
                    first = find(nonEmpty, 1);
                    sz = size(dataCell{first});
                end
            end
        end

        function data = cellCast(data, newClass)
            for i = 1:numel(data)
                data{i} = ChannelImpl.icast(data{i}, newClass);
            end
        end
    end
end