classdef ChannelDescriptor < matlab.mixin.Heterogeneous
    % Manually specified meta data
    properties
        name = ''; % short name, must be valid field name

        groupName = ''; % name of group to which this channel belongs 

        description = ''; % extended description 

        units = ''; % string describing units

        meta % anything you'd like

        special = false; % whether this channel is a "special" identifier channel used by TrialData

        datenum = false; % if true, treat as datenum value 

        %Fs = NaN; % sampling frequency in Hz 
    end

    % Inferred from data by inferAttributesFromData below
    properties
        storageDataClass = 'double'; % original class name, for storage purposes

        dataClass % class used in memory
        
        numeric % the data are numeric

        logical % the data are logical

        string % the data are strings

        vector % all values are vector
        
        scalar % all values are scalar 

        missingValue = []; % value inserted when the original data was missing

        contentsAreCell % each trial has a cell array inside it
        
        collectAsCell = true; % use a cell to aggregate the data rather than a numeric matrix

        ndims % number of dimensions occupied by this signal

        sizeConsistent = false; % is the size the same in every trial?
        size  % size of the signal along each dimension or NaN if variable

        concatenateAlongDim % when dimension along which trials would be concatenated, or NaN if this won't work
    end

    methods(Abstract)
        % return a string with this channels type
        type = getType(cdesc);

        str = describe(cdesc);

        dataFields = getDataFields(cdesc);
    end

    methods
        function cd = ChannelDescriptor(varargin)
            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.parse(varargin{:});

            cd.name = p.Results.name;
        end

        function str = getAxisLabel(cd)
            if isempty(cd.units)
                str = sprintf('%s', cd.name);
            else
                str = sprintf('%s (%s)', cd.name, cd.units);
            end
        end

        function cd = inferAttributesFromData(cd, dataCell)
            assert(nargout > 0, 'ChannelDescriptor is not a handle class. If the return value is not stored this call has no effect');

            % infer data class
            if ~iscell(dataCell)
                if ischar(dataCell)
                    dataCell = {dataCell};
                else
                    dataCell = num2cell(dataCell);
                end
            end
            
            emptyMask = cellfun(@isempty, dataCell);
            if all(emptyMask)
                warning('Data cell contains no non-empty values for channel %s, assuming class double', cd.name);
                cd.storageDataClass = 'double';
            else
                % determine original class in R struct
                classes = unique(cellfun(@class, dataCell(~emptyMask), 'UniformOutput', false));
                if numel(classes) > 1
                    error('Data cell struct contains multiple classes for %s : %s\n', cd.name, strjoin(classes, ', ')); 
                end
                cd.storageDataClass = classes{1};
            end

            % determine if values are strings
            cd.string = all(cellfun(@(x) ischar(x), dataCell(~emptyMask)));

            % determine if all values are scalar or vector, logical, numeric, or neither
            cd.contentsAreCell = any(cellfun(@iscell, dataCell(~emptyMask)));
            cd.scalar = ~cd.string && all(cellfun(@(x) isscalar(x), dataCell(~emptyMask)));
            cd.vector = ~cd.string && all(cellfun(@(x) isvector(x), dataCell(~emptyMask))); 
            cd.numeric = ismember(cd.storageDataClass, {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'});
            cd.logical = strcmp(cd.storageDataClass, 'logical');

            % compute the size and size consistency of the signal
            ndimsEach = cellfun(@ndims, dataCell(~emptyMask));
            cd.ndims = max(ndimsEach);

            sizeEach = cell(cd.ndims, 1);
            [sizeEach{:}] = cellfun(@size, dataCell(~emptyMask));
            sizeEach = cat(1, sizeEach{:})';
            sizeUnique = unique(sizeEach, 'rows');
            if size(sizeUnique, 1) == 1
                cd.sizeConsistent = true;
                cd.size = sizeUnique;
            else
                cd.size = nanvec(cd.ndims);
                for iDim = 1:cd.ndims
                    if numel(unique(sizeUnique(:, iDim))) == 1
                        cd.size(iDim) = sizeUnique(1, iDim);
                    end
                end
                        
                cd.sizeConsistent = false;
            end

            % aggregate in matrix or cell?
            cd.collectAsCell = ~cd.scalar;

            % compute in memory data type
            if cd.string
                cd.dataClass = 'char';
            elseif cd.numeric
                cd.dataClass = 'double';
            elseif cd.logical
                cd.dataClass = 'logical';
            end
            
            % compute concatenation dimension
            if cd.vector
                cd.concatenateAlongDim = 1;
            elseif cd.sizeConsistent
                cd.concatenateAlongDim = cd.ndims + 1; % all dims are consistent along a new dimension
            elseif~any(isnan(cd.size(1:end-1)))
                % all dims are consistent except last, concatenate along
                % that one
                cd.concatenateAlongDim = cd.ndims;
            else
                cd.concatenateAlongDim = NaN;
            end
            
            % compute missing value to fill empty trials with
            if cd.string
                cd.missingValue = '';
            elseif cd.scalar
                if cd.logical
                    cd.missingValue = false;
                else
                    cd.missingValue = NaN;
                end
            else
                cd.missingValue = [];
            end
        end

        % look at the data struct and check that the correct data fields
        % are present in the R struct
        function [ok, msg] = checkData(cd, R)
            fields = cd.getDataFields();
            if ~all(isfield(R, fields))
                ok = false;
                msg = sprintf('Missing fields %s', strjoin(fields(~isfield(R, fields)), ', '));
            else
                ok = true;
                msg = '';
            end
        end

        % do any replacement of missing values, etc.
        function R = repairData(cd, R)
            % replace empty values with default in all data fields
            fields = cd.getDataFields();
            for i = 1:length(fields)
                fld = fields{i};
                if ~isempty(cd.missingValue)
                    emptyMask = cellfun(@isempty, {R.(fld)});
                    [R(emptyMask).(fld)] = deal(cd.missingValue);
                end
            end
        end

        function data = convertData(cd, data)
            fields = cd.getDataFields();
            data = structConvertFieldValues(data, 'double', fields);
        end 
    end
end
