classdef SaveArrayIndividualized < handle

    methods(Static)            
        function saveArray(locationName, S, varargin)
        % saveStructArrayIndividualized(fname, S, name, varargin)
        % save a struct array S(:) as name1, name2, name3, ...
        % allowing individual elements to be loaded quickly

            p = inputParser();
            p.addParameter('message', '', @ischar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('partitionFieldLists', struct(), @isstruct);
            p.addParameter('partitionMeta', struct(), @isstruct);
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            partitionFieldLists = p.Results.partitionFieldLists;
            partitionMeta = p.Results.partitionMeta;
            partitionNames = fieldnames(partitionFieldLists);
            keepfields = @(s, flds) rmfield(s, setdiff(fieldnames(s), flds));
            
            assert(isvector(S));
            N = numel(S);
            
            % create the directory as path/name
            fullPath = GetFullPath(locationName);
            mkdirRecursive(fullPath);
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Saving to %s/', fullPath);
            end
            prog = ProgressBar(N, str);
            for i = 1:N 
                element = S(i);
                
                for p = 1:numel(partitionNames)
                    % strip off this partition's data and save it into a separate file
                    flds = partitionFieldLists.(partitionNames{p});
                    
                    % split the partition data from the struct array
                    partData = keepfields(element, flds); %#ok<NASGU>
                    element = rmfield(element, flds);
                    
                    if isempty(callbackFn)
                        fname = sprintf('el%06d_partition_%s.mat', i, partitionNames{p});
                        file = fullfile(fullPath, fname);
                        save(file, '-v6', 'partData'); % assumes less than 2 GB per element, but much faster
                    else
                        % pass this element, the location, the id number, and the partition to the callback
                        callbackFn(partName, fullPath, i, partitionNames{p});
                    end
                end
                
                % save this element
                if isempty(callbackFn)
                    fname = sprintf('el%06d.mat', i);
                    file = fullfile(fullPath, fname);
                    save(file, '-v6', 'element'); % assumes less than 2 GB per element, but much faster
                else
                    % pass this element, the location, the id number, and indicate not being a partition to the callback
                    callbackFn(element, fullPath, i, '');
                end
                prog.update(i);
            end
            prog.finish();
            
            % save partition meta to separate files
            for p = 1:numel(partitionNames)
                partitionMeta = partitionMeta.(partitionNames{p});

                fname = sprintf('partitionMeta_%s.mat', partitionNames{p});
                file = fullfile(fullPath, fname);
                save(file, '-v6', 'partitionMeta'); % assumes less than 2 GB per element, but much faster
            end
            
            if ~isempty(partitionNames)
                % list partitions in partitions.txt
                partitionFid = fopen(fullfile(locationName, 'partitions.txt'), 'w');
                fprintf(partitionFid, '%s\n', partitionNames{:});
                fclose(partitionFid);
            end
                
            % create count file containing just N
            countFid = fopen(fullfile(locationName, 'count.txt'), 'w');
            fprintf(countFid, '%d\n', N); 
            fclose(countFid);
        end
        
        function tf = isValidLocation(locationName)
            locationName = GetFullPath(locationName);
            if ~exist(locationName, 'dir')
                tf = false; return;
            end
            if ~exist(fullfile(locationName, 'count.txt'), 'file')
                tf = false; return;
            end
            fname = sprintf('el%06d.mat', 1);
            if ~exist(fullfile(locationName, fname), 'file')
                tf = false; return;
            end
            tf = true;
        end
        
        function [S, partitionMeta] = loadArray(locationName, varargin)
            p = inputParser();
            p.addParameter('message', '', @ischar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('partitions', {}, @(x) ischar(x) || iscellstr(x));
            p.addParameter('loadAllPartitions', false, @islogical);
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            
            locationName = GetFullPath(locationName);
            
            function s = structMerge(s, s2)
                flds = fieldnames(s2);
                for f = 1:numel(flds)
                    s.(flds{f}) = s2.(flds{f});
                end
            end
            
            % check that partitions are found
            partitionsAvailable = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(locationName);
            if p.Results.loadAllPartitions
                partitions = partitionsAvailable;
            else
                partitions = p.Results.partitions;
                if ischar(partitions)
                    partitions = {partitions};
                end
                found = ismember(partitions, partitionsAvailable);
                if any(~found)
                    error('Partitions %s not found. Partitions found: %s', strjoin(partitions(~found), ', '), strjoin(partitionsAvailable, ', '));
                end
            end
            
            % load partition meta
            partitionMeta = struct();
            for iP = 1:numel(partitions)
                fname = sprintf('partitionMeta_%s.mat', partitions{iP});
                file = fullfile(locationName, fname);
                loaded = load(file); % assumes less than 2 GB per element, but much faster
                partitionMeta.(partitions{iP}) = loaded.partitionMeta;
            end
            
            N = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationName);
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Loading from %s/', locationName);
            end
            prog = ProgressBar(N, str);
            for i = N:-1:1 % reverse order to preallocate array
                if isempty(callbackFn)
                    fname = sprintf('el%06d.mat', i);
                    file = fullfile(locationName, fname);
                    loaded = load(file, 'element');
                    element = loaded.element;
                else
                    element = callbackFn(locationName, i, '');
                end
                
                % load in any partition parts
                for iP = 1:numel(partitions)
                    if isempty(callbackFn)
                        fname = sprintf('el%06d_partition_%s.mat', i, partitions{iP});
                        file = fullfile(locationName, fname);
                        loaded = load(file); % assumes less than 2 GB per element, but much faster
                        partData = loaded.partData;
                    else
                        % pass this element, the location, the id number, and the partition to the callback
                        partData = callbackFn(locationName, i, partitions{iP});
                    end
                    element = structMerge(element, partData);
                end
                
                S(i) = element;
                prog.update(N-i+1);
            end
            prog.finish();
        end
        
        function N = getArrayCount(locationName)
            locationName = GetFullPath(locationName);
            
            % get element count from count.txt file
            countFname = fullfile(locationName, 'count.txt');
            countFid = fopen(countFname, 'r');
            if countFid == -1
                error('Could not find count file %s', countFname);
            end
            tokens = textscan(countFid, '%d', 1);
            fclose(countFid);
            N = double(tokens{1});
            if isempty(N)
                error('Could not read count file %s', countFname);
            end
        end
        
        function list = listPartitions(locationName)
            locationName = GetFullPath(locationName);
            
            % get element count from count.txt file
            partitionFname = fullfile(locationName, 'partitions.txt');
            partitionFid = fopen(partitionFname, 'r');
            if partitionFid == -1
                error('Could not find partition file %s', partitionFname);
            end
            tokens = textscan(partitionFid, '%s');
            fclose(partitionFid);
            if isempty(tokens)
                error('Could not read partition file %s', partitionFname);
            end
            list = tokens{1};
        end
    end
end
