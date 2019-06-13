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
            partitionMetaStruct = p.Results.partitionMeta;
            partitionNames = fieldnames(partitionFieldLists);
            keepfields = @(s, flds) rmfield(s, setdiff(fieldnames(s), flds));
            
            assert(isvector(S));
            N = numel(S);
            
            % create the directory as path/name
            fullPath = GetFullPath(char(locationName));
            
            if exist(fullPath, 'dir')
                TrialDataUtilities.Data.SaveArrayIndividualized.clearLocationContents(fullPath);
            else
                mkdirRecursive(fullPath);
            end
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Saving to %s/', fullPath);
            end
            prog = ProgressBar(N, str);
            for i = 1:N 
                element = S(i);
                
                for iP = 1:numel(partitionNames)
                    % strip off this partition's data and save it into a separate file
                    flds = partitionFieldLists.(partitionNames{iP});
                    
                    % split the partition data from the struct array
                    partData = keepfields(element, flds); %#ok<NASGU>
                    element = rmfield(element, flds);
                    
                    if isempty(callbackFn)
                        file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationName, N, i, partitionNames{iP});
                        save(file, '-v6', 'partData'); % assumes less than 2 GB per element, but much faster
                    else
                        % pass this element, the location, the id number, and the partition to the callback
                        callbackFn(partName, fullPath, i, partitionNames{iP});
                    end
                end
                
                % save this element
                if isempty(callbackFn)
                    file = TrialDataUtilities.Data.SaveArrayIndividualized.generateElementFileName(locationName, N, i);   
                    save(file, '-v6', 'element'); % assumes less than 2 GB per element, but much faster
                else
                    % pass this element, the location, the id number, and indicate not being a partition to the callback
                    callbackFn(element, fullPath, i, '');
                end
                prog.update(i);
            end
            prog.finish();
            
            % save partition meta to separate files
            for iP = 1:numel(partitionNames)
                partitionMeta = partitionMetaStruct.(partitionNames{iP}); %#ok<NASGU>
                file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionMetaFileName(locationName, partitionNames{iP});
                save(file, '-v6', 'partitionMeta'); % assumes less than 2 GB per element, but much faster
            end
            
            % list partitions in partitions.txt
            TrialDataUtilities.Data.SaveArrayIndividualized.writePartitionList(locationName, partitionNames);
                
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
        
        function assertValidLocation(locationName)
            if ~TrialDataUtilities.Data.SaveArrayIndividualized.isValidLocation(locationName)
                error('%s is not a valid location saved with SaveArrayIndividualized', locationName);
            end
        end
        
        function [S, partitionMeta] = loadArray(locationName, varargin)
            p = inputParser();
            p.addParameter('message', '', @ischar);
            p.addParameter('maxElements', Inf, @isscalar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('partitions', {}, @(x) ischar(x) || iscellstr(x));
            p.addParameter('loadAllPartitions', false, @islogical);
            p.addParameter('ignoreMissingPartitions', false, @islogical);
            
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
                    if ~p.Results.ignoreMissingPartitions
                        error('Partitions %s not found. Partitions found: %s', strjoin(partitions(~found), ', '), strjoin(partitionsAvailable, ', '));
                    else
                        partitions = partitions(found);
                    end
                end
            end
            
            % load partition meta
            partitionMeta = struct();
            for iP = 1:numel(partitions)
                file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionMetaFileName(locationName, partitions{iP});
                loaded = load(file); % assumes less than 2 GB per element, but much faster
                partitionMeta.(partitions{iP}) = loaded.partitionMeta;
            end
            
            N = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationName);
            N = min(N, p.Results.maxElements);
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Loading from %s/', locationName);
            end
            prog = ProgressBar(N, str);
            for i = N:-1:1 % reverse order to preallocate array
                if isempty(callbackFn)
                    file = TrialDataUtilities.Data.SaveArrayIndividualized.generateElementFileName(locationName, N, i);   
                    loaded = load(file, 'element');
                    element = loaded.element;
                else
                    element = callbackFn(locationName, i, '');
                end
                
                % load in any partition parts
                for iP = 1:numel(partitions)
                    if isempty(callbackFn)
                        file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationName, N, i, partitions{iP});
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
            countFname = TrialDataUtilities.Data.SaveArrayIndividualized.generateCountFileName(locationName);
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
            
            TrialDataUtilities.Data.SaveArrayIndividualized.assertValidLocation(locationName);
                
            % get partition list as lines of partitions.txt
            partitionFname = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionListFileName(locationName);
            if ~exist(partitionFname, 'file')
                list = {};
                return;
            end
            
            partitionFid = fopen(partitionFname, 'r');
            if partitionFid == -1
                error('Could not open partition file %s', partitionFname);
            end
            tokens = textscan(partitionFid, '%s');
            fclose(partitionFid);
            if isempty(tokens)
                error('Could not read partition file %s', partitionFname);
            end
            list = tokens{1};
        end
        
        function linkPartitionFromOtherLocation(locationNameRef, locationNameSave, varargin)
            % symlinkPartitionFromOtherLocation(locationNameRef, locationNameSave, 'partitions', {'partitionName'}, 'linkAllPartitions', [true/false])
            %
            % symbolically links a partition saved in locationNameRef into the location in locationNameSave

            p = inputParser();
            p.addParameter('mode', 'symlink', @ischar); % 'symlink' or 'copy'
            p.addParameter('overwrite', false, @islogical);
            p.addParameter('partitions', {}, @(x) ischar(x) || iscellstr(x));
            p.addParameter('linkAllPartitions', false, @islogical);
            p.parse(varargin{:});
            
            mode = p.Results.mode;
            switch mode
                case 'symlink'
                    linkFn = @(src, dest) TrialDataUtilities.Save.symlink(src, dest);
                case 'copy'
                    linkFn = @(src, dest) copyfile(src, dest);
                otherwise
                    error('Unknown mode %s', mode);
            end
            
            % check element counts match
            nSave = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationNameSave);
            nRef = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationNameRef);
            assert(nSave == nRef, 'Save location has %d elements but referenced location has %d elements', nSave, nRef);
            N = nSave;
            
            % check partitions
            partitionsSave = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(locationNameSave);
            partitionsRef = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(locationNameRef);
            
            if p.Results.linkAllPartitions
                partitionsLink = partitionsRef;
            else
                partitionsLink = p.Results.partitions;
                
                tf = ismember(partitionsLink, partitionsRef);
                if ~all(tf)
                    error('Partitions %s not found in reference location', strjoin(partitionsLink(~tf), ', '));
                end
            end
            
            if ~p.Results.overwrite
                partitionsLink = setdiff(partitionsLink, partitionsSave);
            end
            
            if isempty(partitionsLink)
                warning('No partitions found in ref that were not already in destination');
                return;
            end
            
            % symlink each meta file
            for iP = 1:numel(partitionsLink)
                fnameRef = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionMetaFileName(locationNameRef, partitionsLink{iP});
                fnameSave = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionMetaFileName(locationNameSave, partitionsLink{iP});
                assert(exist(fnameRef, 'file') == 2, 'Source partition meta file %s not found', fnameRef);
                linkFn(fnameRef, fnameSave);
            end
            
            % do the symlinking for each trial
            prog = ProgressBar(N, 'Symlinking partition files by trial');
            for i = 1:N
                prog.update(i);
                for iP = 1:numel(partitionsLink)
                    fnameRef =  TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationNameRef, N, i, partitionsLink{iP});
                    fnameSave =  TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationNameSave, N, i, partitionsLink{iP});
                    
                    assert(exist(fnameRef, 'file') == 2, 'Source partition file %s not found', fnameRef);
                    linkFn(fnameRef, fnameSave);
                end
            end
            prog.finish();
            
            % update the partitions list
            partitionsNew = union(partitionsSave, partitionsLink);
            TrialDataUtilities.Data.SaveArrayIndividualized.writePartitionList(locationNameSave, partitionsNew);
        end
        
        function clearLocationContents(locationName)
            if exist(locationName, 'dir')
                deleteIfPresent(fullfile(locationName, 'el*.mat'));
                deleteIfPresent(TrialDataUtilities.Data.SaveArrayIndividualized.generateCountFileName(locationName));
                deleteIfPresent(TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionListFileName(locationName));
            end
            
            function deleteIfPresent(file)
                if exist(file, 'file')
                    delete(file);
                end
            end
        end
                
    end
    
    methods(Static, Hidden)
        function fname = generateCountFileName(locationName)
            fname = fullfile(locationName, 'count.txt');
        end
        
        function nzeros = numZerosInFileForCount(count)
            if count >= 10^6
                nzeros = ceil(log10(count + 1));
            else
                nzeros = 6;
            end
        end
        
        function fname = generatePartitionListFileName(locationName)
            fname = fullfile(locationName, 'partitions.txt');
        end
        
        function writePartitionList(locationName, partitions)
            % write list partitions as lines in partitions.txt
            fname = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionListFileName(locationName);
            if isempty(partitions)
                if exist(fname, 'file')
                    delete(fname); % delete any existing partitions file
                end
            else
                partitionFid = fopen(fname, 'w');
                fprintf(partitionFid, '%s\n', partitions{:});
                fclose(partitionFid);
            end
        end
        
        function fname = generatePartitionElementFileName(locationName, N, elementIndex, partition)
            nzeros = TrialDataUtilities.Data.SaveArrayIndividualized.numZerosInFileForCount(N);
            fname = fullfile(locationName, sprintf('el%0*d_partition_%s.mat', nzeros, elementIndex, partition));
        end
        
        function fname = generateElementFileName(locationName, N, elementIndex)
            nzeros = TrialDataUtilities.Data.SaveArrayIndividualized.numZerosInFileForCount(N);
            fname = fullfile(locationName, sprintf('el%0*d.mat', nzeros, elementIndex));
        end
        
        function fname = generatePartitionMetaFileName(locationName, partitionName)
            fname = fullfile(locationName, sprintf('partitionMeta_%s.mat', partitionName));
        end
    end
end
