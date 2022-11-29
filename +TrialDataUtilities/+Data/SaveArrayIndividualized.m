classdef SaveArrayIndividualized < handle

    methods(Static)            
        function saveArray(locationName, S, varargin)
        % saveStructArrayIndividualized(fname, S, name, varargin)
        % save a struct array S(:) as name1, name2, name3, ...
        % allowing individual elements to be loaded quickly

            p = inputParser();
            p.addParameter('message', '', @isstringlike);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('partitionFieldLists', struct(), @isstruct);
            p.addParameter('partitionMeta', struct(), @isstruct);
            p.addParameter('partialLoadData', [], @(x) isempty(x) || isstruct(x));
            p.addParameter('elementsPerChunk', 1, @isscalar); % split elements into files with this many elements per file
            p.addParameter('progress', true, @islogical);
            p.addParameter('use_parallel', ~isempty(gcp('nocreate')), @islogical);
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            partitionFieldLists = p.Results.partitionFieldLists;
            partitionMetaStruct = p.Results.partitionMeta;
            partitionNames = fieldnames(partitionFieldLists);
            elementsPerChunk = p.Results.elementsPerChunk;
            progress = p.Results.progress && exist('ProgressBar', 'class');
            use_parallel = p.Results.use_parallel;
            
            assert(isvector(S));
            S = makecol(S);
            N = numel(S);
            
            % create the directory as path/name
            fullPath = TrialDataUtilities.Path.GetFullPath(char(locationName));
            
            if exist(fullPath, 'dir')
                TrialDataUtilities.Data.SaveArrayIndividualized.clearLocationContents(fullPath);
            else
                TrialDataUtilities.Save.mkdirRecursive(fullPath);
            end
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Saving to %s/', fullPath);
            end
            
            whichChunk = TrialDataUtilities.Data.SaveArrayIndividualized.buildElementChunkLists(N, elementsPerChunk);
            nChunks = whichChunk(end);
            if progress
                prog = ProgressBar(nChunks, str);
            end
            
            nByChunk = histcounts(whichChunk, 0.5:nChunks+0.5)';
            S_by_chunk = mat2cell(S, nByChunk);
            if use_parallel
                if progress, prog.enableParallel(); end
                parfor c = 1:nChunks 
                    element = S_by_chunk{c};
                    TrialDataUtilities.Data.SaveArrayIndividualized.internal_saveChunk(element, locationName, nChunks, c, partitionNames, partitionFieldLists, callbackFn);
                    if progress, prog.update(c); end %#ok<PFBNS>
                end
            else 
                for c = 1:nChunks 
                    element = S_by_chunk{c};
                    TrialDataUtilities.Data.SaveArrayIndividualized.internal_saveChunk(element, locationName, nChunks, c, partitionNames, partitionFieldLists, callbackFn);
                    if progress, prog.update(c); end
                end
            end
            if progress, prog.finish(); end
            
%             for c = 1:nChunks 
%                 element = S(whichChunk == c);
% 
%                 for iP = 1:numel(partitionNames)
%                     % strip off this partition's data and save it into a separate file
%                     flds = partitionFieldLists.(partitionNames{iP});
% 
%                     % split the partition data from the struct array
%                     partData = keepfields(element, flds);
%                     element = rmfield(element, flds);
% 
%                     if isempty(callbackFn)
%                         file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationName, nChunks, c, partitionNames{iP});
%                         save(file, '-v6', 'partData'); % assumes less than 2 GB per element, but much faster
%                     else
%                         % pass this element, the location, the id number, and the partition to the callback
%                         callbackFn(partName, fullPath, c, partitionNames{iP});
%                     end
%                 end
% 
%                 % save this element
%                 if isempty(callbackFn)
%                     file = TrialDataUtilities.Data.SaveArrayIndividualized.generateElementFileName(locationName, nChunks, c);   
%                     save(file, '-v6', 'element'); % assumes less than 2 GB per element, but much faster
%                 else
%                     % pass this element, the location, the id number, and indicate not being a partition to the callback
%                     callbackFn(element, fullPath, c, '');
%                 end
%                 if progress, prog.update(c); end
%             end
%             if progress, prog.finish(); end
                
            % save partition meta to separate files
            for iP = 1:numel(partitionNames)
                partitionMeta = partitionMetaStruct.(partitionNames{iP});
                file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionMetaFileName(locationName, partitionNames{iP});
                save(file, '-v6', 'partitionMeta'); % assumes less than 2 GB per element, but much faster
            end
            
            % list partitions in partitions.txt
            TrialDataUtilities.Data.SaveArrayIndividualized.writePartitionList(locationName, partitionNames);
                
            % create count file containing just N
            countFid = fopen(fullfile(locationName, 'count.txt'), 'w');
            fprintf(countFid, '%d %d\n', N, elementsPerChunk); 
            fclose(countFid);
            
            % save partial load data to partialLoad.mat
            partialLoadData = p.Results.partialLoadData;
            if ~isempty(partialLoadData) && ~isempty(fieldnames(partialLoadData))
                fname = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartialLoadDataFileName(locationName);
                save(fname, 'partialLoadData');
            end
        end
        
        function internal_saveChunk(element, locationName, nChunks, chunk_ind, partitionNames, partitionFieldLists, callbackFn)
            keepfields = @(s, flds) rmfield(s, setdiff(fieldnames(s), flds));
            
            for iP = 1:numel(partitionNames)
                % strip off this partition's data and save it into a separate file
                flds = partitionFieldLists.(partitionNames{iP});

                % split the partition data from the struct array
                partData = keepfields(element, flds);
                element = rmfield(element, flds);

                if isempty(callbackFn)
                    file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationName, nChunks, chunk_ind, partitionNames{iP});
                    save(file, '-v6', 'partData'); % assumes less than 2 GB per element, but much faster
                else
                    % pass this element, the location, the id number, and the partition to the callback
                    callbackFn(partName, fullPath, c, partitionNames{iP});
                end
            end

            % save this element
            if isempty(callbackFn)
                file = TrialDataUtilities.Data.SaveArrayIndividualized.generateElementFileName(locationName, nChunks, chunk_ind);   
                save(file, '-v6', 'element'); % assumes less than 2 GB per element, but much faster
            else
                % pass this element, the location, the id number, and indicate not being a partition to the callback
                callbackFn(element, fullPath, chunk_ind, '');
            end
        end
        
        function tf = isValidLocation(locationName)
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
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
        
        function [S, partitionMeta, partialLoadMask] = loadArray(locationName, varargin)
            p = inputParser();
            p.addParameter('progress', true, @islogical);
            p.addParameter('message', '', @isstringlike);
            p.addParameter('maxElements', Inf, @isscalar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('partitions', {}, @(x) TrialDataUtilities.String.isstringlike(x));
            p.addParameter('loadAllPartitions', false, @islogical);
            p.addParameter('ignoreMissingPartitions', false, @islogical);
            p.addParameter('partialLoadSpec', [], @(x) isempty(x) || isstruct(x));
            p.addParameter('use_parallel', ~isempty(gcp('nocreate')), @islogical);
            
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            progress = p.Results.progress && exist('ProgressBar', 'class');
            use_parallel = p.Results.use_parallel;
            
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            
            % check that partitions are found
            partitionsAvailable = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(locationName);
            if p.Results.loadAllPartitions
                partitions = partitionsAvailable;
            else
                partitions = string(p.Results.partitions);
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
                if ~exist(file, 'file')
                    warning('Partition meta file %s not found', file);
                else
                    loaded = load(file); % assumes less than 2 GB per element, but much faster
                    partitionMeta.(partitions{iP}) = loaded.partitionMeta;
                end
            end
            partitions = fieldnames(partitionMeta); % drop the ones we couldn't load
            
            [N, elementsPerChunk] = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationName);
            
            partialLoadSpec = p.Results.partialLoadSpec;
            if ~isempty(partialLoadSpec) && ~isempty(fieldnames(partialLoadSpec))
                partialLoadDataPath = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartialLoadDataFileName(locationName);
                if exist(partialLoadDataPath, 'file') > 0
                    ld = load(partialLoadDataPath, 'partialLoadData');
                    partialLoadData = ld.partialLoadData;
                    partialLoadMask = TrialDataUtilities.Data.SaveArrayIndividualized.buildPartialLoadMask(N, partialLoadData, partialLoadSpec);
                else
                    partialLoadMask = true(N, 1);
                    warning('No partialLoadData found in folder, ignoring partialLoadSpec');
                end
            else
                partialLoadMask = true(N, 1);
            end
            
            if p.Results.maxElements < N
                keep = find(partialLoadMask, p.Results.maxElements, 'first');
                partialLoadMask = false(N, 1);
                partialLoadMask(keep) = true;
            end
            
            if ~isempty(p.Results.message)
                str = p.Results.message;
            else
                str = sprintf('Loading from %s/', locationName);
            end
            
            % figure out which chunks to load
            whichChunk = TrialDataUtilities.Data.SaveArrayIndividualized.buildElementChunkLists(N, elementsPerChunk);
            nChunks = whichChunk(end);
            chunksToLoad = unique(whichChunk(partialLoadMask)); % critical that this ends up in sorted order
            nChunksToLoad = numel(chunksToLoad);
            
            partialLoadMaskByChunk = cell(nChunksToLoad, 1);
            for iC = 1:nChunksToLoad
                partialLoadMaskByChunk{iC} = partialLoadMask(whichChunk == chunksToLoad(iC));
            end
            
            if progress
                prog = ProgressBar(nChunksToLoad, str);
            else
                prog.update = @(ind) true;
            end
            loadedChunks = cell(nChunksToLoad, 1);
            
            if use_parallel
                if progress, prog.enableParallel(); end
                parfor iC = 1:nChunksToLoad
                    chunk_ind = chunksToLoad(iC);
                    loadedChunks{iC} = TrialDataUtilities.Data.SaveArrayIndividualized.internal_loadChunk(locationName, nChunks, chunk_ind, callbackFn, partitions, partialLoadMaskByChunk{iC});
                    if progress, prog.update(iC); end %#ok<PFBNS>
                end
            else
                for iC = 1:nChunksToLoad
                    chunk_ind = chunksToLoad(iC);
                    loadedChunks{iC} = TrialDataUtilities.Data.SaveArrayIndividualized.internal_loadChunk(locationName, nChunks, chunk_ind, callbackFn, partitions, partialLoadMaskByChunk{iC});
                    if progress, prog.update(iC); end
                end
            end
            if progress, prog.finish(); end
           
            S = cat(1, loadedChunks{:});
        end
        
        function element = internal_loadChunk(locationName, nChunks, chunk_ind, callbackFn, partitions, partialLoadMask_this_chunk)
            if isempty(callbackFn)
                file = TrialDataUtilities.Data.SaveArrayIndividualized.generateElementFileName(locationName, nChunks, chunk_ind);   
                loaded = load(file, 'element');
                element = loaded.element;
            else
                element = callbackFn(locationName, chunk_ind, '');
            end
                
            nPartitions = numel(partitions);
            % load in any partition parts
            if nPartitions > 0
                partData = cell(nPartitions, 1);
                for iP = 1:nPartitions
                    if isempty(callbackFn)
                        file = TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionElementFileName(locationName, nChunks, chunk_ind, partitions{iP});
                        loaded = load(file); % assumes less than 2 GB per element, but much faster
                        partData{iP} = loaded.partData;
                    else
                        % pass this element, the location, the id number, and the partition to the callback
                        partData{iP} = callbackFn(locationName, chunk_ind, partitions{iP});
                    end
                end
                element = structMerge(element, partData{:});
            end
                
            % apply partial load mask
            element = TrialDataUtilities.makecol(element(partialLoadMask_this_chunk));
                
            function s = structMerge(varargin)
                s2c = cellfun(@struct2cell, varargin, 'UniformOutput', false);
                flds = cellfun(@fieldnames, varargin, 'UniformOutput', false);
                s = cell2struct(cat(1, s2c{:}), cat(1, flds{:}), 1);
            end
        end
        
        function partitionMeta = loadPartitionMeta(locationName, varargin)
            p = inputParser();
            p.addParameter('partitions', {}, @(x) TrialDataUtilities.String.isstringlike(x));
            p.addParameter('loadAllPartitions', [], @(x) isempty(x) || islogical(x));
            p.addParameter('ignoreMissingPartitions', false, @islogical);
            p.parse(varargin{:});
            
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            
            % figure out which partitions to load, loadAll only if none specified manually
            partitionsAvailable = TrialDataUtilities.Data.SaveArrayIndividualized.listPartitions(locationName);
            loadAllPartitions = p.Results.loadAllPartitions;
            partitions = string(p.Results.partitions);
            if isempty(loadAllPartitions)
                loadAllPartitions = isempty(partitions);
            end
            if loadAllPartitions
                partitions = partitionsAvailable;
            else
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
                if ~exist(file, 'file')
                    warning('Could not find meta info file %s\n', file);
                else
                    loaded = load(file); % assumes less than 2 GB per element, but much faster
                    partitionMeta.(partitions{iP}) = loaded.partitionMeta;
                end
            end
        end
        
        function [whichChunk, elementsByChunk] = buildElementChunkLists(N, elementsPerChunk)
            whichChunk = floor((0:N-1)' / elementsPerChunk) + 1; 
            nChunks = whichChunk(end);
            elementsByChunk = arrayfun(@(c) find(whichChunk == c), (1:nChunks)', 'UniformOutput', false);
        end
        
        function mask = buildPartialLoadMask(N, partialLoadData, partialLoadSpec)
            mask = true(N, 1);
            flds = fieldnames(partialLoadSpec);
            for iF = 1:numel(flds)
                fld = flds{iF};
                keepVals = partialLoadSpec.(fld);
                if ~isfield(partialLoadData, fld)
                    warning('Partial load field %s not found in cached partialLoad data on disk, ignoring', fld);
                    continue;
                end
                vals = partialLoadData.(fld);
                assert(numel(vals) == N, 'Number of trials in partialLoadData.%s does not match number of saved elements (%d)', fld, N);
                
                tf = ismember(vals, keepVals);
                mask = mask & tf;
            end
        end         
        
        function [N, elementsPerChunk] = getArrayCount(locationName)
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            
            % get element count from count.txt file
            countFname = TrialDataUtilities.Data.SaveArrayIndividualized.generateCountFileName(locationName);
            countFid = fopen(countFname, 'r');
            if countFid == -1
                error('Could not find count file %s', countFname);
            end
            tokens = textscan(countFid, '%d', 2);
            fclose(countFid);
            vals = double(tokens{1});
            if isempty(vals)
                error('Could not read count file %s', countFname);
            end
            N = vals(1);
            if numel(vals) > 1
                elementsPerChunk = vals(2);
            else
                % wasn't saved in older versions
                elementsPerChunk = 1;
            end
        end
        
        function list = listPartitions(locationName)
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            
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
            list = string(list);
        end
        
        function fields = listPartialLoadFields(locationName)
           file =  TrialDataUtilities.Data.SaveArrayIndividualized.generatePartialLoadDataFileName(locationName);
           if exist(file, 'file') > 0
               ld = load(file, 'partialLoadData');
               partialLoadData = ld.partialLoadData;
               fields = string(fieldnames(partialLoadData));
           else
               fields = strings(0);
           end
        end
        
        function linkPartitionFromOtherLocation(locationNameRef, locationNameSave, varargin)
            % symlinkPartitionFromOtherLocation(locationNameRef, locationNameSave, 'partitions', {'partitionName'}, 'linkAllPartitions', [true/false])
            %
            % symbolically links a partition saved in locationNameRef into the location in locationNameSave

            p = inputParser();
            p.addParameter('mode', 'symlink', @ischar); % 'symlink' or 'copy'
            p.addParameter('overwrite', false, @islogical);
            p.addParameter('partitions', {}, @TrialDataUtilities.String.isstringlike);
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
            [nSave, nPerChunkSave] = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationNameSave);
            [nRef, nPerChunkRef] = TrialDataUtilities.Data.SaveArrayIndividualized.getArrayCount(locationNameRef);
            assert(nSave == nRef && nPerChunkSave == nPerChunkRef, 'Save location has %d (%d) elements but referenced location has %d (%d) elements', ...
                nSave, nPerChunkSave, nRef, nPerChunkRef);
            N = nSave;
            C = ceil(nSave / nPerChunkSave);
            
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
            
            % do the symlinking for each chunk
            prog = ProgressBar(C, 'Symlinking partition files by trial');
            for i = 1:C
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
                deleteIfPresent(fullfile(locationName, 'partialLoad.mat'));
                deleteIfPresent(fullfile(locationName, 'partitionMeta_*.mat'));
                deleteIfPresent(TrialDataUtilities.Data.SaveArrayIndividualized.generateCountFileName(locationName));
                deleteIfPresent(TrialDataUtilities.Data.SaveArrayIndividualized.generatePartitionListFileName(locationName));
            end
            
            function deleteIfPresent(file)
                if any(file == '*') || exist(file, 'file')
                    delete(file);
                end
            end
        end
    
        % general use meta, keyed by suffix with default suffix=""
        function meta = loadMeta(locationName, varargin)
            p = inputParser();
            p.addParameter('suffix', "", @(x) TrialDataUtilities.String.isstringlike(x));
            p.parse(varargin{:});
            
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            file = TrialDataUtilities.Data.SaveArrayIndividualized.generateMetaFileName(locationName, p.Results.suffix);
                
            if ~exist(file, 'file')
                error('Could not find meta file %s\n', file);
            else
                loaded = load(file); % assumes less than 2 GB per element, but much faster
                meta = loaded.meta;
            end
        end

        function saveMeta(locationName, meta, varargin)
            p = inputParser();
            p.addParameter('suffix', "", @(x) TrialDataUtilities.String.isstringlike(x));
            p.parse(varargin{:});
            
            locationName = TrialDataUtilities.Path.GetFullPath(locationName);
            file = TrialDataUtilities.Data.SaveArrayIndividualized.generateMetaFileName(locationName, p.Results.suffix);
                
            TrialDataUtilities.Save.mkdirRecursive(locationName);
            save(file, 'meta');
        end
    end

    % numerical array tools
    methods(Static)
        function saveTensorSplit(locationName, S, varargin)
            p = inputParser();
            p.addParameter('dim', 1, @isscalar); % which dimension to split along
            p.addParameter('elementsPerChunk', 1, @isscalar); % split elements into files with this many elements per file
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            dim = p.Results.dim;
            nPer = p.Results.elementsPerChunk;
            sz = size(S, dim);
            nTotal = ceil(sz/nPer);
            
            nEach = repmat(nPer, nTotal, 1);
            nEach(end) = sz - nPer * (nTotal - 1);

            split = TensorUtils.splitAlongDimension(S, dim, nEach);

            TrialDataUtilities.Data.SaveArrayIndividualized.saveArray(locationName, split, 'elementsPerChunk', 1, p.Unmatched);

            meta.dim = dim;
            meta.elementsPerChunk = nPer;
            TrialDataUtilities.Data.SaveArrayIndividualized.saveMeta(locationName, meta, 'suffix', "tensorSplit");
        end

        function [S, partitionMeta, partialLoadMask] = loadTensorSplit(locationName, varargin)
            [split, partitionMeta, partialLoadMask] = TrialDataUtilities.Data.SaveArrayIndividualized.loadArray(locationName, varargin{:});
            meta = TrialDataUtilities.Data.SaveArrayIndividualized.loadMeta(locationName, 'suffix', "tensorSplit");

            S = cat(meta.dim, split{:});
        end
    end
    
    methods(Static, Hidden)
        function fname = generateCountFileName(locationName)
            fname = fullfile(locationName, 'count.txt');
        end
        
        function fname = generatePartialLoadDataFileName(locationName)
            fname = fullfile(locationName, 'partialLoad.mat');
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
        
        % for individualized files
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

        function fname = generateMetaFileName(locationName, suffix)
            if nargin < 2
                suffix = "";
            end
            if strcmp(suffix, "")
                fname = fullfile(string(locationName), "meta.mat");
            else
                fname = fullfile(string(locationName), sprintf("meta_%s.mat", suffix));
            end
        end
    end
end
