classdef SaveArrayIndividualized < handle

    methods(Static)            
        function saveArray(locationName, S, varargin)
        % saveStructArrayIndividualized(fname, S, name, varargin)
        % save a struct array S(:) as name1, name2, name3, ...
        % allowing individual elements to be loaded quickly

            p = inputParser();
            p.addParameter('message', '', @ischar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            
            assert(isempty(S) || isvector(S));
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
                % save this element
                if isempty(callbackFn)
                    fname = sprintf('el%06d.mat', i);
                    file = fullfile(fullPath, fname);
                    element = S(i); %#ok<NASGU>
                    save(file, '-v6', 'element'); % assumes less than 2 GB per element, but much faster
                else
                    % pass this element, the location, and the id number to the callback
                    element = S(i);
                    callbackFn(element, fullPath, i);
                end
                prog.update(i);
            end
            prog.finish();
            
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
        
        function S = loadArray(locationName, varargin)
            p = inputParser();
            p.addParameter('message', '', @ischar);
            p.addParameter('callbackFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.parse(varargin{:});
            
            callbackFn = p.Results.callbackFn;
            
            locationName = GetFullPath(locationName);
            
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
                    S(i) = loaded.element;
                else
                    S(i) = callbackFn(locationName, i);
                end
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
    end
end
