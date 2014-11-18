classdef SaveStructArrayIndividualized < handle

    methods(Static)            
        function saveStructArray(locationName, S, varargin)
        % saveStructArrayIndividualized(fname, S, name, varargin)
        % save a struct array S(:) as name1, name2, name3, ...
        % allowing individual elements to be loaded quickly

            assert(isvector(S));
            N = numel(S);
            
            % create the directory as path/name
            fullPath = locationName;
            mkdirRecursive(fullPath);
            
            prog = ProgressBar(N, 'Saving to %s/', fullPath);
            for i = 1:N 
                % save this element
                fname = sprintf('el%06d.mat', i);
                file = fullfile(fullPath, fname);
                element = S(i); %#ok<NASGU>
                savefast(file, 'element');
                prog.update(i);
            end
            prog.finish();
            
            % create count file containing just N
            countFid = fopen(fullfile(locationName, 'count.txt'), 'w');
            fprintf(countFid, '%d\n', N); 
            fclose(countFid);
        end
        
        function S = loadStructArray(locationName, varargin)
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
            
            prog = ProgressBar(N, 'Loading from %s/', locationName);
            for i = N:-1:1 % reverse order to preallocate array
                fname = sprintf('el%06d.mat', i);
                file = fullfile(locationName, fname);
                loaded = load(file, 'element');
                S(i) = loaded.element;
                prog.update(N-i+1);
            end
            prog.finish();
            
        end
        
    end
end
