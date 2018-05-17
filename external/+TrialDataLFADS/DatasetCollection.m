classdef DatasetCollection < LFADS.DatasetCollection
    properties(Dependent)
        trialDataSet % loaded from database
    end
    
    methods
        function ds = DatasetCollection(path)
            ds = ds@LFADS.DatasetCollection(path);
        end

        function autoDetectDatasets(dc)
            dc.clearDatasets;
            
            dsClassName = strrep(class(dc), '.DatasetCollection', '.Dataset');
            dsClass = str2func(dsClassName);
            
            % automatically find all files within and build datasets
            files = dir(dc.path);
            for iF = 1:numel(files)
                if strncmp(files(iF).name, '.', 1), continue, end % skip . directories
                if ~files(iF).isdir, continue, end % skip non-directories
                ff = fullfile(dc.path, files(iF).name); 
                if ~TrialData.loadFastIsValidLocation(ff), continue, end % test for TrialData.loadFast compatible
                
                ds = dsClass(dc, files(iF).name); %#ok<NASGU>
            end
        end
        
        function addDataset(dc, ds)
            assert(isa(ds, 'TrialDataLFADS.Dataset'), 'Must be TrialDataLFADS.Dataset instance');
            addDataset@LFADS.DatasetCollection(dc, ds);
        end
    end
    
    methods
        function tdSet = get.trialDataSet(dc)
            tdSet = cellvec(dc.nDatasets);
            for iDS = 1:dc.nDatasets
                tdSet{iDS} = dc.datasets(iDS).trialData;
            end
        end
            
        function tdSet = loadTrialDataSet(dc, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical);
            p.parse(varargin{:});
            
            if ~isempty(dc.trialDataSet) && ~p.Results.reload
                tdSet = dc.trialDataSet;
                maskLoad = cellfun(@isempty, tdSet);
            else
                tdSet = cell(dc.nDatasets, 1);
                maskLoad = true(dc.nDatasets, 1);
            end
            
            if any(maskLoad)
                prog = ProgressBar(dc.nDatasets, 'Loading trialData for each dataset');
                for i = 1:dc.nDatasets
                    prog.update(i);
                    tdSet{i} = dc.datasets(i).loadData(maskLoad(i)); % will cache in their trialData field
                end
                prog.finish();
            end
        end
    end
end
