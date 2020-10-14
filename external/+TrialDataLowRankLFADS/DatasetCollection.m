classdef DatasetCollection < LowRankLFADS.DatasetCollection
    properties(Dependent)
        trialDataSet % loaded from database
    end
    
    methods
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
            assert(isa(ds, 'TrialDataLowRankLFADS.Dataset'), 'Must be TrialDataLowRankLFADS.Dataset instance');
            addDataset@LowRankLFADS.DatasetCollection(dc, ds);
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
            % pierreSettingsEric
            % db = LFADS_PierreExport.Database.loadDb();
            
            p = inputParser();
            p.addOptional('reload', false, @islogical);
            p.addParameter('datasetIdx', 1:dc.nDatasets, @isvector);
            p.parse(varargin{:});
            datasetIdx = LFADS.Utils.vectorMaskToIndices(p.Results.datasetIdx);
            
            if ~isempty(dc.trialDataSet) && ~any(cellfun(@isempty, dc.trialDataSet)) && ~p.Results.reload
                tdSet = dc.trialDataSet(datasetIdx);
                return;
            end
            
            tdSet = cell(numel(datasetIdx), 1);
            prog = ProgressBar(numel(datasetIdx), 'Loading trialData for each dataset');
            for iiDS = 1:numel(datasetIdx)
                prog.update(iiDS);
                tdSet{iiDS} = dc.datasets(datasetIdx(iiDS)).loadData(p.Results.reload); % will cache in their trialData field
            end
            prog.finish();

            % note: we don't set trialDataSet here, since it retrieves from
            % the datasets already as a Dependent property
        end
    end
end
