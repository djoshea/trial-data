classdef Dataset < LFADS.Dataset
    
    properties
        trialData
    end
    
    methods
        function ds = Dataset(collection, tdPathOrShortName, td) 
            ds = ds@LFADS.Dataset(collection, tdPathOrShortName);
            if nargin > 2 && ~isempty(td)
                ds.trialData =  TrialDataConditionAlign(td);
                ds.trialData = ds.trialData.setMetaKey('datasetPath', '');
            end 
        end
        
        function nChannels = getNumChannelsForCounts(ds, td) %#ok<INUSL>
            nChannels = numel(td.listSpikeChannels());
        end
        
        function loadInfo(ds)
            % obtains metadata quickly from meta keys of trial data object
            
            if ds.infoLoaded, return; end

            tdMeta = TrialData.loadFastMetaOnly(ds.path);
            
            ds.subject = tdMeta.getMetaKey('subject');
            ds.saveTags = tdMeta.getMetaKey('saveTags');
            ds.datenum = tdMeta.getMetaKey('datenum');
            ds.nChannels = ds.getNumChannelsForCounts(tdMeta);
            ds.nTrials = tdMeta.nTrials;
            
            ds.infoLoaded = true;
        end
        
        function td = loadData(ds, reload)
            if nargin < 2
                reload = false;
            end
            
            if isempty(ds.trialData) || reload
                td = TrialData.loadFast(ds.path);
                td = TrialDataConditionAlign(td);
                td.datasetName = ds.name;
                td = td.setMetaKey('datasetPath', ds.path);
                ds.trialData = td;
            end
                
            td = ds.trialData;
        end
    end
end
