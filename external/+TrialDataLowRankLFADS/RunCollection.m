classdef RunCollection < LowRankLFADS.RunCollection
    methods
        function rc = RunCollection(varargin)
            rc@LowRankLFADS.RunCollection(varargin{:});
        end
        
        function loadTrialDataFromDatasetCollection(rc)
            % Call `loadPosteriorMeans` on each run in this collection
            prog = ProgressBar(rc.nRunsTotal, 'Loading trialData from datasetCollection');
            for i = 1:rc.nRunsTotal
                prog.update(i);
                rc.runs(i).loadTrialDataFromDatasetCollection();
            end
            prog.finish();
        end
        
         function addPosteriorMeansToTrialData(rc, varargin)
            % Call `loadPosteriorMeans` on each run in this collection
            prog = ProgressBar(rc.nRunsTotal, 'Adding PosteriorMeans to trialData');
            for i = 1:rc.nRunsTotal
                prog.update(i);
                rc.runs(i).addPosteriorMeansToTrialData(varargin{:});
            end
            prog.finish();
        end
    end
end