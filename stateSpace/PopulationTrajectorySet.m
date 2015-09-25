classdef PopulationTrajectorySet
% class for storing a set of population trajectories, one for each of multiple conditions
% or possibly one for each of many trials
%
% The bases may be units (e.g. trial-averaged firing rates by condition by unit), or 
% components of a low-dimensional projection (e.g. PCA), or analog channels
% (i.e. the results of some other computation)
%
% The set of conditions is specified by a ConditionDescriptor
% Common time alignment is specified by an AlignDescriptor or a set of
% AlignDescriptors for multiple alignments
%
% This class uses a design pattern which utilizes the "on demand cache".
% The outer class (PTS) is a value class, but the odc is a handle class.
% This inner class exists to allow PTS to compute properties dynamically on
% the fly (like any dependent property), but to cache the results
% persistently within the odc, so that they need not be recomputed each
% time they are accessed. Manual overrides are also permitted for some
% properties, such that a stored property value is accessed rather than
% computing a dynamic property value on access.
%
% This pattern requires a coordinated set of .get, .set and .build methods 
% for each property. The set method (protected) will store the result
% inside the odc, which must be copied (cloned) before each write, so that
% different PTS instances share no common resources and thus act like value
% classes ought to behave. The get method checks whether the odc's stored 
% value is empty. If it is not empty, it returns that value. If it is empty
% it delegates to the .build method to compute the value of that property.
% Note that .build methods may be overriden in subclassses, even though
% .get accessors may not be. To flush or invalidate the value of a property
% simply set its value to []. In some cases, a property will have a manual
% override setting, which takes precedence over the computed value. The
% .get accessor will check the corresponding manual property, named
% .manualProperty, and return that value if not empty.

    properties
        datasetName
        
        dataUnits = '';
    end

    properties(SetAccess=protected, Hidden)
        % odc is an instance of PopulationTrajectorySetOnDemandCache
        % which is a handle class. Properties which derive from 
        odc
    end
    
    % properties which control the behavior of the pset and will invalidate
    % computed values
    % see .initialize() for default values
    properties
        % The following parameters affect data extraction:
        
        % numeric scalar indicating spacing between successive time points
        % all channels timeseries will be interpolated to a time vector
        % with this spacing
        timeDelta
        
        % The following parameters affect trial-averaging:
        
        % The minimum number of trials over which to compute a trial
        % average. This parameter determines the valid time windows for
        % trial-averaged data (e.g. dataMean)
        minTrialsForTrialAveraging
        
        % The minimum fraction of trials in a given condition over which to
        % compute a trial average, relative to the the total number of trials
        % in that condition. This parameter determines the valid time
        % windows for trial-averaged data (e.g. dataMean)
        minFractionTrialsForTrialAveraging 
        
        % quantile of trial data to use for dataIntervalLow
        dataIntervalQuantileLow

        % quantile of trial data to use for dataIntervalHigh
        dataIntervalQuantileHigh
    end
    
    % alignDescriptors, conditionDescriptor, and translationNormalization
    % which align, group, and translate/normalize the data generated from
    % the dataSources. These may be set using eponymous methods prefixed with set*
    properties(SetAccess=?PopulationTrajectorySetBuilder)
        % nAlign x 1 cell of alignDescriptors
        alignDescriptorSet = {};
        
        % nAlign-1 x 1 vector of time gaps between successive alignments
        % initially empty, indicating that no gap has been specified
        interAlignGaps

        % ConditionDescriptor instance describing condition information 
        conditionDescriptor
        
        % StateSpaceTranslationNormalization instance describing the
        % translation and normalization to apply to each basis.
        % This will be applied during buildDataByTrial for dataSourceManual
        % = false psets (thus being reflected in the trial-averages automatically)
        % or applied manually to dataByTrial (if non-empty) and dataMean
        % (if non-empty) for dataSourceManual = false
        translationNormalization
        
        % SpikeFilter instance to use when converting spiking units to
        % firing rate channels
        spikeFilter
    end
    
    % randomized data we don't want to clear with the ODC because it's so
    % costly to compute
    properties(SetAccess=?PopulationTrajectorySetBuilder)  
        % ConditionDescriptor for dataMeanRandomized describing the
        % randomization method used
        conditionDescriptorRandomized
        
        % nAlign x 1 cell with nBases x nCondition x nTime(iAlign) x nRandomSamples numeric
        % tensors containing randomly generated dataMean
        dataMeanRandomized
        dataSemRandomized
    end

    % These properties store raw data sources (TrialDataConditionAlign instances)
    % from which data is extracted, as well as track from where each basis
    % originates. Along with other misc settings
    properties(SetAccess=?PopulationTrajectorySetBuilder)
        timeUnitName % string name of common time units
        
        timeUnitsPerSecond % scalar conversion factor
        
        % boolean scalar
        % is trial averaged data computed from byTrial data (false) or
        % specified manually (true). This flag controls whether property
        % values are computed dynamically from the dataSources and stored
        % in the .odc, or stored persistently in .manualData
        dataSourceManual
        
        % TrialData data sources which source all data for the trajectories
        % this may be a single trial data object or many. If there is only
        % one, all bases are considered simultaneous.
        dataSources
        
        % nBases x 1 index into dataSourceSet.
        basisDataSourceIdx
        
        % nBases x 1 cellstr indicating which channel name to extract data
        % from
        basisDataSourceChannelNames
        
        % random seed used as initial seed when generating random data
        % sets
        randomSeed
                        
        % nConditions x 1 logical: conditions not listed valid will be
        % ignored (and thus NaN) in .dataMean when it is computed
        % this can be used if some bases but not all bases have data for
        % that condition, and the condition is messing up the time vector
        % limits
        conditionIncludeMask
    end
    
    % Properties whose values are computed dynamically and persist within odc
    % or are specified manually and persist within manualData
    properties(Dependent, Transient, SetAccess=?PopulationTrajectorySetBuilder)
        % Alignment summary statistics by basis by align
        
        % nAlignSummaryData (originally nDataSources) x nAlign cell containing AlignSummary instances for each
        % basis. built by buildAlignSummary
        alignSummaryData
        
        % nAlign x 1 cell containing an alignSummary instance for
        % each align, aggregated over all bases
        alignSummaryAggregated
        
        % nBases x 1 vector of indices into alignSummaryData indicating
        % the AlignSummary instance which corresponds to each basis
        basisAlignSummaryLookup
        
        % Basis info 
        
        % nBases x 1 cellstr: names for each basis 
        basisNames = {};
        
        % nBases x 1 cellstr: units for each basis 
        basisUnits = {}
        
        % nBases x 1 logical: which bases are considered included in the
        % present analyses and time-vector computations. updating this will
        % invalidate data means. it is added to make basis filtering
        % operations simpler (since all psets can maintain the same logical
        % vector
        basisValid = []
        
        basisInvalidCause = {};
        
        %% BELOW ARE FOR INDIVIDUAL TRIAL DATA
        
        % data by trial cells contain common time vectors across trials,
        % which differ across bases and aligns
        
        % nBases x nAlign cell containing ordered data by trial
        % each cell contains nTrials x nTime analog data for that basis,
        % for that alignment, in order by trial number. The time vector
        % along the columns is given by 
        % nanmin(tMinByTrial{..}) : nanmax(tMaxByTrial{..})
        dataByTrial
        
        % nBases x nAlign numeric matrix containing the start and stop
        % times for the time vector indicating the time along the columns
        % of the corresponding cell of dataByTrial
        tMinForDataByTrial
        tMaxForDataByTrial
        
        % nBases x nAlign cell containing a nTrials x 1 vector indicating
        % whether each trial was considered valid for the given alignment
        alignValidByTrial
        
        % nBasis x nAlign cell arrays indicating the start and stop
        % timepoints for each trial in dataByTrial
        tMinByTrial
        tMaxByTrial
        
        % Below are for finding time windows to use for trial averaging
        % which occurs just before
        
        % nAlign x nBases x nConditions matrices containing the start and stop times for
        % which sufficient trials exist to compute a trial-average for each
        % align, basis, and condition
        tMinValidByAlignBasisCondition
        tMaxValidByAlignBasisCondition
        
        % nConditions x 1 logical vector indicating whether the condition
        % has trial average data for EVERY basis on ALL aligns
        conditionHasValidTrialAverageAllAlignsBases
        
        % nAlign x nConditions matrix containing the start and stop times
        % for which sufficient trials exist to compute a trial-average for 
        % ALL valid bases, on each align and condition. only Conditions for
        % which conditionHasValidTrialAverageAllAlignsBases is true are
        % considered. Only bases which are valid are considered
        tMinValidAllBasesByAlignCondition
        tMaxValidAllBasesByAlignCondition
 
        %% BELOW ARE FOR TRIAL-AVERAGED DATA WITHIN CONDITION
        
        % trial averaged data contains common time vectors across
        % conditions and bases, but which differ across aligns. Invalid
        % portions of the bases with smaller windows will be filled with
        % NaNs.
        
        % nAlign x 1 numeric vectors indicating the time window used for each
        % alignmnet for the dataMean and dataIntervalHigh/Low cells
        % (buildDataMean)
        tMinForDataMean
        tMaxForDataMean
        
        % trial-averaged data within each condition
        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign) 
        % numeric array of single trial-averaged traces
        % (buildDataMean)
        dataMean

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign) 
        % numeric array of single trial-averaged trace standard error
        % (buildDataMean)
        dataSem 
        
        % nAlign x nBases x nCondition scalar array indicating how many
        % trials contributed to data in each cell  (buildDataNTrials)
        dataNTrials
          
        % nBases x nConditions cell of trial (into dataByTrial) for each
        % condition (buildDataNTrials)
        trialLists
        
        % nAlign x nBases x nCondition logical array indicating whether there is  
        % valid data in the corresponding cell, which is required when the conditions
        % matrix isn't complete (all conditions valid) or only some alignments are valid
        % for a given condition, etc. (buildDataNTrials)
        dataValid

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the lower confidence interval value for
        % dataMean (buildDataRandomizedIntervals)
        dataIntervalLow
        
        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the upper confidence interval value for
        % dataMean (buildDataRandomizedIntervals)
        dataIntervalHigh
        
        % nBases x sum(nTimeDataMean) (aka TA) x nConditions numeric array
        % containing randomized samples of differences of pairs of trials,
        % scaled by 1/sqrt(2*nTrials). These serve as a sample from the
        % centered distribution of the estimation noise on the dataMean
        % traces. These are used by StateSpaceProjectionStatistics in
        % estimating a noise floor (buildDifferenceOfTrials)
        dataDifferenceOfTrialsScaledNoiseEstimate % built by buildDifferenceOfTrials
    end
    
    % generated on the fly on get.property. Best to retrieve this value
    % once rather than access it directly every iteration in a loop
    properties(Transient, Dependent, SetAccess=protected)
        % nAlign x 1 numeric vectors indicating the number of timepoints
        % used for each alignment in the dataMean cell
        nTimeDataMean
        
        % nAlign x 1 cell array of time vectors corresponding to the 3rd
        % dimension of dataMean{iAlign}
        tvecDataMean
        
        % nAlign x nBases cell of time vectors for each .dataByTrial matrix
        tvecDataByTrial
           
        % nAlign x nBases x nConditions numeric of number of non-nan
        % timepoints in each data mean
        nTimeValidByAlignBasisCondition
        
        % is dataByTrial empty?
        hasDataByTrial
        
        % has one of the storeDataRandomized* methods been called to
        % populate dataMeanRandomized?
        hasDataRandomized
        
        % number of randomized samples drawn when generating
        % dataMeanRandomized
        nRandomSamples
        
        % C x 1 logical indicating which conditions have trials for each
        % basis and alignment. Other conditions will not have any
        % trial-averaged data in .dataMean
        conditionsWithTrialsAllBasesAligns
    end
    
    % Properties within *Manual properties store manually-specified values for each of the
    % above properties. These are used to store persistent copies of
    % the data when .dataSourceManual is true
    properties(Hidden, Access=protected) 
        basisNamesManual
        basisUnitsManual
        basisValidManual
        basisInvalidCauseManual
        conditionIncludeMaskManual
        alignSummaryDataManual
        alignSummaryAggregatedManual
        basisAlignSummaryLookupManual
        dataByTrialManual
        tMinForDataByTrialManual
        tMaxForDataByTrialManual
        alignValidByTrialManual
        tMinValidByAlignBasisConditionManual
        tMaxValidByAlignBasisConditionManual
        tMinByTrialManual
        tMaxByTrialManual
        tMinForDataMeanManual
        tMaxForDataMeanManual
        nTimeDataMeanManual % this could easily be computed from tMin/tMaxForDataMean
        dataMeanManual
        dataSemManual
        dataIntervalHighManual
        dataIntervalLowManual
        dataNTrialsManual
        trialListsManual
        dataValidManual
        dataDifferenceOfTrialsScaledNoiseEstimateManual
    end
    
    % Dependent properties which we compute on the fly rather than cache
    properties(Dependent)
        % indicates how the elements in data are computed from individual trials
        % true means all trials correspond one-to-one with each other
        % false means trials were not simultaneous, implying that only
        % trial-averages should be compared
        simultaneous
        
        nTrials % only valid when simultaneous is true
        
        nTrialsValid % only valid when simultaneous is true
        
        % number of data sources used across all bases
        nDataSources
        
        nBases % number of bases (e.g. units, analog channels)
        
        nBasesValid % number of bases currently marked valid

        nConditions % number of conditions in conditionDescriptor
        
        nAxes % conditionDescriptor.nAxes
        
        nAlign % number of alignments in alignDescriptorSet
        
        conditions % shortcut to conditionDescriptor.conditions
        
        conditionsSize % pass-thru to .conditionDescriptor
        
        conditionsSizeNoExpand
        
        alignNames % names pulled from the alignDescriptors

        conditionNames % condition names pulled from conditionDescriptor
        
        dataIntervalQuantilesAsString % dataIntervalQuantileLow/High summarized as a string
        
        nAlignSummaryData % number of unique alignSummaryData instances (for each align)

        % nAlign x nCondition logical indicating which condition/alignments
        % have any trials for any bases
        alignConditionsWithTrials 
        
        % nAlign x nBases x nConditions logical indicating which bases have valid
        % trial averages on each condition
        alignBasisConditionsWithValidTrialAverage
        
        % nAlign x nConditions logical indicating which align/conditions
        % have trial averages on at least one basis
        alignConditionsWithTrialAverage
        
        % nBases x nConditions logical indicating which bases have valid
        % trial averages on each condition on ALL aligns
        basisConditionsWithValidTrialAverage
        
        % nBases x 1 logical array indicating bases which have no valid trials on some condition for which other bases have trials
        basesMissingTrialAverageForNonEmptyConditionAligns
        
        % nBases x 1 logical indicating which valid basess have at least one valid
        % trial average on some condition, on some alignment
        basesNonEmpty
        
        % nConditions x 1 logical indicating which conditions are present
        % in all nonEmptyBases
        conditionsWithValidTrialAverageOnNonEmptyBases
    end
   
    % Constructor
    methods
        function pset = PopulationTrajectorySet()
            pset = pset.initialize();
        end
    end
    
    % loadobj custom loading
    methods(Static)
        function pset = loadobj(pset)
            pset = builtin('loadobj', pset);
            if isstruct(pset)
                error('Pset loaded as struct, cannot load');
            end
            pset = pset.initialize();
        end
    end
    
    % saveobj custom saving
    methods
        function pset = saveobj(pset)

        end
    end
    
    % Saving and loading piecemeal
    methods
        function saveFast(pset, location, varargin)
            p = inputParser();
            p.addParameter('recursive', false, @islogical); % calls saveFast on each source too
            p.parse(varargin{:});
            
            sources = pset.dataSources;
            
            if ~pset.dataSourceManual
                pset.dataSources = 'saved separately, use loadFast';
            end
            
            % comment this out to not clear this to save the alignment we've done
            pset.odc = [];

            mkdirRecursive(location);
            savefast(fullfile(location, 'pset.mat'), 'pset');

            % save elements of sources
            if ~pset.dataSourceManual
                msg = sprintf('Saving PopulationTrajectorySet to %s', location);
                if p.Results.recursive
                    TrialDataUtilities.Data.SaveArrayIndividualized.saveArray(location, sources, 'message', msg, 'callbackFn', @saveCallback);
                else
                    TrialDataUtilities.Data.SaveArrayIndividualized.saveArray(location, sources, 'message', msg);
                end
            end
            
            function saveCallback(tdcaCell, location, i)
                sub = fullfile(location, sprintf('source%06d', i));
                tdcaCell{1}.saveFast(sub);
            end
        end
    end
    
    % loadFast
    methods(Static)
        function pset = loadFast(location, varargin)
            p = inputParser();
            p.addParameter('recursive', false, @islogical); % calls saveFast on each source too
            p.parse(varargin{:});
            
            if ~exist(location, 'dir')
                error('Directory %s not found. Did you save with saveFast?', location);
            end
            loaded = load(fullfile(location, 'pset.mat'));
            pset = loaded.pset;
            
            % load elements of sources
            if ~pset.dataSourceManual
                msg = sprintf('Loading PopulationTrajectorySet from %s', location);
                if p.Results.recursive
                    sources = TrialDataUtilities.Data.SaveArrayIndividualized.loadArray(location, 'message', msg, 'callbackFn', @loadCallback);
                else
                    sources = TrialDataUtilities.Data.SaveArrayIndividualized.loadArray(location, 'message', msg);
                end
                pset.dataSources = sources;
            end
            
            pset = pset.initialize();
            
            function tdcaCell = loadCallback(location, i)
                % need to wrap in cell because of the way save/load array
                % work
                sub = fullfile(location, sprintf('source%06d', i));
                tdcaCell = { TrialData.loadFast(sub) };
            end
        end
    end
    
    % save / load wrappers for CacheCustomSaveLoad
    methods
        function tf = getUseCustomSaveLoad(pset, info) %#ok<INUSD>
            % use save custom only if I have data sources that could be
            % saved separately. If I only carry pre-extracted data, we'll
            % save as one entity using the normal save process
            tf = ~pset.dataSourceManual;
        end
        
        function token = saveCustomToLocation(pset, location)
            pset.saveFast(location);
            token = [];
        end
    end
    
    % load for CacheCustomSaveLoad
    methods(Static)
        function data = loadCustomFromLocation(location, token) %#ok<INUSD>
            data = PopulationTrajectorySet.loadFast(location);
        end
    end
    
    % initialization, cache invalidation
    methods
        function pset = initialize(pset)
            pset.warnIfNoArgOut(nargout);
            
            % for all empty properties which must be initialized, assign a
            % reasonable initial value. However, do not overwrite existing
            % values since this is a ConstructOnLoad class
            if isempty(pset.odc)
                pset.odc = PopulationTrajectorySetOnDemandCache();
            end
            
            if isempty(pset.alignDescriptorSet)
                pset.alignDescriptorSet = {AlignDescriptor()};
                pset = pset.applyAlignDescriptorSet();
            end
            
            if isempty(pset.conditionDescriptor)
                pset.conditionDescriptor = ConditionDescriptor();
                pset = pset.applyConditionDescriptor();
            end
            
%             if isempty(pset.translationNormalization) || pset.translationNormalization.nBases ~= pset.nBases
%                 pset.translationNormalization = ...
%                     StateSpaceTranslationNormalization.buildIdentityForPopulationTrajectorySet(pset);
%             end
            
            if isempty(pset.timeDelta)
                pset.timeDelta = 1;
            end 
            
            if isempty(pset.spikeFilter)
                pset.spikeFilter = SpikeFilter.getDefaultFilter();
            end
            
            if isempty(pset.minTrialsForTrialAveraging)
                pset.minTrialsForTrialAveraging = 1; % default to requiring only 1 trial
            end
            
            if isempty(pset.minFractionTrialsForTrialAveraging)
                pset.minFractionTrialsForTrialAveraging = 0; % default to requiring only 1 trials
            end
            
            if isempty(pset.dataIntervalQuantileLow)
                pset.dataIntervalQuantileLow = 0.025;
            end
            
            if isempty(pset.dataIntervalQuantileHigh)
                pset.dataIntervalQuantileHigh = 0.975;
            end
            
            if isempty(pset.dataSourceManual)
                pset.dataSourceManual = false;
            end
            
            if isempty(pset.dataSources)
                pset.dataSources = {};
            end
            
            if isempty(pset.basisDataSourceChannelNames)
                pset.basisDataSourceChannelNames = {};
            end
            
            pset = pset.updateValid();
            
            pset = pset.invalidateCache('clearDataRandom', false);
        end
            
        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function pset = invalidateCache(pset, varargin)
            p = inputParser();
            p.addParameter('clearDataRandom', true, @islogical);
            p.parse(varargin{:});
            
            pset.warnIfNoArgOut(nargout);
            % copy before writing to odc!
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flush();
            end
            if p.Results.clearDataRandom
                pset = pset.invalidateRandomizedTrialAveragedData();
            end
        end
        
        function pset = invalidateTrialAveragedData(pset)
            pset.warnIfNoArgOut(nargout);
            % copy before writing to odc!
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushTrialAveragedData();
            end
            pset = pset.invalidateRandomizedTrialAveragedData();
        end

        function pset = invalidateRandomizedTrialAveragedData(pset)
            pset.warnIfNoArgOut(nargout);
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushRandomizedTrialAveragedData();
            end
            pset.dataMeanRandomized = {};
            pset.dataSemRandomized = {};
        end
        
        function pset = invalidateAlignSummaryData(pset)
            pset.warnIfNoArgOut(nargout);
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushAlignSummaryData();
            end
        end
    end
    
    % Display / description
    methods 
       function printDescription(pset)
           if pset.dataSourceManual
               dataSourceStr = 'manual stored data';
           else
               dataSourceStr = sprintf('%d data sources', pset.nDataSources);
           end
            tcprintf('inline', '{yellow}%s: {bright white}%d bases {red}(%d valid){bright white}, %d conditions, %d alignments, %s\n', ...
                class(pset), pset.nBases, pset.nBasesValid, pset.nConditions, pset.nAlign, dataSourceStr);
            tcprintf('inline', '{yellow}Dataset: {none}%s\n', pset.datasetName);
            
            if pset.simultaneous
                tcprintf('inline', '{yellow}Simultaneous: {none}%d trials {red}(%d valid)\n', pset.nTrials, pset.nTrialsValid);
            end
            tcprintf('inline', '{yellow}Trial Averaging: {none}at least %d and %g%% of trials for average\n', ...
                pset.minTrialsForTrialAveraging, pset.minFractionTrialsForTrialAveraging*100);
            tcprintf('inline', '{yellow}Time Delta: {none}%g ms, {yellow}Spike Filter: {none}%s\n', ...
                pset.timeDelta, pset.spikeFilter.getDescription);
            
            if pset.hasDataRandomized
                tcprintf('inline', '{yellow}Data Randomized: {none}%d random samples, {red}%s\n', ...
                    pset.nRandomSamples, pset.conditionDescriptorRandomized.randomizationDescription);
            end
            
            fprintf('\n');
            pset.conditionDescriptor.printDescription();
            for i = 1:pset.nAlign
                pset.alignDescriptorSet{i}.printDescription();
            end
            if ~isempty(pset.translationNormalization)
                pset.translationNormalization.printOneLineDescription();
            end
        end
        
        function disp(pset)
            pset.printDescription();
            fprintf('\n');
            %builtin('disp', pset);
        end 
    end
    
    % these methods are setters for property values which change the
    % behavior of the PTS. They automatically invalidate downstream cached
    % values that depend on the value of the property being set.
    methods
        function pset = setSpikeFilter(pset, f)
            % changing spikeFilter invalidates everything and we reapply
            % all alignDescriptors to get the padding right
            pset.warnIfNoArgOut(nargout);
            pset.spikeFilter = f;
            pset = pset.invalidateCache();
            pset = pset.setAlignDescriptorSet(pset.alignDescriptorSet);
        end
        
        function pset = set.timeDelta(pset, v)
            % changing timeDelta invalidates everything
            pset.timeDelta = v;
            pset = pset.invalidateCache();
        end
        
        function pset = set.minTrialsForTrialAveraging(pset, v)
            % only affects trial averaging
            pset.minTrialsForTrialAveraging = v;
            pset = pset.invalidateTrialAveragedData();
        end
        
        function pset = set.minFractionTrialsForTrialAveraging(pset, v)
            % only affects trial averaging
            pset.minFractionTrialsForTrialAveraging = v;
            pset = pset.invalidateTrialAveragedData();
        end
        
        function pset = set.dataIntervalQuantileLow(pset, v)
            pset.dataIntervalQuantileLow = v;
            pset = pset.invalidateRandomizedTrialAveragedData();
        end

        function pset = set.dataIntervalQuantileHigh(pset, v)
            pset.dataIntervalQuantileHigh = v;
            pset = pset.invalidateRandomizedTrialAveragedData();
        end
        
        function pset = set.conditionIncludeMask(pset, v)
            assert(islogical(v) && isvector(v) && numel(v) == pset.nConditions, ...
                'conditionIncludeMask must be logical vector with length nConditions');
            pset.conditionIncludeMaskManual = v;
        end
        
        function v = get.conditionIncludeMask(pset)
            if isempty(pset.conditionIncludeMaskManual)
                v = truevec(pset.nConditions);
            else
                v = pset.conditionIncludeMaskManual;
            end
        end
    end
    
    % General utility methods
    methods(Access=protected)
        function warnIfNoArgOut(obj, nargOut)
        % call using obj.warnIfNoArgOut(nargout);
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('%s is not a handle class. If the instance returned by this method is not stored, this call has no effect.', ...
                    class(obj));
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end
    
    % get and set properties from odc, deferring to build* methods for
    % initial computation and storage in the odc
    methods
        function v = get.alignSummaryData(pset)
            if ~pset.dataSourceManual
                v = pset.odc.alignSummaryData;
                if isempty(v)
                    pset.buildAlignSummaryData();
                    v = pset.odc.alignSummaryData;
                end
            else
                v = pset.alignSummaryDataManual;
            end
        end
        
        function pset = set.alignSummaryData(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.alignSummaryData = v;
            else
                pset.alignSummaryDataManual = v;
            end
        end
        
        function v = get.basisAlignSummaryLookup(pset)
            if ~pset.dataSourceManual
                v = pset.odc.basisAlignSummaryLookup;
                if isempty(v)
                    pset.buildAlignSummaryData();
                    v = pset.odc.basisAlignSummaryLookup;
                end
            else
                v = pset.basisAlignSummaryLookupManual;
            end
        end 
        
        function pset = set.basisAlignSummaryLookup(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.basisAlignSummaryLookup = v;
            else
                pset.basisAlignSummaryLookupManual = v;
            end
        end
        
        function v = get.alignSummaryAggregated(pset)
            %if ~pset.dataSourceManual
                v = pset.odc.alignSummaryAggregated;
                if isempty(v)
                    pset.buildAlignSummaryAggregated();
                    v = pset.odc.alignSummaryAggregated;
                end
%             else
%                 v = pset.alignSummaryAggregatedManual;
%             end
        end 
        
        function pset = set.alignSummaryAggregated(pset, v)
%             if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.alignSummaryAggregated = v;
%             else
%                pset.alignSummaryAggregatedManual = v;
%            end
        end
        
        function v = get.dataByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.dataByTrial;
                end
            else
                v = pset.dataByTrialManual;
            end
        end 
        
        function pset = set.dataByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataByTrial = v;
            else
                pset.dataByTrialManual = v;
            end
        end
        
        function v = get.tMinForDataByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMinForDataByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.tMinForDataByTrial;
                end
            else
                v = pset.tMinForDataByTrialManual;
            end
        end 
        
        function pset = set.tMinForDataByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMinForDataByTrial = v;
            else
                pset.tMinForDataByTrialManual = v;
            end
        end
        
        function v = get.tMaxForDataByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMaxForDataByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.tMaxForDataByTrial;
                end
            else
                v = pset.tMaxForDataByTrialManual;
            end
        end 
        
        function pset = set.tMaxForDataByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMaxForDataByTrial = v;
            else
                pset.tMaxForDataByTrialManual = v;
            end
        end
        
        function v = get.tvecDataByTrial(pset)
            % generate on the fly, no caching
            v = cell(pset.nBases,pset.nAlign);
            for iAlign = 1:pset.nAlign
                for iBasis = 1:pset.nBases
                    v{iBasis, iAlign} = makecol(pset.tMinForDataByTrial(iBasis,iAlign):pset.timeDelta:pset.tMaxForDataByTrial(iBasis,iAlign));
                end
            end
        end 
        
        function v = get.alignValidByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.alignValidByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.alignValidByTrial;
                end
            else
                v = pset.alignValidByTrialManual;
            end
        end 
        
        function pset = set.alignValidByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.alignValidByTrial = v;
            else
                pset.alignValidByTrialManual = v;
            end
        end
        
        function v = get.tMinByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMinByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.tMinByTrial;
                end
            else
                v = pset.tMinByTrialManual;
            end
        end 
        
        function pset = set.tMinByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMinByTrial = v;
            else
                pset.tMinByTrialManual = v;
            end
        end
        
        function v = get.tMaxByTrial(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMaxByTrial;
                if isempty(v)
                    pset.buildDataByTrial();
                    v = pset.odc.tMaxByTrial;
                end
            else
                v = pset.tMaxByTrialManual;
            end
        end 
        
        function pset = set.tMaxByTrial(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMaxByTrial = v;
            else
                pset.tMaxByTrialManual = v;
            end
        end
        
        function v = get.basisNames(pset)
            if ~pset.dataSourceManual
                v = pset.odc.basisNames;
                if isempty(v)
                    pset.buildBasisNamesUnits();
                    v = pset.odc.basisNames;
                end
            else
                v = pset.basisNamesManual;
            end
        end 
        
        function pset = set.basisNames(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.basisNames = v;
            else
                pset.basisNamesManual = v;
            end
        end
        
        function v = get.basisUnits(pset)
            if ~pset.dataSourceManual
                v = pset.odc.basisUnits;
                if isempty(v)
                    pset.buildBasisNamesUnits();
                    v = pset.odc.basisUnits;
                end
            else
                v = pset.basisUnitsManual;
            end
        end 
        
        function pset = set.basisUnits(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.basisUnits = v;
            else
                pset.basisUnitsManual = v;
            end
        end
        
        function v = get.basisValid(pset)
            v = pset.odc.basisValid;            
            if isempty(v)
                pset.buildBasisValid();
                v = pset.odc.basisValid;
            end
        end
        
        function pset = set.basisValid(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.basisValid= v;
            else
                pset.basisValidManual = v;
            end
        end
        
        
        function v = get.basisInvalidCause(pset)
            v = pset.odc.basisInvalidCause;            
            if isempty(v)
                pset.buildBasisValid();
                v = pset.odc.basisInvalidCause;
            end
        end
        
        function pset = set.basisInvalidCause(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.basisInvalidCause = v;
            else
                pset.basisInvalidCauseManual = v;
            end
        end
        
        function v = get.tMinValidByAlignBasisCondition(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMinValidByAlignBasisCondition;
                if isempty(v)
                    pset.buildTimeWindowsByAlignBasisCondition();
                    v = pset.odc.tMinValidByAlignBasisCondition;
                end
            else
                v = pset.tMaxValidByAlignBasisConditionManual;
            end
        end
        
        function pset = set.tMinValidByAlignBasisCondition(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMinValidByAlignBasisCondition = v;
            else
                pset.tMaxValidByAlignBasisConditionManual = v;
            end
        end
        
        function v = get.tMaxValidByAlignBasisCondition(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMaxValidByAlignBasisCondition;
                if isempty(v)
                    pset.buildTimeWindowsByAlignBasisCondition();
                    v = pset.odc.tMaxValidByAlignBasisCondition;
                end
            else
                v = pset.tMaxValidByAlignBasisConditionManual;
            end
        end
        
        function pset = set.tMaxValidByAlignBasisCondition(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMaxValidByAlignBasisCondition = v;
            else
                pset.tMaxValidByAlignBasisConditionManual = v;
            end
        end
        
        function v = get.conditionHasValidTrialAverageAllAlignsBases(pset)
            hasAvg = pset.tMinValidByAlignBasisCondition <= pset.tMaxValidByAlignBasisCondition;
            v = makecol(TensorUtils.allMultiDim(hasAvg(:, pset.basisValid, :), [1 2]));
        end
        
        function v = get.tMinValidAllBasesByAlignCondition(pset)
            % generate on the fly, no caching
            % only consider conditions that have the potential to
            % contribute data to the trial averages on all bases
            cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            
            if ~any(cMask)
                v = nan(pset.nAlign, pset.nConditions);
            else
                v = nanmax(pset.tMinValidByAlignBasisCondition(:, pset.basisValid, cMask), [], 2);
                % but then reexpand this to have the full complement of
                % conditions, using NaNs for non-contributing conditions
                v = TensorUtils.squeezeDims(TensorUtils.inflateMaskedTensor(v, 3, cMask, NaN), 2);
            end
        end
        
        function v = get.tMaxValidAllBasesByAlignCondition(pset)
            % generate on the fly, no caching
            cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            if ~any(cMask)
                v = nan(pset.nAlign, pset.nConditions);
            else
                v = nanmin(pset.tMaxValidByAlignBasisCondition(:, pset.basisValid, cMask), [], 2);
                v = TensorUtils.squeezeDims(TensorUtils.inflateMaskedTensor(v, 3, cMask, NaN), 2);
            end
        end
        
        function v = get.tvecDataMean(pset)
            % generate on the fly, no caching
            v = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                v{iAlign} = makecol(pset.tMinForDataMean(iAlign):pset.timeDelta:pset.tMaxForDataMean(iAlign));
            end
        end 
        
        function v = get.nTimeValidByAlignBasisCondition(pset)
            v = (pset.tMaxValidByAlignBasisCondition - pset.tMinValidByAlignBasisCondition) / pset.timeDelta + 1;
        end
        
        function v = get.nTimeDataMean(pset)
            v = makecol(cellfun(@numel, pset.tvecDataMean));
        end 
        
        function pset = set.nTimeDataMean(pset, v)
            pset.odc = pset.odc.copy();
            pset.odc.nTimeDataMean = v;
        end
        
        function v = get.tMinForDataMean(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMinForDataMean;
                if isempty(v)
                    pset.buildDataMean();
                    v = pset.odc.tMinForDataMean;
                end
            else
                v = pset.tMinForDataMeanManual;
            end
        end 
        
        function pset = set.tMinForDataMean(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMinForDataMean = v;
            else
                pset.tMinForDataMeanManual = v;
            end
        end
        
        function v = get.tMaxForDataMean(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMaxForDataMean;
                if isempty(v)
                    pset.buildDataMean();
                    v = pset.odc.tMaxForDataMean;
                end
            else
                v = pset.tMaxForDataMeanManual;
            end
        end 
        
        function pset = set.tMaxForDataMean(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.tMaxForDataMean = v;
            else
                pset.tMaxForDataMeanManual = v;
            end
        end
        
        function v = get.dataValid(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataValid;
                if isempty(v)
                    pset.buildDataNTrials();
                    v = pset.odc.dataValid;
                end
            else
                v = pset.dataValidManual;
            end
        end 
        
        function pset = set.dataValid(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataValid = v;
            else
                pset.dataValidManual = v;
            end
        end
        
        function v = get.dataMean(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataMean;
                if isempty(v)
                    pset.buildDataMean();
                    v = pset.odc.dataMean;
                end
            else
                v = pset.dataMeanManual;
            end
        end
        
        function pset = set.dataMean(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataMean = v;
            else
                pset.dataMeanManual = v;
            end
        end
        
        function v = get.dataSem(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataSem;
                if isempty(v)
                    pset.buildDataMean();
                    v = pset.odc.dataSem;
                end
            else
                v = pset.dataSemManual;
            end
        end
        
        function pset = set.dataSem(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataSem = v;
            else
                pset.dataSemManual = v;
            end
        end

        function v = get.dataNTrials(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataNTrials;
                if isempty(v)
                    pset.buildDataNTrials();
                    v = pset.odc.dataNTrials;
                end
            else
                v = pset.dataNTrialsManual;
            end
        end 
        
        function pset = set.dataNTrials(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataNTrials = v;
            else
                pset.dataNTrialsManual = v;
            end
        end
        
        
        function v = get.trialLists(pset)
            if ~pset.dataSourceManual
                v = pset.odc.trialLists;
                if isempty(v)
                    pset.buildDataNTrials();
                    v = pset.odc.trialLists;
                end
            else
                v = pset.trialListsManual;
            end
        end 
        
        function pset = set.trialLists(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.trialLists = v;
            else
                pset.trialListsManual = v;
            end
        end
        
        function v = get.dataIntervalLow(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataIntervalLow;
                if isempty(v) && pset.hasDataRandomized
                     pset.buildDataRandomizedIntervals();
                     v = pset.odc.dataIntervalLow;
                end
            else
                v = pset.dataIntervalLowManual;
            end
        end
        
        function pset = set.dataIntervalLow(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataIntervalLow = v;
            else
                pset.dataIntervalLowManual = v;
            end
        end

        function v = get.dataIntervalHigh(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataIntervalHigh;
                if isempty(v) && pset.hasDataRandomized
                    pset.buildDataRandomizedIntervals();
                    v = pset.odc.dataIntervalHigh;
                end
            else
                v = pset.dataIntervalHighManual;
            end
        end
        
        function pset = set.dataIntervalHigh(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataIntervalHigh = v;
            else
                pset.dataIntervalHighManual = v;
            end
        end
        
        
        
        function v = get.dataDifferenceOfTrialsScaledNoiseEstimate(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataDifferenceOfTrialsScaledNoiseEstimate;
                if isempty(v)
                    pset.buildDataNoiseEstimate();
                    v = pset.odc.dataDifferenceOfTrialsScaledNoiseEstimate;
                end
            else
                v = pset.dataDifferenceOfTrialsScaledNoiseEstimateManual;
            end
        end
        
        function pset = set.dataDifferenceOfTrialsScaledNoiseEstimate(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataDifferenceOfTrialsScaledNoiseEstimate = v;
            else
                pset.dataDifferenceOfTrialsScaledNoiseEstimateManual = v;
            end
        end
        
        function tf = get.hasDataByTrial(pset)
            tf = ~isempty(pset.dataByTrial);
        end
        
        % has one of the storeDataRandomized* methods been called to
        % populate dataMeanRandomized?
        function tf = get.hasDataRandomized(pset)
            tf = ~isempty(pset.dataMeanRandomized);
        end
        
        % number of randomized samples drawn when generating
        % dataMeanRandomized
        function n = get.nRandomSamples(pset)
            if ~pset.hasDataRandomized
                n = 0;
            else
                n = size(pset.dataMeanRandomized{1}, 4);
            end
        end
        
        function c = get.conditionsWithTrialsAllBasesAligns(pset)
            c = squeeze(all(all(pset.dataNTrials(:, pset.basisValid, :), 1), 2));
        end
    end
  
    % methods which set and apply the conditionDescriptor
    methods 
        function pset = setConditionDescriptor(pset, cd)
            % update the conditionDescriptor for all bases within
            % and clear any generated condition average data
            pset.warnIfNoArgOut(nargout);
            assert(ismember(class(cd), {'ConditionDescriptor', 'ConditionInfo'}), ...
                'Must be ConditionDescriptor instance');
            
            if isa(cd, 'ConditionInfo')
                cd = cd.fixAllValueLists().getConditionDescriptor();
            end
            
            assert(cd.allAxisValueListsManual, 'ConditionDescriptor must have manual axis value lists. Use .setAxisValueList or .fixValueListsByApplyingToTrialData');
            
            if pset.dataSourceManual
                assert(pset.conditionDescriptor.nConditions == cd.nConditions, ...
                    'Pset has manual data source, conditionDescriptor can only be replaced so as to preserve nConditions');
                warning('Replacing ConditionDescriptor for PopulationTrajectorySet with manual data source. This will only change the labeling of conditions, not the extracted data');
            end
            
            pset.conditionDescriptor = cd.getConditionDescriptor();
            
            pset = pset.applyConditionDescriptor();
        end
        
        function pset = applyConditionDescriptor(pset)
            % apply the condition descriptor to each trial data now and
            % store the resulting tdca
            pset.warnIfNoArgOut(nargout);
            
            prog = ProgressBar(pset.nDataSources, 'Grouping trials in data sources');
            conditionDescriptor = pset.conditionDescriptor;
            for iSrc = 1:pset.nDataSources
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.setConditionDescriptor(conditionDescriptor);
                prog.update(iSrc);
            end
            prog.finish();
            
            % only trial averaging and align summary needs to be done again
            pset = pset.invalidateTrialAveragedData();
            pset = pset.invalidateAlignSummaryData();
        end
    end
    
    % alignment: methods which directly update the align descriptor set of this pset
    methods
        function pset = setAlignDescriptorSet(pset, adSet)
            pset.warnIfNoArgOut(nargout);
            
            assert(~pset.dataSourceManual, 'PopulationTrajectorySets with manual data source cannot be realigned'); 
            
            if ~iscell(adSet)
                adSet = {adSet};
            end
            for i = 1:numel(adSet)
                if isa(adSet{i}, 'AlignInfo')
                    adSet{i} = AlignDescriptor.fromAlignDescriptor(adSet{i});
                end
            end
            assert(all(cellfun(@(x) isequal(class(x), 'AlignDescriptor') || ischar(x), ...
                adSet)), 'Must be AlignDescriptor instance or string');
            
            for i = 1:numel(adSet)
                if ischar(adSet{i})
                    adSet{i} = AlignDescriptor(adSet{i});
                end
            end
            
            % pre-pad all align descriptors based on spikeFilter for speed
            sf = pset.spikeFilter;
            for iA = 1:numel(adSet)
                adSet{iA} = adSet{iA}.pad([sf.preWindow sf.postWindow]);
            end
            
            pset.alignDescriptorSet = adSet;
            
            pset = pset.applyAlignDescriptorSet();
        end
        
        function pset = applyAlignDescriptorSet(pset)
            pset.warnIfNoArgOut(nargout);
            
            % align all data sources to the FIRST alignDescriptor
            % since we have to hold onto one anyway
            prog = ProgressBar(pset.nDataSources, 'Aligning data sources to each alignDescriptor');
            dataSources = pset.dataSources;
            nDataSources = pset.nDataSources;
            alignDescriptorSet = pset.alignDescriptorSet;
            
            for iSrc = 1:nDataSources
                dataSources{iSrc} = dataSources{iSrc}.align(alignDescriptorSet{:});
                prog.update(iSrc);
            end
            prog.finish();
            
            pset.dataSources = dataSources;
            
            % changing alignments invalidates everything
            pset = pset.invalidateCache();
        end
    end
    
    % manual time reslicing: converts to manual data source pset and
    % selects along time axis
    methods
        function psetManual = getAsManual(pset)
            if pset.dataSourceManual
                psetManual = pset;
            else
                psetManual = PopulationTrajectorySetBuilder.convertToManualWithSingleTrialData(pset);
            end
        end
        
        function psetSliced = manualSliceTimeWindow(pset, tMinByAlign, tMaxByAlign)
            % convert to a manual pset, then manually slice timepoints 
            % TODO adjust alignSummary data to reflect the slicing?
            assert(numel(tMinByAlign) == pset.nAlign, 'tMinByAlign must be vector with length nAlign');
            assert(numel(tMaxByAlign) == pset.nAlign, 'tMaxByAlign must be vector with length nAlign');

            pb = PopulationTrajectorySetBuilder.copyFromPopulationTrajectorySet(pset);

            % manually alter the align descriptors to truncate around zero
            for iAlign = 1:pset.nAlign
              pb.alignDescriptorSet{iAlign} = pb.alignDescriptorSet{iAlign}.windowAroundZero(...
                tMinByAlign(iAlign), tMaxByAlign(iAlign));
            end

            if pset.hasDataByTrial
              prog = ProgressBar(pset.nBases, 'Slicing dataByTrial');
              for iBasis = 1:pset.nBases
                prog.update(iBasis);
                for iAlign = 1:pset.nAlign

                  % slice data by trial
                  tvec = pset.tvecDataByTrial{iBasis, iAlign};
                  tmask = tvec >= tMinByAlign(iAlign) & tvec <= tMaxByAlign(iAlign);
                  tvecNew = tvec(tmask);
                  pb.tMinForDataByTrial(iBasis, iAlign) = tvecNew(1);
                  pb.tMaxForDataByTrial(iBasis, iAlign) = tvecNew(end);
                  pb.dataByTrial{iBasis, iAlign} = pb.dataByTrial{iBasis, iAlign}(:, tmask);

                  % also adjust the tMinByTrial and tMaxByTrial
                  pb.tMinByTrial{iBasis, iAlign} = max(tMinByAlign(iAlign), ...
                    pb.tMinByTrial{iBasis, iAlign});
                  pb.tMaxByTrial{iBasis, iAlign} = min(tMaxByAlign(iAlign), ...
                    pb.tMaxByTrial{iBasis, iAlign});
                end
              end
              prog.finish();
            end

            % adjust time stats for data mean
            for iAlign = 1:pset.nAlign
              pb.tMinValidByAlignBasisCondition(iAlign, :, :) = max(tMinByAlign(iAlign), ...
              pb.tMinValidByAlignBasisCondition(iAlign, :, :));
              pb.tMaxValidByAlignBasisCondition(iAlign, :, :) = min(tMaxByAlign(iAlign), ...
              pb.tMaxValidByAlignBasisCondition(iAlign, :, :));

              pb.tMinForDataMean(iAlign) = max(tMinByAlign(iAlign), pb.tMinForDataMean(iAlign));
              pb.tMaxForDataMean(iAlign) = min(tMaxByAlign(iAlign), pb.tMaxForDataMean(iAlign));
            end

            % slice dataByMean
            tmaskByAlign = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
              % build masks used below to slice dataMean
              tvec = pset.tvecDataMean{iAlign};
              tmaskByAlign{iAlign} = tvec >= tMinByAlign(iAlign) & tvec <= tMaxByAlign(iAlign);

              pb.dataMean{iAlign} = pb.dataMean{iAlign}(:, :, tmask);
              pb.dataSem{iAlign} = pb.dataSem{iAlign}(:, :, tmask);

              if pset.hasDataRandomized
                pb.dataMeanRandomized{iAlign} = pb.dataMeanRandomized{iAlign}(:, :, tmask, :);
                pb.dataSemRandomized{iAlign} = pb.dataSemRandomized{iAlign}(:, :, tmask, :);
              end
            end

            % slice dataDifferenceOfTrialsScaledNoiseEstimate
            tmaskCombined = cat(1, tmaskByAlign{:});
            pb.dataDifferenceOfTrialsScaledNoiseEstimate = pb.dataDifferenceOfTrialsScaledNoiseEstimate(:, tmaskCombined, :);

            psetSliced = pb.buildManualWithSingleTrialData();
        end
    end
    
    methods % Update ConditionDescriptor and AlignDescriptor appearances without invalidating everything
        function pset = setConditionAppearanceFn(pset, fn)
            % update the conditionDescriptor appearanceFn 
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is the essential part
            pset.conditionDescriptor.appearanceFn = fn;
            
            % this is for convenience / avoiding confusion
            prog = ProgressBar(pset.nDataSources, 'Updating condition appearanceFn in data sources');
            for iSrc = 1:pset.nDataSources
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.setConditionAppearanceFn(fn);
                prog.update(iSrc);
            end  
            prog.finish();
        end
        
        function pset = setInterAlignGap(pset, gaps)
            % set .interAlignGaps, which represent the time gaps between
            % successive alignments, mainly when plotting
            pset.warnIfNoArgOut(nargout);
            
            if pset.nAlign < 2
                error('Inter alignment gap not valid when only one align present');
            end
            if isscalar(gaps)
                gaps = repmat(gaps, pset.nAlign - 1, 1);
            else
                assert(numel(gaps) == pset.nAlign - 1, 'Gaps must be scalar or be length nAlign-1');
            end
            
            pset.interAlignGaps = gaps;
            
            % this is for convenience / avoiding confusion
            prog = ProgressBar(pset.nDataSources, 'Updating condition appearanceFn in data sources');
            for iSrc = 1:pset.nDataSources
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.setInterAlignGap(gaps);
                prog.update(iSrc);
            end  
            prog.finish();
        end 
        
        function pset = setStartAppearanceForAlign(pset, alignInd, varargin)
            % update this AppearanceSpec in alignDescriptorSet{ind}
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is essential
            pset.alignDescriptorSet{alignInd} = ...
                pset.alignDescriptorSet{alignInd}.setStartAppearance(varargin{:});
            
            % this is for convenience
            prog = ProgressBar(pset.nDataSources, 'Updating align appearance in each data source');
            for iSrc = 1:pset.nDataSources
                % update the align appearance but keep the active alignment 1
                % to maintain consistency
                prog.update(iSrc);
                pset.dataSources{iSrc} = ...
                    pset.dataSources{iSrc}.useAlign(alignInd).setStartAppearance(varargin{:}).useAlign(1);
            end
            prog.finish();
        end
        
        function pset = setStopAppearanceForAlign(pset, alignInd, varargin)
            % update this AppearanceSpec in alignDescriptorSet{ind}
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is essential
            pset.alignDescriptorSet{alignInd} = ...
                pset.alignDescriptorSet{alignInd}.setStopAppearance(varargin{:});
            
            % this is for convenience
            prog = ProgressBar(pset.nDataSources, 'Updating align appearance in each data source');
            for iSrc = 1:pset.nDataSources
                % update the align appearance but keep the active alignment 1
                % to maintain consistency
                prog.update(iSrc);
                pset.dataSources{iSrc} = ...
                    pset.dataSources{iSrc}.useAlign(alignInd).setStopAppearance(varargin{:}).useAlign(1);
            end
            prog.finish();
        end
        
        function pset = setZeroAppearanceForAlign(pset, alignInd, varargin)
            % update this AppearanceSpec in alignDescriptorSet{ind}
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is essential
            pset.alignDescriptorSet{alignInd} = ...
                pset.alignDescriptorSet{alignInd}.setZeroAppearance(varargin{:});
            
            % this is for convenience
            prog = ProgressBar(pset.nDataSources, 'Updating align appearance in each data source');
            for iSrc = 1:pset.nDataSources
                % update the align appearance but keep the active alignment 1
                % to maintain consistency
                prog.update(iSrc);
                pset.dataSources{iSrc} = ...
                    pset.dataSources{iSrc}.useAlign(alignInd).setZeroAppearance(varargin{:}).useAlign(1);
            end
            prog.finish();
        end
        
        function pset = setMarkAppearanceForAlign(pset, alignInd, varargin)
            % update this AppearanceSpec in alignDescriptorSet{ind}
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is essential
            pset.alignDescriptorSet{alignInd} = ...
                pset.alignDescriptorSet{alignInd}.setMarkAppearance(varargin{:});
            
            % this is for convenience
            prog = ProgressBar(pset.nDataSources, 'Updating align appearance in each data source');
            for iSrc = 1:pset.nDataSources
                % update the align appearance but keep the active alignment 1
                % to maintain consistency
                prog.update(iSrc);
                pset.dataSources{iSrc} = ...
                    pset.dataSources{iSrc}.useAlign(alignInd).setMarkAppearance(varargin{:}).useAlign(1);
            end
            prog.finish();
        end
        
        function pset = setIntervalAppearanceForAlign(pset, alignInd, varargin)
            % update this AppearanceSpec in alignDescriptorSet{ind}
            % we apply this to all bases within because it is fast, but
            % this is mostly for convenience for the user who grabs a trial
            % data from within the pset. we don't use the appearances from
            % within the trialData directly
            pset.warnIfNoArgOut(nargout);
            
            % this is essential
            pset.alignDescriptorSet{alignInd} = ...
                pset.alignDescriptorSet{alignInd}.setIntervalAppearance(varargin{:});
            
            % this is for convenience
            prog = ProgressBar(pset.nDataSources, 'Updating align appearance in each data source');
            for iSrc = 1:pset.nDataSources
                % update the align appearance but keep the active alignment 1
                % to maintain consistency
                prog.update(iSrc);
                pset.dataSources{iSrc} = ...
                    pset.dataSources{iSrc}.useAlign(alignInd).setIntervalAppearance(varargin{:}).useAlign(1);
            end
            prog.finish();
        end
    end
    
    methods % Translation / normalization
        function pset = applyTranslationNormalization(pset, trNorm)
            pset.warnIfNoArgOut(nargout);
            
            % only allow a single translation / normalization for
            % simplicity.
            if ~isempty(pset.translationNormalization)
                % Translation / normalization already applied, combine it
                % with this new one
                combined = pset.translationNormalization.combineWith(trNorm);
                pset.translationNormalization = combined;
            else
                pset.translationNormalization = trNorm;
            end
            
            cellApplyToDataFn = @(data) cellfun(@trNorm.applyTranslationNormalizationToData, ...
                data, 'UniformOutput', false);
            
            cellApplyNormOnlyToDataFn = @(data) cellfun(@trNorm.applyNormalizationToData, ...
                data, 'UniformOutput', false);
            
            if pset.dataSourceManual
                % for manual data, we apply this now to all data stored
                % manually since it will not be regenerated later
                if ~isempty(pset.dataByTrial)
                    pset.dataByTrial = trNorm.applyTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.dataMean)
                    pset.dataMean = cellApplyToDataFn(pset.dataMean);
                end
                if ~isempty(pset.dataSem)
                    pset.dataSem = cellApplyNormOnlyToDataFn(pset.dataSem);
                end
                if ~isempty(pset.dataIntervalHigh)
                    pset.dataIntervalHigh = cellApplyToDataFn(pset.dataIntervalHigh);
                end
                if ~isempty(pset.dataIntervalLow)
                    pset.dataIntervalLow = cellApplyToDataFn(pset.dataIntervalLow);
                end
                if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimate)
                    % since this is a difference between pairs of trials,
                    % we normalize here only
                    pset.dataDifferenceOfTrialsScaledNoiseEstimate = trNorm.applyNormalizationToData(pset.dataDifferenceOfTrialsScaledNoiseEstimate);
                end
            else
                % for auto-computed data, we can check whether the data has
                % been already computed (by checking the odc properties).
                % if so we can do the transformation to save time.
                % otherwise, we can defer until the buildData* methods
                % apply the translation / normalization
                %
                % It's not necessary to invalidate everything if we're careful about making
                % updates to all derived quantities here, but invalidating
                % is less prone to bugs
                
                if ~isempty(pset.odc.dataByTrial)
                    pset.dataByTrial = trNorm.applyTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.odc.dataMean)
                    pset.dataMean = cellApplyToDataFn(pset.dataMean);
                end
                if ~isempty(pset.odc.dataSem)
                    pset.dataSem = cellApplyNormOnlyToDataFn(pset.dataSem);
                end
                if ~isempty(pset.odc.dataIntervalHigh)
                    pset.dataIntervalHigh = cellApplyToDataFn(pset.dataIntervalHigh);
                end
                if ~isempty(pset.odc.dataIntervalLow)
                    pset.dataIntervalLow = cellApplyToDataFn(pset.dataIntervalLow);
                end
                if ~isempty(pset.odc.dataDifferenceOfTrialsScaledNoiseEstimate)
                    % since this is a difference between pairs of trials,
                    % we normalize here only
                    pset.dataDifferenceOfTrialsScaledNoiseEstimate = trNorm.applyNormalizationToData(pset.dataDifferenceOfTrialsScaledNoiseEstimate);
                end

            end
        end
        
        function pset = clearTranslationNormalization(pset)
             pset.warnIfNoArgOut(nargout);
            
            if isempty(pset.translationNormalization)
                % No translation / normalization applied
                return;
            end
            
            trNorm = pset.translationNormalization; 
            
            if pset.dataSourceManual
                % for manual data, we reverse this now to all data stored
                % manually since it will not be regenerated later
                if ~isempty(pset.dataByTrial)
                    pset.dataByTrial = trNorm.undoTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.dataMean)
                    pset.dataMean = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataMean, 'UniformOutput', false);
                end
                if ~isempty(pset.dataSem)
                    pset.dataSem = cellfun(@trNorm.undoNormalizationToData, pset.dataSem, 'UniformOutput', false);
                end
                if ~isempty(pset.dataIntervalHigh)
                    pset.dataIntervalHigh = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalHigh, 'UniformOutput', false);
                end
                if ~isempty(pset.dataIntervalLow)
                    pset.dataIntervalLow = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalLow, 'UniformOutput', false);
                end
                if ~isempty(pset.dataDifferenceOfTrialsScaledNoiseEstimate)
                    % since this is a difference between pairs of trials,
                    % we normalize here only
                    pset.dataDifferenceOfTrialsScaledNoiseEstimate = trNorm.undoNormalizationToData(pset.dataDifferenceOfTrialsScaledNoiseEstimate);
                end
            else
                % for auto-computed data, we can check whether the data has
                % been already computed (by checking the odc properties).
                % if so we can do the transformation to save time later.
                % otherwise, we can defer until the buildData* methods
                % apply the translation / normalization
                if ~isempty(pset.odc.dataByTrial)
                    pset.dataByTrial = trNorm.undoTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.odc.dataMean)
                    pset.dataMean = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataMean, 'UniformOutput', false);
                end
                if ~isempty(pset.odc.dataSem)
                    pset.dataSem = cellfun(@trNorm.undoNormalizationToData, pset.dataSem, 'UniformOutput', false);
                end
                if ~isempty(pset.odc.dataIntervalHigh)
                    pset.dataIntervalHigh = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalHigh, 'UniformOutput', false);
                end
                if ~isempty(pset.odc.dataIntervalLow)
                    pset.dataIntervalLow = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalLow, 'UniformOutput', false);
                end
                if ~isempty(pset.odc.dataDifferenceOfTrialsScaledNoiseEstimate)
                    % since this is a difference between pairs of trials,
                    % we normalize here only
                    pset.dataDifferenceOfTrialsScaledNoiseEstimate = trNorm.undoNormalizationToData(pset.dataDifferenceOfTrialsScaledNoiseEstimate);
                end
                
                % not necessary to invalidate if we're careful about making
                % updates to all derived quantities here.
                %pset = pset.invalidateTrialAveragedData();
            end
            
            pset.translationNormalization = [];
        end
        
        function pset = translateNormalize(pset, offsetByBasis, normalizationByBasis, varargin)
            pset.warnIfNoArgOut(nargout);
           
            if any(normalizationByBasis == 0)
                warning('Replacing normalization by 0 with 1');
                normalizationByBasis(normalizationByBasis == 0) = 1;
            end
            tr = StateSpaceTranslationNormalization.buildManual(offsetByBasis, normalizationByBasis, varargin{:});
            pset = pset.applyTranslationNormalization(tr);
        end
        
        function pset = translate(pset, offsetByBasis, varargin)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translateNormalize(offsetByBasis, onesvec(pset.nBases), varargin{:});
        end
        
        function pset = normalize(pset, normalizationByBasis, varargin)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translateNormalize(zerosvec(pset.nBases), normalizationByBasis, varargin{:});
        end
        
        function pset = meanSubtractBases(pset, varargin)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translate(-pset.computeMeanByBasis(varargin{:}), 'translationDescription', 'mean-subtracted');
        end
        
        function pset = normalizeBasesByStd(pset, varargin)
            p = inputParser();
            p.addParameter('denominatorOffset', 0, @isscalar); % x = x / (std(x) + offset)
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pset.warnIfNoArgOut(nargout);
            pset = pset.normalize(pset.computeStdByBasis(p.Unmatched) + p.Results.denominatorOffset, 'normalizationDescription', 'std-normalized');
        end
        
        function pset = zscoreByBasis(pset, varargin)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translateNormalize(-pset.computeMeanByBasis(varargin{:}), pset.computeStdByBasis(varargin{:}), ...
                'translationDescription', 'mean-subtracted', ...
                'normalizationDescription', 'std-normalized');
        end
    end
    
    % build methods for the odc properties, each must store results 
    % directly into the odc (without copying first, which allows the
    % results to persist)
    methods 
        function [src, chName] = getDataSourceForBasis(pset, iBasis)
            % return the data source which provided that data for basis
            % iBasis
            assert(isscalar(iBasis));
            src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
            chName = pset.basisDataSourceChannelNames{iBasis};
        end
        
        function idxBases = getBasisIdxForDataSource(pset, iSrc)
            assert(isscalar(iSrc));
            idxBases = find(pset.basisDataSourceIdx == iSrc);
        end
        
        function buildBasisNamesUnits(pset)
            % generate basis names by concatenating the trial data source
            % name and channel name
            
            [basisNames, basisUnits] = deal(cellvec(pset.nBases));
            for iBasis = 1:pset.nBases
                td = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                tdName = td.datasetName;
                chName = pset.basisDataSourceChannelNames{iBasis};
                
                if ~isempty(tdName)
                    basisNames{iBasis} = sprintf('%s %s', tdName, chName);
                else
                    basisNames{iBasis} = chName;
                end
                
                % this call works for unit names as well
                basisUnits{iBasis} = td.getChannelUnitsPrimary(chName);
            end
            
            % convert the basis units as specified by the
            % translationNormalization
            if ~isempty(pset.translationNormalization)
                basisUnits = pset.translationNormalization.convertBasisUnits(basisUnits);
            end
            
            % we're modifying the odc handle class here
            c = pset.odc; 
            c.basisNames = basisNames; 
            c.basisUnits = basisUnits;
        end
        
        function buildBasisValid(pset)
            % stores basisValid in odc
            % right now we do nothing except copy the manual invalid vector
            validManual = pset.basisValidManual;
            if isempty(validManual)
                validManual = truevec(pset.nBases);
            end
            
            causeManual = pset.basisInvalidCauseManual;
            if isempty(causeManual)
                causeManual = cellstrvec(pset.nBases);
            end
            hasManualCause = ~cellfun(@isempty, causeManual);
            
            % compute causes why each basis is considered invalid
            cause = cellvec(pset.nBases);
            explained = falsevec(pset.nBases);
            
            % manually marked invalid with cause
            mask = ~validManual & ~explained & hasManualCause;
            cause(mask) = causeManual(mask);
            explained(mask) = true;
            
            % manually marked invalid without cause
            mask =  ~validManual & ~explained;
            cause(mask) = {'marked invalid manually'};
            explained(mask) = true; %#ok<NASGU>
            
            valid = validManual;
            
            c = pset.odc; % this notation makes explicit that we're manipulating the internal handle
            c.basisValid = valid;
            c.basisInvalidCause = cause;
        end
        
        function buildDataByTrial(pset)
            % stores dataByTrial, alignValidByTrial, tMinByTrial, and tMaxByTrial in odc
            % This method fetches the aligned data for EVERY trial in 
            % each data source, for each alignment. It does NOT consider
            % the condition grouping, allowing this to be handled later.
            % It stores this data as a matrix, whose time vector is given
            % by the widest time vector along any trial. Missing samples in
            % this matrix are NaNs.
            
            if pset.dataSourceManual
                return;
            end
            
            [dataByTrial, tMinByTrial, tMaxByTrial, alignValidByTrial] = ...
                deal(cell(pset.nBases, pset.nAlign));
            
            [tMinForDataByTrial, tMaxForDataByTrial] = deal(nan(pset.nBases, pset.nAlign));
           
            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
%             basisAlignSummaryLookup = pset.basisDataSourceIdx;
%             alignSummaryData = cell(pset.nDataSources, pset.nAlign);
            
            dataSourcesByBasis = pset.dataSources(pset.basisDataSourceIdx);
            basisDataSourceChannelNames = pset.basisDataSourceChannelNames; %#ok<*PROP>
            nAlign = pset.nAlign;
            nBases = pset.nBases;
            timeDelta = pset.timeDelta;
            spikeFilter = pset.spikeFilter;
            
            isSpikeChannel = false(pset.nBases, 1);
            
            [dataCell, timeCell, timeInfoCell, unitsPerSecondCell] = deal(cell(pset.nBases, pset.nAlign));
            
            prog = ProgressBar(pset.nBases, 'Extracting aligned data by basis');
            for iBasis = 1:nBases    
                prog.update(iBasis);
                % request the specified aligned analog channel from the
                % specified data source.spikeFilter
                src = dataSourcesByBasis{iBasis};
                chName = basisDataSourceChannelNames{iBasis};
                
                % unapply the condition descriptor so that we can grab
                % all trials in this call, even the ones that would be
                % marked invalid by this condition info. Manually
                % invalid trials will still not be considered.
                src = src.resetConditionInfo();
                
                for iAlign = 1:nAlign
                    % mark this align as active
                    src = src.useAlign(iAlign);

                    % currently will request either analog trials or
                    % filtered spike rates channel
                    if src.hasAnalogChannel(chName)
                        [dataCell{iBasis, iAlign}, timeCell{iBasis, iAlign}] = src.getAnalog(chName);
                        isSpikeChannel(iBasis) = false;
                        
                    elseif src.hasSpikeChannel(chName)
                        %src = src.padForSpikeFilter(spikeFilter); % should have been done in applyAlignInfoSet already
                        dataCell{iBasis, iAlign} = src.getSpikeTimes(chName, true); % include padding
                        timeInfoCell{iBasis, iAlign} = src.alignInfoActive.timeInfo;
                        unitsPerSecondCell{iBasis, iAlign} = src.timeUnitsPerSecond;
                        isSpikeChannel(iBasis) = true;
                    else
                        error('Unknown channel type');
                    end
                    
                    % essential that we use computedValid here since .valid
                    % will also reflect trials which are invalid based on
                    % the current conditionInfo, which we don't want to
                    % consider here.
                    alignValidByTrial{iBasis, iAlign} = src.alignInfoActive.computedValid;
                    
                    % also store the precise time starts and stops for EACH
                    % trial that comprises that matrix. Essential that all
                    % padding be done to src before this call to ensure
                    % that the tvec returned above matches these numbers
                    [tMinByTrial{iBasis, iAlign}, ...
                     tMaxByTrial{iBasis, iAlign}] = src.getTimeStartStopEachTrial();
                end
            end
            prog.finish();
               
            prog = ProgressBar(pset.nBases, 'Resampling, filtering data by basis');
            
            for iBasis = 1:nBases
                prog.update(iBasis);
                for iAlign = 1:nAlign
                    if ~isSpikeChannel(iBasis)
                        % from TDCA.getAnalogAsMatrix
                        [dataByTrial{iBasis, iAlign}, tvec] = ...
                            TrialDataUtilies.Data.embedTimeseriesInMatrix(dataCell{iBasis, iAlign}, timeCell{iBasis, iAlign}, ...
                            'timeDelta', timeDelta, 'timeReference', 0, 'interpolate', true);
                    
                    else
                        % from TDCA.getSpikeRateFilteredAsMatrix
                        
                        % convert to .zero relative times since that's what spikeCell
                        % will be in (when called in this class)
                        timeInfo = timeInfoCell{iBasis, iAlign};
                        
                        [dataByTrial{iBasis, iAlign}, tvec] = ...
                            spikeFilter.filterSpikeTrainsWindowByTrialAsMatrix(dataCell{iBasis, iAlign}, ...
                            timeInfo.start - timeInfo.zero, timeInfo.stop - timeInfo.zero, ...
                            unitsPerSecondCell{iBasis, iAlign}, ...
                            'timeDelta', timeDelta);
                    end
                    
                    % store the time limits used in the time vector for 
                    % the dataByTrial{..} matrix
                    tMinForDataByTrial(iBasis, iAlign) = min(tvec);
                    tMaxForDataByTrial(iBasis, iAlign) = max(tvec);
                 
                    % bring trial start stop in within the limits of tvec,
                    % to deal with any rounding issues
                    tMinByTrial{iBasis, iAlign}(tMinByTrial{iBasis, iAlign} < ...
                        tMinForDataByTrial(iBasis, iAlign)) = tMinForDataByTrial(iBasis, iAlign);
                    tMaxByTrial{iBasis, iAlign}(tMaxByTrial{iBasis, iAlign} > ...
                        tMaxForDataByTrial(iBasis, iAlign)) = tMaxForDataByTrial(iBasis, iAlign);
                 
%                     assert(min(tvec) == nanmin(tMinByTrial{iBasis, iAlign}) && ...
%                            max(tvec) == nanmax(tMaxByTrial{iBasis, iAlign}), ...
%                         'Time vector returned by TrialData has invalid limits');
                end
            end
            prog.finish();            
            
            % apply translation / normalization to data
            if ~isempty(pset.translationNormalization)
                dataByTrial = pset.translationNormalization.applyTranslationNormalizationToData(dataByTrial);
            end
            
            % store the results in the odc without copying
            c = pset.odc;
%             c.alignSummaryData = alignSummaryData;
%             c.basisAlignSummaryLookup = basisAlignSummaryLookup;
            c.dataByTrial = dataByTrial;
            c.tMinForDataByTrial = tMinForDataByTrial;
            c.tMaxForDataByTrial = tMaxForDataByTrial;
            c.alignValidByTrial = alignValidByTrial;
            c.tMinByTrial = tMinByTrial;
            c.tMaxByTrial = tMaxByTrial;
        end
        
        function buildDataNTrials(pset)
            % computes and stores dataNTrials and dataValid into odc
            
            trialLists = cell(pset.nBases, pset.nConditions);
            [dataValid, dataNTrials] = deal(nan(pset.nAlign, pset.nBases, pset.nConditions));
            
            prog = ProgressBar(pset.nBases, 'Computing trial-counts by condition');
            for iBasis = 1:pset.nBases
                prog.update(iBasis); 
                % note, this src will not be aligned to this iAlign,
                % but this isn't necessary since we've already
                % extracted the aligned data
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                
                trialLists(iBasis, :) = src.conditionInfo.listByCondition(:);
                nTrialsByCondition = cellfun(@numel, src.conditionInfo.listByCondition);
                nTrialsByCondition = nTrialsByCondition(:);
                
                for iAlign = 1:pset.nAlign
                    dataNTrials(iAlign, iBasis, :) = nTrialsByCondition(:);
                    dataValid(iAlign, iBasis, :) = nTrialsByCondition > 0;
                end
            end
            prog.finish();
            
            c = pset.odc;
            c.trialLists = trialLists;
            c.dataNTrials = dataNTrials;
            c.dataValid = dataValid;
        end
        
        function buildTimeWindowsByAlignBasisCondition(pset)
            % computes and stores tMin/MaxValidByBasisAlignCondition, the
            % time windows for each basis over which enough trials exist to
            
            if pset.dataSourceManual
                return;
            end
            
            % first, compute the all-inclusive time window for each basis,
            % for each alignment, using only condition and align valid trials
            [tMinValidByAlignBasisCondition, tMaxValidByAlignBasisCondition] = ...
                deal(nan(pset.nAlign, pset.nBases, pset.nConditions));

            % do this first to force computation of data by trial at the beginning, 
            % rather than having it happen on the first loop iteration
            temp = pset.tMinByTrial; %#ok<NASGU>
            temp = pset.dataNTrials; %#ok<NASGU>
            
            prog = ProgressBar(pset.nBases, 'Computing trial-averaged time windows by basis/align/condition');
            
            for iBasis = 1:pset.nBases
                if ~pset.basisValid(iBasis);
                    continue;
                end
                for iAlign = 1:pset.nAlign
                    prog.update(iBasis);
                    
                    % note, this src will not be aligned to this iAlign,
                    % but this isn't necessary since we've already
                    % extracted the aligned data
                    src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                    tMinByTrial = pset.tMinByTrial{iBasis, iAlign};
                    tMaxByTrial = pset.tMaxByTrial{iBasis, iAlign};
                    
                    % group the condition windows by conditionLists
                    [tMinByTrialGrouped, tMaxByTrialGrouped] = src.conditionInfo.groupElements(tMinByTrial, tMaxByTrial);
                    
                    % for each basis, align, take the widest window we can that 
                    % is valid for a sufficient number of trials on this basis 
                    for iCondition = 1:pset.nConditions
                        nTrialsThis = pset.dataNTrials(iAlign, iBasis, iCondition);
                        if nTrialsThis == 0, continue, end
                        trialCountThresh = max(pset.minTrialsForTrialAveraging, ...
                            ceil(pset.minFractionTrialsForTrialAveraging*nTrialsThis));
                        if nTrialsThis >= trialCountThresh
                            tMinSorted = sort(removenan(tMinByTrialGrouped{iCondition}), 1, 'ascend');
                            tMinValidByAlignBasisCondition(iAlign, iBasis, iCondition) = tMinSorted(trialCountThresh);
                            tMaxSorted = sort(removenan(tMaxByTrialGrouped{iCondition}), 1, 'descend');
                            tMaxValidByAlignBasisCondition(iAlign, iBasis, iCondition) = tMaxSorted(trialCountThresh);
                        end
                    end
                end
            end
            prog.finish();

            c = pset.odc;
            c.tMinValidByAlignBasisCondition = tMinValidByAlignBasisCondition;
            c.tMaxValidByAlignBasisCondition = tMaxValidByAlignBasisCondition;
        end
        
        function buildDataMean(pset)
            % computes and stores dataMean, dataIntervalHigh/Low, and dataNTrials into odc
            % this function computes summary statistics across trials,
            % especially dataMean. It selects time windows which are
            % consistent across all bases and conditions for a specific
            % align to simplify subsequent operations.
            
            % warn about valid bases missing trials in at least one condition for
            % which other bases have trials. This is because the existence
            % of these bases will invalidate ALL data on those conditions
            % since there is no time window with valid trial averaged data
            % across ALL bases
            pset.warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns();
                       
            % first, we need to figure out the time windows we'll compute
            % the trial average within. For each alignment, this window will be chosen to be
            % the smallest window valid across (i.e. completely spanned by) 
            % all bases. The time window for each basis is the largest window across
            % all valid trials in all conditions, which implies that we must consider only
            % trials considered valid by the conditionInfo and alignInfo.
            % We consider all VALID trials for a given conditionDescriptor and 
            % alignDescriptor so as to not generate excessively large traces when 
            % some long trials have been filtered out. However, we do not consider
            % the particular trials which have been placed in each
            % condition, which ensures that these time windows will not
            % change if shuffling or resampling is applied.
            
            % for each align/condition, compute the widest window
            % that is valid for ALL valid bases (there's no point in taking a
            % window wider than that, for a given align/condition).
            
            % note that this creates a dependency on .basisValid for the averaged
            % data, so we invalidate the averaged data whenever the
            % basisValid is updated. We use min/max instead of
            % nanmin/nanmax here, since a NaN for a valid basis here means
            % that the basis has no trials for this condition, so the
            % condition should not be included in the trial average
%             basisMask = pset.basisValid;
%             
%             conditionsWithTrialsAllBasesAligns = squeeze(all(all(pset.dataNTrials(:, basisMask, :), 1), 2));
%             if ~any(conditionsWithTrialsAllBasesAligns)
%                 error('No conditions have trial averages for all bases on all aligns. Try lowering minTrialsForTrialAveraging or minFractionTrialsForTrialAveraging?');
%             end
%             
%             tMinValidByAlignCondition = TensorUtils.squeezeDims(max(pset.tMinValidByAlignBasisCondition(:, basisMask, conditionsWithTrialsAllBasesAligns), [], 2), 2);
%             tMaxValidByAlignCondition = TensorUtils.squeezeDims(min(pset.tMaxValidByAlignBasisCondition(:, basisMask, conditionsWithTrialsAllBasesAligns), [], 2), 2);
            
            % for each align, compute the widest window valid for ANY
            % condition for ALL bases
            tMinForDataMean = makecol(nanmin(pset.tMinValidAllBasesByAlignCondition, [], 2));
            tMaxForDataMean = makecol(nanmax(pset.tMaxValidAllBasesByAlignCondition, [], 2));
            
            if any(isnan(tMinForDataMean)) || any(isnan(tMaxForDataMean))
                error('No valid time window is valid across all bases across all conditions for some alignment');
            end
            
            nTimeByAlign = nanvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                % realign time vector to run through t=0 with increments of timeDelta
                tempVec = TrialDataUtilities.Data.linspaceIntercept(tMinForDataMean(iAlign), pset.timeDelta, tMaxForDataMean(iAlign), 0);
                tMinForDataMean(iAlign) = min(tempVec);
                tMaxForDataMean(iAlign) = max(tempVec);
            
                % number of time points for each alignment
                nTimeByAlign(iAlign) = numel(tempVec);
            end
            
            assert(pset.nBasesValid > 0, 'No valid bases found to compute trial-averaged data. Check .basisValid and .basisInvalidCause');
            assert(~isempty(tMinForDataMean) && ~isempty(tMaxForDataMean) && all(tMinForDataMean <= tMaxForDataMean), ...
                'No time window is valid across all bases. Try using .setBasesInvalidMissingTrialsOnNonEmptyConditions()');
            
            % now that we've determined the time window, we can compute the
            % trial average using data from these windows
            [dataMean, dataSem] = deal(cellvec(pset.nAlign));
            for iAlign = 1:pset.nAlign
                [dataMean{iAlign}, dataSem{iAlign}] = ...
                    deal(nan(pset.nBases, pset.nConditions, nTimeByAlign(iAlign)));
            end
           
            % copy temp values for for slicing
            conditionHasValidTrialAverageAllAlignsBases = pset.conditionHasValidTrialAverageAllAlignsBases;
            dataByTrial = pset.dataByTrial;
            tMinForDataByTrial = pset.tMinForDataByTrial;
            tMaxForDataByTrial = pset.tMaxForDataByTrial;
            timeDelta = pset.timeDelta;
            minTrialsForTrialAveraging = pset.minTrialsForTrialAveraging;
            minFractionTrialsForTrialAveraging = pset.minFractionTrialsForTrialAveraging;
            
            trialLists = pset.trialLists;
            
            for iAlign = 1:pset.nAlign
                prog = ProgressBar(pset.nBases, 'Computing trial-averaged data for align %d', iAlign);
                
                for iBasis = 1:pset.nBases
                    prog.update(iBasis);
                    % don't process invalid bases, leave these as NaNs
                    if ~pset.basisValid(iBasis);
                        continue;
                    end
                    
                    % pull the by-trial data from .dataByTrial, and use the
                    % conditioned tdca source to group the trials by
                    % condition
                    byTrial = dataByTrial{iBasis, iAlign};
                    
                    % lookup the time limits which describe the byTrial
                    % matrix
                    tMinAll = tMinForDataByTrial(iBasis, iAlign);
                    tMaxAll = tMaxForDataByTrial(iBasis, iAlign);
                    
                    if isnan(tMinAll) || isnan(tMaxAll)
                        error('Basis %d has no valid timepoints', iBasis);
                    end
                    tvecAll = tMinAll:timeDelta:tMaxAll;
                    
                    % lookup the new time limits which we'll compute the
                    % trial average within
                    tMinValid = tMinForDataMean(iAlign); 
                    tMaxValid = tMaxForDataMean(iAlign);
                    
                    % IGNORE THIS FOR NOW, CHANGED THIS ABOVE
                    % it is possible for the data by trial time vector to be
                    % smaller than the data mean time vector when some
                    % neurons are missing data on some conditions. (This
                    % will generate a warning above). For instance, if
                    % condition 1 has some bases with really long trials,
                    % and basis 2 is missing condition 1, then it might
                    % have all shorter trials, even though 
                    
                    tMaskValid = tvecAll >= tMinValid & tvecAll <= tMaxValid;
                    
                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    byTrialValid = byTrial(:, tMaskValid);
                    
                    byCondition = cellfun(@(idx) byTrialValid(idx,:), trialLists(iBasis, :)', ...
                        'UniformOutput', false);
                    
                    for iCondition = 1:pset.nConditions
                        if ~conditionHasValidTrialAverageAllAlignsBases(iCondition), continue, end
                        mat = byCondition{iCondition};
                        nTrials = size(mat, 1);
                        if nTrials == 0, continue, end;
                        
                        % minimum trial count at each time point needed to
                        % compute an average, otherwise NaN
                        minTrials = max(minTrialsForTrialAveraging, ...
                            ceil(nTrials * minFractionTrialsForTrialAveraging)); 
                        
                        nTrialsByTime = sum(~isnan(mat), 1);
                        
                        % compute mean, sem, and trial count
                        m = nanmean(mat, 1)';
                        m(nTrialsByTime < minTrials) = NaN;
                        
                        se = nansem(mat, 1)';
                        se(nTrialsByTime < minTrials) = NaN;
                        dataMean{iAlign}(iBasis, iCondition, :) = m;
                        dataSem{iAlign}(iBasis, iCondition, :) = se;
                    end
                end

                prog.finish();
            end
            
            % old way of accomplishing the same thing
%            for iAlign = 1:pset.nAlign
%                 % for this align, which conditions have at least one basis
%                 % with valid trials. nConditions x 1
%                 conditionMaskWithTrials = squeeze(any(pset.dataNTrials(iAlign, :, :), 2));
%                 
%                 % dataNtrials is nAlign x nBases x nConditions.  nBases x 1
%                 basesMissingTrials = any(pset.dataNTrials(iAlign, pset.basisValid, conditionMaskWithTrials) == 0, 3)';
%                 
%                 if any(basesMissingTrials)
%                     warning('%d bases have no valid trials for at least one condition on alignment %d for which other bases have valid trials. Use .setBasesInvalidMissingTrialAverageForNonEmptyConditionAligns() to mark these as invalid', ...
%                         nnz(basesMissingTrials), iAlign);
%                 end
%             end
            
            % no need to apply translation / normalization here since it is
            % already applied to dataByTrial!
            
            % store in odc without copying
            c = pset.odc;
            c.dataMean = dataMean;
            c.dataSem = dataSem;
            c.tMinForDataMean = tMinForDataMean;
            c.tMaxForDataMean = tMaxForDataMean;
        end

        function pset = storeDataMeanResampleTrialsWithinConditions(pset, varargin)
            % see computeDataMeanResampleTrialsWithinConditions for args,
            % this version stores the results in pset.dataMeanRandomized
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();
            
            [pset.dataMeanRandomized, pset.dataSemRandomized, pset.conditionDescriptorRandomized] = pset.computeDataMeanResampleTrialsWithinConditions(varargin{:});
        end    
        
        function [dataMeanRandomized, dataSemRandomized, cd] = computeDataMeanResampleTrialsWithinConditions(pset, varargin)
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.parse(varargin{:});
            
            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;
            
            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.resampleTrialsWithinConditions();
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);
            end
            prog.finish();
            
            [dataMeanRandomized, dataSemRandomized] = pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end
                
        function pset = storeDataMeanAxisResampleFromSpecifiedValueListIndices(pset, varargin)
            % see computeDataMeanAxisResampleFromSpecifiedValueListIndices
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace) 
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();
            
            [pset.dataMeanRandomized, pset.dataSemRandomized, pset.conditionDescriptorRandomized] = pset.computeDataMeanAxisResampleFromSpecifiedValueListIndices(varargin{:});
        end
        
        function [dataMeanRandomized, dataSemRandomized, cd] = computeDataMeanAxisResampleFromSpecifiedValueListIndices(pset, ...
                axisIdxOrAttr, resampleFromList, replace, varargin)            
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.parse(varargin{:});
            
            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;
            
            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisResampleFromSpecifiedValueListIndices(axisIdxOrAttr, resampleFromList, replace);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);
            end
            prog.finish();
            
            [dataMeanRandomized, dataSemRandomized] = pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end
        
        function pset = storeDataMeanAxisResampleFromSpecifiedValues(pset, varargin)
            % see computeDataMeanAxisResampleFromSpecifiedValuesfal
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace) 
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();
            
            [pset.dataMeanRandomized, pset.dataSemRandomized, pset.conditionDescriptorRandomized] = pset.computeDataMeanAxisResampleFromSpecifiedValues(varargin{:});
        end
        
        function [dataMeanRandomized, dataSemRandomized, cd] = computeDataMeanAxisResampleFromSpecifiedValues(pset, ...
                axisIdxOrAttr, resampleFromList, replace, varargin)            
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.parse(varargin{:});
            
            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;
            
            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisResampleFromSpecifiedValues(axisIdxOrAttr, resampleFromList, replace);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);
            end
            prog.finish();
            
            [dataMeanRandomized, dataSemRandomized] = pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function pset = storeDataMeanAxisShuffle(pset, varargin)
            % see computeDataMeanAxisShuffled
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace) 
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();
            
            [pset.dataMeanRandomized, pset.dataSemRandomized, pset.conditionDescriptorRandomized] = ...
                pset.computeDataMeanAxisShuffle(varargin{:});
        end
        
        function [dataMeanRandomized, dataSemRandomized, cd] = computeDataMeanAxisShuffle(pset, ...
                axisIdxOrAttr, varargin)        
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.parse(varargin{:});
            
            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;
            
            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisShuffle(axisIdxOrAttr);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);
            end
            prog.finish();
            
            [dataMeanRandomized, dataSemRandomized] = pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function [dataMeanRandomized, dataSemRandomized] = computeDataMeanUsingMultipleListByCondition(pset, listByConditionCell)
            % a utility function for recomputing dataMean using a different
            % conditionInfo.listByCondition than the one already present in
            % the data source. This is typically used when shuffling or resampling
            % condition labels for statistical bootstraps.
            %
            % listByConditionCell is a nBases x nConditions x nSamples
            
            tMinForDataByTrial = pset.tMinForDataByTrial;
            tMaxForDataByTrial = pset.tMaxForDataByTrial;
            tMinForDataMean = pset.tMinForDataMean;
            tMaxForDataMean = pset.tMaxForDataMean;
            timeDelta = pset.timeDelta;
            nAlign = pset.nAlign;
            nBases = pset.nBases;
            nConditions = pset.nConditions;
            nTimeDataMean = pset.nTimeDataMean;
            nRandomSamples = size(listByConditionCell, 3);
            minTrialsForTrialAveraging = pset.minTrialsForTrialAveraging;
            minFractionTrialsForTrialAveraging = pset.minFractionTrialsForTrialAveraging;
            T = sum(pset.nTimeDataMean);
            dataByTrial = pset.dataByTrial;
            
            % new matrix will contain all data, concatenated by alignment
            [dataMeanRandomizedCat, dataSemRandomizedCat] = deal(nan(nBases, nConditions, T, nRandomSamples));
            
            prog = ProgressBar(pset.nBases, 'Computing re-conditioned trial-averaged data');
            for iBasis = 1:nBases    
                
                % don't process invalid bases, leave these as NaNs
                if ~pset.basisValid(iBasis);
                    continue;
                end
                
                % compute a align-concatenated time mask for indexing into
                % data by trial
                tMaskValidByAlign = cell(nAlign, 1);
                for iAlign = 1:nAlign
                    % lookup the time limits which describe the byTrial
                    % matrix
                    tMinAll = tMinForDataByTrial(iBasis, iAlign);
                    tMaxAll = tMaxForDataByTrial(iBasis, iAlign);
                    tvecAll = (tMinAll:timeDelta:tMaxAll)';

                    % lookup the new time limits which we'll compute the
                    % trial average within
                    tMinValid = tMinForDataMean(iAlign);
                    tMaxValid = tMaxForDataMean(iAlign);
                    tMaskValidByAlign{iAlign} = tvecAll >= tMinValid & tvecAll <= tMaxValid;
                end
                tMaskValidCat = cat(1, tMaskValidByAlign{:});

                % grab the valid time portion of the nTrials x
                % nTime data matrix
                % pull the by-trial data from .dataByTrial
                byTrialCat = cat(2, dataByTrial{iBasis, :});
                byTrialCatValid = byTrialCat(:, tMaskValidCat);

                % nConditions x nSamples
                listByConditionSamples = squeeze(listByConditionCell(iBasis, :, :));
                
                % figure out how many trials for averagin we need by condition, based on
                % the number of trials each condition (assuming all lists
                % maintain the same number of trials by condition)
                nTrialsByCondition = cellfun(@numel, listByConditionSamples(:, 1));
                minTrialsByCondition = max(minTrialsForTrialAveraging, ...
                    ceil(nTrialsByCondition * minFractionTrialsForTrialAveraging));
                       
                [dataMeanRandomizedThisBasis, dataSemRandomizedThisBasis] = deal(nan(nConditions, sum(nTimeDataMean), nRandomSamples));
                for iSample = 1:nRandomSamples
                    byCondition = cellfun(@(idx) byTrialCatValid(idx,:), ...
                        listByConditionSamples(:, iSample), 'UniformOutput', false);

                    for iCondition = 1:numel(byCondition)
                        % n by t
                        mat = byCondition{iCondition};
                        nanMask = isnan(mat);
                        nTrialsByTime = sum(~nanMask, 1);
                        
                        % compute mean
                        m = nanmean(mat, 1);
                        m(nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;

                        dataMeanRandomizedThisBasis(iCondition, :, iSample) = m;

                        se = nansem(mat, 1)';
                        se(nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;
                        dataSemRandomizedThisBasis(iCondition, :, iSample) = se;
                    end
                end
                
                dataMeanRandomizedCat(iBasis, :, :, :) = dataMeanRandomizedThisBasis;
                dataSemRandomizedCat(iBasis, :, :, :) = dataSemRandomizedThisBasis;
                prog.update(iBasis);
            end
            prog.finish();

            % segregate into different alignments again
            if nRandomSamples == 1
                % avoid warning about trailing singleton dimensions
                dataMeanRandomized = squeeze(mat2cell(dataMeanRandomizedCat, nBases, nConditions, ...
                    pset.nTimeDataMean));
                dataSemRandomized = squeeze(mat2cell(dataSemRandomizedCat, nBases, nConditions, ...
                    pset.nTimeDataMean));
            else
                dataMeanRandomized = squeeze(mat2cell(dataMeanRandomizedCat, nBases, nConditions, ...
                    pset.nTimeDataMean, nRandomSamples));
                dataSemRandomized = squeeze(mat2cell(dataSemRandomizedCat, nBases, nConditions, ...
                    pset.nTimeDataMean, nRandomSamples));
            end
        end
        
        function buildDataRandomizedIntervals(pset)
            % compute the quantiles to use as intervals
            [dataIntervalHigh, dataIntervalLow] = deal(cellvec(pset.nAlign));
            qLow = pset.dataIntervalQuantileLow;
            qHigh = pset.dataIntervalQuantileHigh;
            for iAlign = 1:pset.nAlign
                dataIntervals = quantile(pset.dataMeanRandomized{iAlign}, [qLow qHigh], 4);
                dataIntervalLow{iAlign} = dataIntervals(:, :, :, 1);
                dataIntervalHigh{iAlign} = dataIntervals(:, :, :, 2);
            end

            c = pset.odc;
            c.dataIntervalHigh = dataIntervalHigh;
            c.dataIntervalLow = dataIntervalLow;  
        end
        
        function pset = withDataRandomSampleAsData(pset, dataRandomIndex, varargin)
            % take random sample idx from dataMeanRandomized as the new
            % dataMean, so that the randomization can be studied as a
            % normal pset            
            pset.warnIfNoArgOut(nargout);
            cdr = pset.conditionDescriptorRandomized;
            % the ith sample used seed initialSeed+i-1
            cdr = cdr.setRandomSeed(cdr.randomSeed+dataRandomIndex-1);
            pset = pset.setConditionDescriptor(cdr);
        end
        
        function buildAlignSummaryData(pset)
            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
            basisAlignSummaryLookup = pset.basisDataSourceIdx;
            alignSummaryData = cell(pset.nDataSources, pset.nAlign);
            
            % copy the align summary data from each data source, which may
            % take time since it is typically computed on the fly
            prog = ProgressBar(pset.nDataSources, 'Computing alignment summary statistics by data source');
            dataSources = pset.dataSources;
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc); 
                alignSummaryData(iSrc, :) = dataSources{iSrc}.alignSummarySet;
            end
            prog.finish();
            
            c = pset.odc;
            c.basisAlignSummaryLookup = basisAlignSummaryLookup;
            c.alignSummaryData = alignSummaryData;
        end

        function buildAlignSummaryAggregated(pset)
            % build the aggregated data across all bases too, for each alignment
            alignSummaryAggregated = cell(pset.nAlign, 1);
            alignSummaryData = pset.alignSummaryData;
            prog = ProgressBar(pset.nAlign, 'Computing aggregate alignment summary statistics');
            for iAlign = 1:pset.nAlign
                prog.update(iAlign);
                alignSummaryAggregated{iAlign} = AlignSummary.buildByAggregation(alignSummaryData(:, iAlign));
            end 
            prog.finish();
            
            c = pset.odc;
            c.alignSummaryAggregated = alignSummaryAggregated;
        end
        
        function buildDataNoiseEstimate(pset)
            % builds .dataDifferenceOfTrialsScaledNoiseEstimate used for
            % estimating signal vs. noise variance in StateSpaceProjections
            % This involves taking differences of pairs of random trials
            % and dividing by sqrt(2*numTrials) to give the difference
            % trace the same variance as the estimation noise in the mean
            
            % note, no need to worry about translation normalization here,
            % since this is based on dataByTrial which is already
            % normalized
            
            if ~pset.hasDataByTrial
                return;
            end
            
            % grab random pair of trials for each basis, each condition,
            % concatenate aligns in time
            NbyTAbyCby2 = pset.buildNbyTAbyCbyTrials('maxTrials', 2, 'minimizeMissingSamples', true, ...
                'message', 'Building difference of trials noise estimate');
                
            % grab trial counts
            nTrials_NbyC = pset.buildTrialCountsNbyC();
            
            % diff be N x TA x C
            dif_NbyTAbyC = diff(NbyTAbyCby2, 1, 4);

            % scale appropriately each timeseries (running along dim 2) taking into account trial counts
            nTrials_Nby1byC = permute(nTrials_NbyC, [1 3 2]);
            
            % N x TA x C
            dif_NbyTAbyC = bsxfun(@rdivide, dif_NbyTAbyC, sqrt(2*nTrials_Nby1byC));

            c = pset.odc;
            c.dataDifferenceOfTrialsScaledNoiseEstimate = dif_NbyTAbyC;
        end
    end
    
    methods % Filtering bases
        function pset = filterBases(pset, mask)
            % keep only bases listed in or selected by mask
            % only filter fields in odc if they are non-empty, implying that have already
            % been computed
            assert(isvector(mask));
            assert(islogical(mask) || all(mask >= 1 & mask <= pset.nBases));

            pset.warnIfNoArgOut(nargout);
            
            maskDim1 = {'basisNames', 'basisUnits', 'basisValid', 'basisInvalidCause', 'dataByTrial', 'tMinForDataByTrial', ...
                'tMaxForDataByTrial', 'tMinByTrial', 'tMaxByTrial', 'alignValidByTrial'};
            maskDim2 = {'dataValid', 'dataNTrials'};
            %maskCellDim1 = {'dataMean', 'dataSem', 'dataMeanRandomized', 'dataIntervalHigh', 'dataIntervalLow'}; 

            % essential that we copy before write
            c = pset.odc.copy();

            for i = 1:numel(maskDim1)
                fld = maskDim1{i};
                if ~isempty(c.(fld))
                    % some have 3 dimensions
                    c.(fld) = c.(fld)(mask, :, :);
                end
            end

            for i = 1:numel(maskDim2)
                fld = maskDim2{i};
                if ~isempty(c.(fld))
                    c.(fld) = c.(fld)(:, mask, :);
                end
            end
            
            % we need to recompute the trial-averaged data again from
            % scratch, since the time windows will be shifted to reflect
            % the bases
            c.flushTrialAveragedData();
        
%             for i = 1:numel(maskCellDim1)
%                 fld = maskCellDim1{i};
%                 if ~isempty(c.(fld))
%                     for j = 1:numel(c.(fld))
%                         % some have 4 dimensions
%                         c.(fld){j} = c.(fld){j}(mask, :, :, :);
%                     end
%                 end
%             end

            % filter alignSummaryData
            if ~isempty(c.alignSummaryData)
                alignSummaryKeep = false(size(c.alignSummaryData, 1), 1);
                alignSummaryKeep(c.basisAlignSummaryLookup(mask)) = true;

                % updating the lookup table is tricky, since the lookup
                % indices change
                newIdxForOldSummary = cumsum(alignSummaryKeep) .* alignSummaryKeep;
                newSummaryIdx = newIdxForOldSummary(pset.basisAlignSummaryLookup);
                pset.basisAlignSummaryLookup = makecol(newSummaryIdx(mask));
                
                c.alignSummaryData = c.alignSummaryData(alignSummaryKeep, :);
                c.alignSummaryAggregated = [];
            end

            % filter dataSources
            if ~isempty(pset.dataSources)
                dataSourcesKeep = false(numel(pset.dataSources), 1);
                dataSourcesKeep(pset.basisDataSourceIdx(mask)) = true;
                
                % updating the lookup table is tricky, since the lookup
                % indices change
                newIdxForOldSrc = cumsum(dataSourcesKeep) .* dataSourcesKeep;
                newSourceIdx = newIdxForOldSrc(pset.basisDataSourceIdx);
                pset.basisDataSourceIdx = makecol(newSourceIdx(mask));
                
                pset.basisDataSourceChannelNames = pset.basisDataSourceChannelNames(mask);
                
                pset.dataSources = pset.dataSources(dataSourcesKeep);
            end
            
            % filter translationNormalization
            if ~isempty(pset.translationNormalization)
                pset.translationNormalization = pset.translationNormalization.filterBases(mask);
            end
            
            if ~isempty(pset.basisValidManual)
                pset.basisValidManual = pset.basisValidManual(mask);
            end
            if ~isempty(pset.basisInvalidCauseManual)
                pset.basisInvalidCauseManual = pset.basisInvalidCauseManual(mask);
            end
            
            pset.odc = c;
                
            % force time windows to become updated
            pset = pset.updateValid();
        end
        
        function [pset, mask] = filterBasesMissingTrialAverageForNonEmptyConditionAligns(pset)  
            pset.warnIfNoArgOut(nargout);
            mask = ~pset.basesMissingTrialAverageForNonEmptyConditionAligns;
            pset = pset.filterBases(mask);
        end

        function [pset, maskC] = selectConditions(pset, varargin)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);
            
            [cd, maskC] = pset.conditionDescriptor.selectConditions(varargin{:});
            pset = pset.setConditionDescriptor(cd);
        end
        
        function [pset, maskC] = selectConditionsAlongAxis(pset, varargin)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);
            
            [cd, maskC] = pset.conditionDescriptor.selectConditionsAlongAxis(varargin{:});
            pset = pset.setConditionDescriptor(cd);
        end
        
        function [pset, maskC] = matchSelectConditionsAlongAxis(pset, varargin)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);
            
            [cd, maskC] = pset.conditionDescriptor.matchSelectConditionsAlongAxis(varargin{:});
            pset = pset.setConditionDescriptor(cd);
        end
        
%         function filterConditionsByAttribute(pset, attributeName, valueList)
%             % keep only bases listed in or masked by idx
%             p = inputParser;
%             p.addRequired('attributeName', @ischar);
%             p.addRequired('valueList', @(x) true);
%             p.parse(attributeName, valueList);
% 
%             % let the conditionDescriptor do the work
%             [pset.conditionDescriptor mask] = pset.conditionDescriptor.filteredByAttribute(attributeName, valueList, ...
%                 'removeFromGroupBy', true);
% 
%             pset.data = pset.data(:, mask, :);
%             pset.timeData = pset.timeData(:, mask, :);
%             pset.alignTimeInfoData = pset.alignTimeInfoData(:, mask, :);
%             pset.dataValid = pset.dataValid(:, mask, :);
%             pset.dataNTrials = pset.dataNTrials(:, mask, :);
% 
%             pset.tMinDataManual = pset.tMinDataManual(:, mask, :);
%             pset.tMaxDataManual = pset.tMaxDataManual(:, mask, :);
%         end
        
%         function updateConditionDescriptor(pset)
%             pset.conditionDescriptor = pset.conditionDescriptor.updateCache();
%         end
    end
    
     methods % Marking bases as invalid
        function pset = updateValid(pset)
            pset.warnIfNoArgOut(nargout);
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushValid();
                pset.odc.flushTrialAveragedData();
            end
        end
        
        function pset = resetBasisValid(pset)
            pset.warnIfNoArgOut(nargout);
            pset.basisValidManual = truevec(pset.nBases);
            pset = pset.updateValid();
        end
        
        function pset = setBasesInvalid(pset, mask, cause)
            pset.warnIfNoArgOut(nargout);
            if nargin < 3
                cause = 'unspecified';
            end
            assert(iscellstr(cause) || ischar(cause), 'Cause must be a cellstr or a string');
            
            mask = makecol(TensorUtils.vectorIndicesToMask(mask, pset.nBases));
            
            if ischar(cause)
                cause = repmat({cause}, nnz(mask), 1);
            end
            cause = makecol(cause);
            assert(numel(cause) == nnz(mask), 'Length of cellstr cause must match nnz(mask)');
            
            v = pset.basisValidManual;
            if isempty(v)
                v = truevec(pset.nBases);
            end
            v = v & ~mask;
            pset.basisValidManual = v;
            
            if isempty(pset.basisInvalidCauseManual)
                pset.basisInvalidCauseManual = cellstrvec(pset.nBases);
            end
            pset.basisInvalidCauseManual(mask) = cause;
            pset = pset.updateValid();
        end
        
        function [pset, mask] = setBasesInvalidMissingTrialAverageForNonEmptyConditionAligns(pset)  
            pset.warnIfNoArgOut(nargout);
            mask = pset.basesMissingTrialAverageForNonEmptyConditionAligns;
            pset = pset.setBasesInvalid(mask, 'missing trial average for non empty condition aligns');
        end
        
        function warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns(pset)
            N = nnz(pset.basesMissingTrialAverageForNonEmptyConditionAligns);
            if N > 0
                warning('%d bases are missing valid trial averages on condition/aligns where other bases have valid trial averages. Use .findBasesMissingTrialAverageForNonEmptyConditionAligns() to identify these and .setBasesInvalidMissingTrialAverageForNonEmptyConditionAligns() to mark these invalid', N);
            end
        end
        
        function idx = findBasesMissingTrialAverageForNonEmptyConditionAligns(pset)
            idx = find(pset.basesMissingTrialAverageForNonEmptyConditionAligns);
        end
    end
    
    methods(Static) % Static methods which take multiple psets as args 
        function [varargout] = equalizeBasesInvalid(varargin)
            assert(nargin == nargout, 'Number of input and output arguments must match');
            N = nargin;
            
            if N == 1
                varargout{1} = varargin{1};
                return;
            end
            
            nBasesVec = cellfun(@(pset) pset.nBases, varargin);
            assert(numel(unique(nBasesVec)) == 1, 'All inputs must have the same nBases');
            nBases = nBasesVec(1);
            
            mask = truevec(nBases);
            cause = cellstrvec(nBases);
            
            for i = 1:N
                pset = varargin{i};
                mask = mask & pset.basisValid;
                cause(~pset.basisValid) = pset.basisInvalidCause(~pset.basisValid);
            end
            
            varargout = cellvec(N);
            for i = 1:N
                pset = varargin{i};
                varargout{i} = pset.setBasesInvalid(~mask, cause(~mask));
            end
        end 
    end

    methods % Compute on-the-fly Dependent properties
        function tf = get.simultaneous(pset)
            tf = pset.nDataSources == 1;
        end
        
        function n = get.nTrials(pset)
            if ~pset.simultaneous
                n = NaN;
            else
                n = pset.dataSources{1}.nTrials;
            end
        end
        
        function n = get.nTrialsValid(pset)
            if ~pset.simultaneous
                n = NaN;
            else
                n = pset.dataSources{1}.nTrialsValid;
            end
        end
        
        function n = get.nDataSources(pset)
            n = numel(pset.dataSources);
        end
        
        function n = get.nBases(pset)
            if pset.dataSourceManual
                n = size(pset.dataMean{1}, 1);
            else
                n = numel(pset.basisDataSourceIdx);
            end
        end
        
        function n = get.nBasesValid(pset)
            n = nnz(pset.basisValid);
        end
        
        function n = get.nConditions(pset)
            if isempty(pset.conditionDescriptor)
                n = 0;
            else
                n = pset.conditionDescriptor.nConditions;
            end
        end
        
        function n = get.nAlign(pset)
            n = length(pset.alignDescriptorSet);
        end
  
        function sz = get.conditionsSize(pset)
            if isempty(pset.conditionDescriptor)
                sz = [0 0];
            else
                sz = pset.conditionDescriptor.conditionsSize;
            end
        end
        
        function sz = get.conditionsSizeNoExpand(pset)
            if isempty(pset.conditionDescriptor)
                sz = 0;
            else
                sz = pset.conditionDescriptor.conditionsSizeNoExpand;
            end
        end
        
        function n = get.nAxes(pset)
            if isempty(pset.conditionDescriptor)
                n = 0;
            else
                n = pset.conditionDescriptor.nAxes;
            end
        end     
        
        function as = lookupAlignSummaryDataForBasis(pset, basisIdx)
            asIdx = pset.basisAlignSummaryLookup(basisIdx);
            as = pset.alignSummaryData(asIdx, :);
        end
        
        function as = lookupAlignSummaryDataForBasisAlign(pset, basisIdx, alignIdx)
            asIdx = pset.basisAlignSummaryLookup(basisIdx);
            as = pset.alignSummaryData{asIdx, alignIdx};
        end
        
        function n = get.nAlignSummaryData(pset)
            n = size(pset.alignSummaryData, 1);
        end
        
        function conditions = get.conditions(pset)
            if isempty(pset.conditionDescriptor)
                conditions = {};
            else
                conditions = pset.conditionDescriptor.conditions;
            end
        end 
   
        function names = get.conditionNames(pset)
            if isempty(pset.conditionDescriptor)
                names = {};
            else
                names = pset.conditionDescriptor.names;
            end
        end

        function names = get.alignNames(pset)
            names = cellfun(@(ad) ad.name, pset.alignDescriptorSet, 'Uniform', false);
        end
        
        function str = get.dataIntervalQuantilesAsString(pset)
            str = sprintf('[%g - %g] %%', 100*pset.dataIntervalQuantileLow, ...
                100*pset.dataIntervalQuantileHigh);
        end
        
        function c = get.alignConditionsWithTrials(pset)
            % dataNTrials is nAlign x nBases x nConditions
            % c is nAlign x nConditions
            c = TensorUtils.squeezeDims(any(pset.dataNTrials(:, pset.basisValid, :), 2), 2);
        end
        
        function c = get.alignConditionsWithTrialAverage(pset)
            % c is nAlign x nConditions
            c = TensorUtils.squeezeDims(any(pset.alignBasisConditionsWithValidTrialAverage(:, pset.basisValid, :), 2), 2);
        end
        
        function c = get.alignBasisConditionsWithValidTrialAverage(pset)
            % c is nAlign x nBases x nConditions
            c = ~isnan(pset.tMinValidByAlignBasisCondition) & ...
                ~isnan(pset.tMaxValidByAlignBasisCondition);
        end
        
        function c = get.basisConditionsWithValidTrialAverage(pset)
            % c is nBases x nConditions
            c = TensorUtils.squeezeDims(all(pset.alignBasisConditionsWithValidTrialAverage, 1), 1);
        end
%         
        function basesMissingTrialAverages = get.basesMissingTrialAverageForNonEmptyConditionAligns(pset)
            % determine which bases missing trials in at least one condition for
            % which other bases have trials in that alignment. basesMask is
            % nBases x 1 logical vector
            
            % nAlign x nBases x nConditions
            shouldHaveTrialAverage = repmat(permute(pset.alignConditionsWithTrialAverage, [1 3 2]), [1 pset.nBasesValid 1]);
            hasTrialAverage = pset.alignBasisConditionsWithValidTrialAverage(:, pset.basisValid, :);
            
            % nBases x 1
            basesValidMissingTrialAverages = makecol(TensorUtils.squeezeDims(any(any(shouldHaveTrialAverage & ~hasTrialAverage, 3), 1), [1 3]));
            basesMissingTrialAverages = TensorUtils.inflateMaskedTensor(basesValidMissingTrialAverages, 1, pset.basisValid, false);
        end
        
        function bMask = get.basesNonEmpty(pset)
            % nBases x 1 logical
            bMask = pset.basisValid & makecol(any(any(pset.alignBasisConditionsWithValidTrialAverage, 1), 3));
        end
        
        function cMask = get.conditionsWithValidTrialAverageOnNonEmptyBases(pset)
            % look at valid bases that have a trial average on at least one
            % condition. which conditions have trial averages on ALL bases
            % satisfying this?
            bMask = pset.basesNonEmpty;
            % any alignment, all bases
            cMask = squeeze(all(any(pset.alignBasisConditionsWithValidTrialAverage(:,bMask,:), 1), 2));
        end
           
    end
    
    methods % Resampling, shuffling : WILL IMPLEMENT PASS-THRUS to
%     condition descriptor
%
%         function shuffleSourcesAlong(pset, compareAlong)
%             assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
%             pset.dataSources = cellfun(@(sr) sr.buildShuffledAlong(compareAlong), pset.dataSourcesOrig, 'UniformOutput', false); 
%             pset.updateFromDataSources();
%         end
% 
%         function resampleSources(pset)
%             assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
%             pset.dataSources = cellfun(@(sr) sr.buildResampled(), pset.dataSourcesOrig, 'UniformOutput', false); 
%             pset.updateFromDataSources();
%         end
% 
%         function resampleSourcesFromSingleAttributeValue(pset, attr, value)
%             assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
%             pset.dataSources = cellfun(@(sr) sr.buildResampledFromSingleAttributeValue(attr, value), pset.dataSourcesOrig, 'UniformOutput', false); 
%             pset.updateFromDataSources();
%         end
    end

    methods % Simple statistics, should be NaN for invalid bases
        function snrByBasis = computeSnrByBasis(pset, varargin)
            semCTAbyN = pset.buildCTAbyN('type', 'sem');
            noiseByBasis = nanmax(semCTAbyN, [], 1)';
            rangeByBasis = pset.computeRangeByBasis();
            
            snrByBasis = rangeByBasis ./ noiseByBasis;
            snrByBasis(~pset.basisValid) = NaN;
        end
        
        function minByBasis = computeMinByBasis(pset, varargin)
            CTAbyN = pset.buildCTAbyN(varargin{:});
            minByBasis = nanmin(CTAbyN, [], 1)';
            minByBasis(~pset.basisValid) = NaN;
        end
        
        function maxByBasis = computeMaxByBasis(pset, varargin)
            CTAbyN = pset.buildCTAbyN(varargin{:});
            maxByBasis = nanmax(CTAbyN, [], 1)';
            maxByBasis(~pset.basisValid) = NaN;
        end
        
        function meanByBasis = computeMeanByBasis(pset, varargin)
            CTAbyN = pset.buildCTAbyN(varargin{:});
            meanByBasis = nanmean(CTAbyN, 1)';
            meanByBasis(~pset.basisValid) = NaN;
        end
        
        function varByBasis = computeVarByBasis(pset, varargin)
            varByBasis = nanvar(pset.buildCTAbyN(varargin{:}), 0, 1)';
            varByBasis(~pset.basisValid) = NaN;
        end
        
        function stdByBasis = computeStdByBasis(pset, varargin)
            stdByBasis = nanstd(pset.buildCTAbyN(varargin{:}), 0, 1)';
            stdByBasis(~pset.basisValid) = NaN;
        end
        
        function rangeByBasis = computeRangeByBasis(pset, varargin)
            CTAbyN = pset.buildCTAbyN(varargin{:});
            rangeByBasis = nanmax(CTAbyN, [], 1)' - nanmin(CTAbyN, [], 1)';
            rangeByBasis(~pset.basisValid) = NaN;
        end
  
%         function [maxValues maxTimes] = getMaximumByBasisEachConditionEachAlign(pset, varargin)
%             % maxValues: nAlign x nBases x nConditions matrices with the maximum value
%             % time of occurrence
% 
%             maxValues = nan(pset.nAlign, pset.nBases, pset.nConditions);
%             for iAlign = 1:pset.nAlign 
%                 [maxValues(iAlign, :, :), maxTimes(iAlign, :, :)] = ...
%                     max();
%             end
% 
%             function [valMax timeMax] = findMaxFn(time,data)
%                 [valMax i] = max(data);
%                 timeMax = time(i);
%             end
%         end
% 
%         function [minValues minTimes] = findMinimum(pset, varargin)
%             % return nBases x nConditions x nAlign matrices with the minimum value and
%             % time of minimum value for each basis x condition x align
% 
%             [minValues minTimes] = pset.alignBasisConditionDataFun(@findMinFn, 'asMatrix', [true true]);
% 
%             function [valMin timeMin] = findMinFn(time,data)
%                 [valMin i] = min(data);
%                 timeMin = time(i);
%             end
%         end
% 
%         function maxByBasis = findMaximumAcrossConditions(pset)
%             % return a nBases x nAlign vector of maximum values across conditions
%             
%             maxValues = pset.findMaximum();
%             maxByBasis = max(maxValues, [], 2);
%         end
% 
%         function minByBasis = findMinimumAcrossConditions(pset)
%             % return a nBases x nAlign vector of minimum values across conditions
%             
%             minValues = pset.findMinimum();
%             minByBasis = min(minValues, [], 2);
%         end
%         
%         function varData = computeVarianceOverTime(pset)
%             % varData is nBases x nConditions x nAlign
%             
%             varData = pset.alignBasisConditionDataFun(@var, 'asMatrix', true);
%         end
% 
%         function [ccvByBasisByAlign tVecByBasisByAlign] = crossConditionVariance(pset)
%             % compute variance across conditions per condition over time for each basis x align
%             % ccvByBasisByAlign is a cell array of nBases x nAlign with a T x 1 vector of ccv over time
%             % tVecByAlign has the same size and carries the time vector associated with each ccv vector
% 
%             [data timeData tVecByBasisByAlign] = pset.getDataTimeWindowedValidAcrossConditions();
% 
%             ccvByBasisByAlign = TensorUtils.mapToSizeFromSubs([pset.nBases pset.nAlign], ...
%                 'contentsFn', @getCCV, 'asCell', true);
% 
%             function ccv = getCCV(iBasis, iAlign)
%                 % time by conditions matrix for this basis x align
%                 tByC = cell2mat(squeezedim(data(iBasis, :, iAlign), [1 3]));
% 
%                 % t by 1 vector of variance across conditions
%                 ccv = makecol(var(tByC, [], 2));
%             end
%         end
%         
%         function [ccvByBasisByAlign tVecByBasisByAlign] = crossConditionStd(pset)
%             % compute variance across conditions per condition over time for each basis x align
%             % ccvByBasisByAlign is a cell array of nBases x nAlign with a T x 1 vector of ccv over time
%             % tVecByAlign has the same size and carries the time vector associated with each ccv vector
% 
%             [data timeData tVecByBasisByAlign] = pset.getDataTimeWindowedValidAcrossConditions();
% 
%             ccvByBasisByAlign = TensorUtils.mapToSizeFromSubs([pset.nBases pset.nAlign], ...
%                 'contentsFn', @getCCS, 'asCell', true);
% 
%             function ccv = getCCS(iBasis, iAlign)
%                 % time by conditions matrix for this basis x align
%                 tByC = cell2mat(squeezedim(data(iBasis, :, iAlign), [1 3]));
% 
%                 % t by 1 vector of variance across conditions
%                 ccv = makecol(std(tByC, [], 2));
%             end
%         end
    end

    methods % Visualization and plotting utilities
%         function alignTimeInfo = drawTimeAxisForAlign(pset, iAlign, alignTimeInfo)
%             ad = pset.alignDescriptorSet{iAlign};
%             if ~exist('alignTimeInfo', 'var')
%                 alignTimeInfo = pset.alignTimeInfoData(:, :, iAlign);
%                 emptyMask = cellfun(@isempty, alignTimeInfo);
%                 alignTimeInfo = alignTimeInfo(~emptyMask);
%                 alignTimeInfo = cell2mat(makecol(alignTimeInfo(:)));
%             end
%             %ad.drawTimeAxis(alignTimeInfo);
%         end
%         
%         function drawTimeAxisForConditionAlign(pset, iCondition, iAlign)
%             ad = pset.alignDescriptorSet{iAlign};
%             alignTimeInfo = pset.alignTimeInfoData(:, iCondition, iAlign);
%             alignTimeInfo = cell2mat(makecol(alignTimeInfo(:)));
%             ad.drawTimeAxis(alignTimeInfo);
%         end

%         function plotConditionPanels(pset, varargin)
%             % draw all bases together, with each condition as a panel
%             %fig();
%             clf
%             p = panel();
%             p.pack(nRow, nCol);
%             p.margin = 10;
%             
%             % determine how to layout the conditions, use a row if 1-d conditions
%             % 2-d if necessary
%             nDims = min(2, pset.conditionDescriptor.nAxes);
%             conditionInds = pset.conditionDescriptor.conditionsAsLinearInds;
%             
%             if pset.ConditionDescriptor.nAxes == 1
%                 nRows = 1;
%                 nCols = 
%                 
%                 
%             for iAlign = 1:pset.nAlign 
%                 if nDims == 1
%                     % if there is only one attribute, plot it along one row
%                     nRow = 1;
% 
%                     % we need a column for every valid condition on this align
%                     nCol = nnz(conditionsValidThisAlign);
%                     
%                     conditionByRowCol = makerow(conditionInds(conditionsValidThisAlign));
%                 else
%                     % one row for each value of the first attribute
%                     nRow = pset.conditionDescriptor.nValuesByAttributeGroupBy(1);
%                     
%                     % reshape conditionsValid into nRow rows and the rest
%                     % as columns
%                     conditionIndsReshaped = reshape(conditionInds(:), nRow, []);
%                     conditionsValidReshaped = reshape(conditionsValidThisAlign, nRow, []);
%                     
%                     % pare down the rows and columns that have at least one
%                     % valid condition in them
%                     rowMask = any(conditionsValidReshaped, 2);
%                     colMask = any(conditionsValidReshaped, 1);
%                     conditionByRowCol = conditionIndsReshaped(rowMask, colMask);
%                     
%                     nRow = nnz(rowMask);
%                     nCol = nnz(colMask);
%                 end
%                 
%                 
% 
%                 % build a nice colormap
%                 cmap = cbrewer('qual', 'Set1', pset.nBases);
%                 cmap = jet(pset.nBases);
%                 
%                 % loop over row and column panels
%                 yl = nan(pset.nConditions, 2);
%                 panelHasData = false(nRow, nCol);
%                 for iCol = 1:nCol
%                     for iRow = 1:nRow
%                         
%                         iCondition = conditionByRowCol(iRow, iCol);
%                         if ~conditionsValidThisAlign(iCondition)
%                             continue;
%                         end
%                         
%                         h(iCondition) = p(iRow, iCol).select();
% 
%                         [tMin, tMax, yMin, yMax] = deal(NaN);
%                        
%                         for iBasis = 1:pset.nBases
%                             dataVec = data{iBasis, iCondition, iAlign};
%                             tvec = time{iBasis, iCondition, iAlign};
%                             
%                             if ~isempty(dataVec) && any(~isnan(dataVec))
%                                 panelHasData(iRow, iCol) = true;
%                             end
%                             
%                             tMin = min([tMin min(tvec)]);
%                             tMax = max([tMax max(tvec)]);
%                             yMin = min([yMin min(dataVec)]);
%                             yMax = max([yMax max(dataVec)]);
%                             
%                             plot(tvec, dataVec, '-', 'LineWidth', 2, 'Color', cmap(iBasis,:));
%                             hold on
%                             
%                             yl(iCondition, :) = get(gca, 'YLim');
%                         end
%                         
%                         xlim([tMin tMax]);
%                         ylim([yMin yMax]);
%                         
%                         title(sprintf('%s (%s)', pset.conditionNames{iCondition}, pset.alignNames{iAlign}));
%                     end
%                 end
%                 
%                 % draw time axes
%                 for iCol = 1:nCol
%                     for iRow = 1:nRow
%                         p(iRow, iCol).select();
%                         iCondition = conditionByRowCol(iRow, iCol);
%                         if ~panelHasData(iRow, iCol)
%                             axis off;
%                         else
%                             pset.drawTimeAxisForConditionAlign(iCondition, iAlign);
%                         end
%                     end
%                 end
%                 %whitebg(gcf, [0 0 0]);
%                 p.refresh();
%             end
%         end
% 
%         function plotBasisPanels(pset, varargin)
%             p = inputParser;
%             p.addParameter('basisIdx', [1:6], @(x) isvector(x) && ...
%                 all(inRange(x, [1 pset.nBases])));
%             p.parse(varargin{:});
% 
%             basisIdx = intersect(p.Results.basisIdx, 1:pset.nBases);
%             nBasesPlot = length(basisIdx);
% 
%             [data, time] = pset.getDataTimeWindowed();
% 
%             for iAlign = 1:pset.nAlign
%                 fig();
%                 clf;
%                 p = panel();
%                 p.pack(nBasesPlot,1);
% 
%                 for iBasisIdx = 1:nBasesPlot
%                     p(iBasisIdx,1).select();
%                     iBasis = basisIdx(iBasisIdx);
% 
%                     for iCondition = 1:pset.nConditions
%                         timeVec = time{iBasis, iCondition, iAlign};
%                         dataVec = data{iBasis, iCondition, iAlign};
% 
%                         appear = pset.conditionDescriptor.appearances(iCondition);
% 
%                         plot(timeVec, dataVec, ...
%                             'Color', appear.color, 'LineWidth', appear.lineWidth);
%                         hold on
%                     end
%                     hold off
%                     box off
%                     title(sprintf('%s (%s)', pset.basisNames{iBasis}, pset.alignNames{iAlign}));
%                     pset.drawTimeAxisForAlign(iAlign);
%                     drawnow;
%                 end
%                 
%                 p.margin = 10;
%             end
%         end

        function addConditionLabels(pset, varargin)
            pset.conditionDescriptor.addColoredLabels(varargin{:});
        end

        function plotBases(pset, varargin)
            % plot bases one above the next 
            p = inputParser;
            p.addParameter('dataRandomIndex', [], @(x) isempty(x) || isscalar(x));
            
            p.addParameter('alignIdx', 1:pset.nAlign, @isnumeric);            
            p.addParameter('basisIdx', [], @(x) isempty(x) || (isvector(x) && ...
                all(inRange(x, [1 pset.nBases]))));
            p.addParameter('validBasesOnly', true, @islogical);
            p.addParameter('conditionIdx', 1:pset.nConditions, @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nConditions])));
            
            % change the y size of each basis to fit a band of height
            % 'scaling' below. if false, neither normalize nor scaling do
            % anything. this is best used when a single basis is being
            % plotted
            p.addParameter('scaleBases', true, @islogical);
            
            % plot each basis at it's original scale or normalized to fit
            % the bands
            p.addParameter('normalize', true, @islogical);
            % each basis starts at y = 1, 2, 3 and data is scaled to fit a
            % band with y-height scaling
            p.addParameter('scaling', 0.8, @isscalar);
            
            p.addParameter('showVerticalScaleBars', true, @islogical);
            p.addParameter('showBasisLabels', true, @islogical);
            
            % usually plot first basis at top, last basis at bottom
            p.addParameter('reverse', false, @islogical);
            
            % when no interAlignGap is specified, this is used to compute
            % the interAlignGap
            p.addParameter('alignGapFraction', 0.02, @isscalar);
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('yOffset', 0, @isscalar);
            
            p.addParameter('lineWidth', 1, @isscalar);
            
            p.addParameter('alpha', 1, @isscalar); % alpha for main traces
            p.addParameter('showSem', false, @islogical); % show standard error traces?
            p.addParameter('errorAlpha', 0.5, @isscalar); % alpha for surrounding error fills
            
            % show less by default than in TDCA
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('markShowOnAxis', true, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 2, @isscalar);
            
            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalShowOnAxis', true, @islogical);
            p.addParameter('intervalAlpha', 0.5, @isscalar);
            
            p.addParameter('showRangesOnData', false, @islogical); % show ranges for marks on traces
            p.addParameter('showRangesOnAxis', true, @islogical); % show ranges for marks below axis
            
            p.addParameter('timeAxisStyle', 'marker', @ischar); % 'tickBridge' or 'marker'
            
            % make room for labels using AutoAxis
            p.addParameter('axisMarginLeft', 2.5, @isscalar);
           
            p.parse(varargin{:});

            axh = newplot;
            wasHolding = ishold(axh);
            hold(axh, 'on');
            
            if p.Results.validBasesOnly
                if isempty(p.Results.basisIdx)
                    basisIdx = find(pset.basisValid, 20, 'first');
                else
                    basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
                    eachValid = pset.basisValid(basisIdx);
                    basisIdx = basisIdx(eachValid);
                end
            else
                if isempty(p.Results.basisIdx);
                    basisIdx = 1:min(pset.nBases, 20);
                else
                    basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
                end
            end
            nBasesPlot = numel(basisIdx);
            conditionIdx = p.Results.conditionIdx;
            nConditionsPlot = numel(conditionIdx);
            
            scaleBases = p.Results.scaleBases;
            scaling = p.Results.scaling;
            normalize = p.Results.normalize;
            reverse = p.Results.reverse;
            alignGapFraction = p.Results.alignGapFraction;
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            
            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            timeWidthByAlign = pset.nTimeDataMean(alignIdx)*pset.timeDelta;
            nAlignUsed = numel(alignIdx);
            
            if isempty(pset.interAlignGaps)
                % compute absolute x-gap between alignments
                alignGaps = repmat(alignGapFraction*sum(timeWidthByAlign) / (1 - alignGapFraction*nAlignUsed), nAlignUsed-1, 1);
            else
                alignGaps = pset.interAlignGaps;
            end 
            
            % keep track of start and stop of each align in time
            tAlignZero = nanvec(nAlignUsed);

            % concatenate the data across alignments in order to compute
            % normalization constants. data is N x CTA
            
            if ~isempty(p.Results.dataRandomIndex)
                % use a sample from dataMeanRandomized instead
                data = pset.buildCTAbyN('type', 'meanRandom', 'dataRandomIndex', p.Results.dataRandomIndex, ...
                    'basisIdx', basisIdx, 'conditionIdx', conditionIdx, 'alignIdx', alignIdx)';
            else
                data = pset.buildCTAbyN('basisIdx', basisIdx, ...
                    'conditionIdx', conditionIdx, 'alignIdx', alignIdx)';
            end
            
            if p.Results.showSem
                if ~isempty(p.Results.dataRandomIndex)
                    % use a sample from dataSemRandomized instead
                    dataSem = pset.buildCTAbyN('type', 'semRandom', 'dataRandomIndex', p.Results.dataRandomIndex, ...
                        'basisIdx', basisIdx, 'conditionIdx', conditionIdx, 'alignIdx', alignIdx)';
                else
                    dataSem = pset.buildCTAbyN('type', 'sem', 'basisIdx', basisIdx, ...
                        'conditionIdx', conditionIdx, 'alignIdx', alignIdx)';
                end
            end
            
            if p.Results.scaleBases
                if normalize
                    % each basis will be independently scaled to [0 1]
                    if p.Results.showSem
                        % include sem in limits
                        offsets = nanmin(data - dataSem, [], 2);
                        norms = nanmax(data + dataSem, [], 2) - offsets;
                    else
                        offsets = nanmin(data, [], 2); % N x 1
                        norms = nanmax(data, [], 2) - offsets;
                    end
                else
                    % data will collectively be scaled to [0 1], but the same
                    % transformation will apply to all bases 
                    if p.Results.showSem
                        % include sem in limits
                        offsets = nanmin(data - dataSem, [], 2); % N x 1

                        %offsets = repmat(m, nBasesPlot, 1);
                        ranges = nanmax(data + dataSem, [], 2) - nanmin(data - dataSem, [], 2);
                        norms = repmat(nanmax(ranges), nBasesPlot, 1);
                    else
                        offsets = nanmin(data, [], 2); % N x 1

                        %offsets = repmat(m, nBasesPlot, 1);
                        ranges = nanmax(data, [], 2) - nanmin(data, [], 2);
                        norms = repmat(nanmax(ranges), nBasesPlot, 1);
                    end
                end
                            
                % apply scaling directly to norms since norms are computed
                % for unit scaling above
                norms = norms / scaling;
            else
                offsets = zerosvec(nBasesPlot);
                norms = onesvec(nBasesPlot);
            end
            
            app = pset.conditionDescriptor.appearances;
            
            hData = cell(nConditionsPlot, nAlignUsed);
            for iAlign = 1:nAlignUsed
                idxAlign = alignIdx(iAlign);
                tvecOrig = pset.tvecDataMean{idxAlign};
                % figure out where zero in tvec should be located
                if iAlign == 1
                    tOffsetCurrent = nanmin(tvecOrig);
                end
                tAlignZero(iAlign) = tOffsetCurrent - nanmin(tvecOrig) + xOffset;
                
                tvecPlot = tvecOrig + tAlignZero(iAlign);
                
                % here we assemble all data across bases, normalize and
                % scale and offset them vertically. This makes plotting
                % much more efficient than looping over bases individually.
                
                % N x C x T matrices
                % apply separate translation /normalization to each basis
                % to bring in range
                if isempty(p.Results.dataRandomIndex)
                    data = pset.dataMean{idxAlign}(basisIdx, conditionIdx, :);
                else
                    data = pset.dataMeanRandomized{idxAlign}(basisIdx, conditionIdx, :, p.Results.dataRandomIndex);
                end
                data = bsxfun(@minus, data, offsets);
                data = bsxfun(@rdivide, data, norms);
                
                if isempty(p.Results.dataRandomIndex)
                    dataSem = pset.dataSem{idxAlign}(basisIdx, conditionIdx, :);
                else
                    dataSem = pset.dataSemRandomized{idxAlign}(basisIdx, conditionIdx, :, p.Results.dataRandomIndex);
                end
                % only apply normalization
                dataSem = bsxfun(@rdivide, dataSem, norms);
                
                % uniformly scale and separate data vertically
                data = data + yOffset;
                data = bsxfun(@plus, data, (nBasesPlot-1:-1:0)');
                
                if reverse
                    data = flipud(data);
                    dataSem = flipud(dataSem);
                end
                
                % draw error bars
                if p.Results.showSem
                    for iCondition = 1:nConditionsPlot
                        idxCondition = conditionIdx(iCondition);
                        dataC = TensorUtils.squeezeDims(data(:, iCondition, :), 2);
                        dataSemC = TensorUtils.squeezeDims(dataSem(:, iCondition, :), 2);
                        for iBasis = 1:nBasesPlot
                            hShade = TrialDataUtilities.Plotting.errorshade(tvecPlot, dataC(iBasis, :), ...                   
                                dataSemC(iBasis, :), app(idxCondition).Color, 'axh', axh, ...
                                'alpha', p.Results.errorAlpha, 'z', 1, 'showLine', false);
                            TrialDataUtilities.Plotting.hideInLegend(hShade);
                        end
                    end
                end
                
                % draw means
                for iCond = 1:nConditionsPlot
                    idxCondition = conditionIdx(iCond);
                    dataC = TensorUtils.squeezeDims(data(:, iCond, :), 2);
                    if p.Results.alpha < 1
                        hData{iCond, iAlign} = TrialDataUtilities.Plotting.patchline(tvecPlot, dataC', ...
                            'EdgeColor', app(idxCondition).Color, 'EdgeAlpha', p.Results.alpha, ...
                            'LineWidth', p.Results.lineWidth, 'Parent', axh);
                    else
                        hData{iCond, iAlign} = plot(tvecPlot, dataC, '-', ...
                            'Parent', axh, 'LineWidth', p.Results.lineWidth, ...
                            'Color', app(idxCondition).Color);
                    end
                end
                
                %draw marks and intervals on data
                if p.Results.markShowOnData || p.Results.intervalShowOnData
                    for iBasis = 1:numel(basisIdx)
                        idxBasis = basisIdx(iBasis);
                        alignSummary = pset.alignSummaryData{pset.basisAlignSummaryLookup(idxBasis), idxAlign};
                        
                        % data needs to be T x D x C, currently N (select 1) x C x T 
                        dataPermute = permute(data(iBasis, :, :), [3 1 2]);
                        alignSummary.drawOnDataByCondition(tvecOrig, ...
                            dataPermute + yOffset, ...  
                            'conditionIdx', conditionIdx, ...
                            'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
                            'tOffsetZero', tAlignZero(iAlign), 'alpha', p.Results.alpha, ...
                            'markAlpha', p.Results.markAlpha, 'markSize', p.Results.markSize, ...
                            'intervalAlpha', p.Results.intervalAlpha, ...
                            'showRanges', p.Results.showRangesOnData, ...
                            'tMin', min(tvecOrig), 'tMax', max(tvecOrig));
                    end
                end
                
                % advance laterally
                if iAlign == nAlignUsed
                    gap = 0;
                else
                    gap = alignGaps(iAlign);
                end
                tOffsetCurrent = tOffsetCurrent + timeWidthByAlign(iAlign) + gap;
            end
            
            box off;
            if scaleBases
                ylim([-0.1, nBasesPlot+0.1]);
            end
            xlim([min(pset.tvecDataMean{alignIdx(1)}) tOffsetCurrent]);
            
            % setup auto axis
            au = AutoAxis();
            
            if p.Results.showBasisLabels
                yloc = yOffset + (nBasesPlot-1:-1:0)' + 0.5;
                ylabel = pset.basisNames(basisIdx);
                au.addTicklessLabels('y', 'tick', yloc, 'tickLabel', ylabel);   
            end
            au.yUnits = pset.dataUnits;
            
            switch p.Results.timeAxisStyle
                case 'tickBridge'
                    for iAlign = 1:nAlignUsed
                        idxAlign = alignIdx(iAlign);
                        if numel(basisIdx) == 1
                            as = pset.alignSummaryData{pset.basisDataSourceIdx(basisIdx), idxAlign};
                        else
                            as = pset.alignSummaryAggregated{idxAlign};
                        end
                        tvec = pset.tvecDataMean{idxAlign};
                        offset = tAlignZero(iAlign);
                        
                        as.setupTimeAutoAxis('axh', gca, 'style', 'tickBridge', ...
                            'tMin', min(tvec), 'tMax', max(tvec), 'tOffsetZero', offset, ...
                            'showRanges', p.Results.showRangesOnAxis);
                    end
                    
                case 'marker'
                    for iAlign = 1:nAlignUsed
                        idxAlign = alignIdx(iAlign);
                        if numel(basisIdx) == 1 && ~isempty(pset.basisDataSourceIdx)
                            as = pset.alignSummaryData{pset.basisDataSourceIdx(basisIdx), idxAlign};
                        else
                            as = pset.alignSummaryAggregated{idxAlign};
                        end
                        tvec = pset.tvecDataMean{idxAlign};
                        offset = tAlignZero(iAlign);
                        
                        as.setupTimeAutoAxis('axh', gca, 'style', 'marker', ...
                            'tMin', min(tvec), 'tMax', max(tvec), 'tOffsetZero', offset, ...
                            'showRanges', p.Results.showRangesOnAxis);
                    end
                    au.xUnits = pset.timeUnitName;
                    au.addAutoScaleBarX();
            end
            
            % add scale bars to right side of axis
            au.yUnits = pset.dataUnits;
            if p.Results.showVerticalScaleBars
                if ~scaleBases
                    au.addAutoScaleBar('y');
                elseif ~normalize
                    % all data scaled by same amount, show one scale bar at
                    % bottom whose length is roughly half the dynamic range of
                    % half the range of the highest-amplitude channel.
                    % all the norms are the same, such that 1 y-unit
                    % corresponds to norms(1).
                    plottedValue = scaling;
                    actualValue = TrialDataUtilities.Plotting.closestNiceNumber(norms(1) * plottedValue, 'down');
                    label = sprintf('%g %s', actualValue, pset.dataUnits);
                    au.addScaleBar('y', 'length', plottedValue, 'manualLabel', label);
                else
                    % all data are scaled independently, draw n scale bars
                    % with the same size. the actual data size of each will
                    % be set by the largest amplitude signal, which has the
                    % largest norm
                    
                    actualValue = TrialDataUtilities.Plotting.closestNiceNumber(min(norms), 'down');
                    for iBasis = 1:numel(basisIdx)
                        plottedValue = actualValue / norms(iBasis);
                        if iBasis == numel(basisIdx)
                            label = sprintf('%g %s', actualValue, pset.dataUnits);
                        else
                            label = '';
                        end
                        if ~isnan(offsets(iBasis)) && ~isnan(plottedValue)
                            au.addScaleBar('y', 'length', plottedValue, 'manualLabel', label, ...
                                'manualPositionAlongAxis', numel(basisIdx) - iBasis);
                        end
                    end
                    
                end
            end
                    
            % make large left side to accommodate labels
            au.axisMarginLeft = p.Results.axisMarginLeft;
            axis off;
            au.update();
            au.installCallbacks();
            
            if ~wasHolding
                hold(axh, 'off');
            end
            
            set(axh, 'SortMethod', 'childorder');
        end

        function plotSingleBasis(pset, basisIdx, varargin)
            assert(isscalar(basisIdx));
            pset.plotBases('basisIdx', basisIdx, 'showSem', true, ...
                'showVerticalScaleBars', false, 'showBasisLabels', false, 'scaleBases', false, varargin{:});
            au = AutoAxis(gca);
            au.ylabel('spikes/s');
            au.addAutoAxisY();
            au.update();
        end
        
        function plotStateSpace(pset, varargin)
            % plot a 2d or 3d basis1 x basis2 x basis3 trajectory plot
            p = inputParser;
            p.addParameter('basisIdx', 1:min(pset.nBases, 3), @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', [], @isnumeric);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('lineWidth', 2, @isscalar);
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 10, @isscalar);
            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('useTranslucentMark3d', true, @islogical);
            p.addParameter('markShowOnData', true, @islogical);
            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);
            
            p.addParameter('showRangesOnData', true, @islogical); % show ranges for marks on traces
            p.parse(varargin{:});
            
            axh = newplot;
            hold(axh, 'on');
            
            basisIdx = p.Results.basisIdx;
            conditionIdx = p.Results.conditionIdx;
            if islogical(conditionIdx)
                conditionIdx = find(conditionIdx);
            end
            
            % figure out which alignments are used
            if isempty(p.Results.alignIdx)
                alignIdx = 1:pset.nAlign;
            else
                alignIdx = 1:pset.nAlign;
                alignIdx = alignIdx(p.Results.alignIdx);
            end
            nAlignUsed = numel(alignIdx);
            
            nConditions = numel(conditionIdx);
            nBasesPlot = numel(basisIdx);
            plotArgs = p.Results.plotOptions;

            switch(nBasesPlot)
                case 2
                    use3d = false;
                case 3
                    use3d = true;
                otherwise
                    error('Number of bases must be 2 or 3. For 1 basis, use plotBases(''basisIdx'', idx)');
            end

            for iAlign = 1:nAlignUsed
                idxAlign = alignIdx(iAlign);
                tvec = pset.tvecDataMean{idxAlign};
                data = pset.dataMean{idxAlign};

                for iCondition = 1:nConditions
                    c = conditionIdx(iCondition);
                    appear = pset.conditionDescriptor.appearances(c);
                    plotArgsC = appear.getPlotArgs();
                    
                    dataMat = squeeze(data(basisIdx, c, :));
                    
                    if p.Results.alpha < 1
                        if use3d
                            TrialDataUtilities.Plotting.patchline(dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                'EdgeColor', appear.Color, 'EdgeAlpha', p.Results.alpha, ...
                                 plotArgs{:});
                        else
                            TrialDataUtilities.Plotting.patchline(dataMat(1, :), dataMat(2, :), ...
                                'LineWidth', p.Results.lineWidth, 'Parent', axh, ...
                                'EdgeColor', appear.Color, 'EdgeAlpha', p.Results.alpha, ...
                                plotArgs{:});
                        end
                    else
                        if use3d
                            plot3(axh, dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                'LineWidth', p.Results.lineWidth, plotArgsC{:}, plotArgs{:});
                        else
%                             dataMat = [dataVec1 dataVec2];
                            plot(axh, dataMat(1, :), dataMat(2, :), ...
                                'LineWidth', p.Results.lineWidth, ...
                                plotArgsC{:}, plotArgs{:});
                        end
                    end
                end
            
                % draw marks and intervals on the data traces
                as = pset.alignSummaryAggregated{idxAlign};
                % data is nBases x C x T; drawOnDataByConditions needs T x nBasesPlot x C
                dataForDraw = permute(data(basisIdx, conditionIdx, :), [3 1 2]);
                as.drawOnDataByCondition(tvec, dataForDraw, ...
                    'conditionIdx', conditionIdx, 'markAlpha', p.Results.markAlpha, ...
                    'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
                    'useTranslucentMark3d', p.Results.useTranslucentMark3d, ...
                    'alpha', p.Results.alpha, ...
                    'intervalAlpha', p.Results.intervalAlpha, ...
                    'showRanges', p.Results.showRangesOnData, ...
                    'markSize', p.Results.markSize); 
            end

            box(axh, 'off')
            axis(axh, 'tight');

             xlabel(pset.basisNames{basisIdx(1)});
             ylabel(pset.basisNames{basisIdx(2)});
             if use3d
                zlabel(pset.basisNames{basisIdx(3)});
                view(axh, [-40 20]);
                axis(axh, 'vis3d');
                axis(axh, 'off');
                ThreeVector(axh);
             else
                axis(axh, 'off');
                ThreeVector(axh);
             end
             
             hold(axh, 'off');
             
        end
        
    end

    methods % Build data matrices
        % Notation:
        % C is nConditions
        % T is nTimepoints for a given alignment, such that T*A really
        %     means T(align1) + T(align2) + T(align3) + ...
        % A is nAlign
        % N is nBases
        %
        % all of these take a few common arguments:
        % conditionIdx, basisIdx, alignIdx - pick which ones you want
        %   included
        
        function [NbyCbyTA, tvec, avec] = buildNbyCbyTA(pset, varargin)
            % [NbyCbyTA, tvec, avec] = buildNbyCbyTA(pset, ...)
            p = inputParser();
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('type', 'mean', @ischar); % mean, sem, meanRandom
            p.addParameter('dataRandomIndex', 1, @isscalar);
            p.parse(varargin{:});
            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            nAlign = numel(alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            conditionIdx = makecol(p.Results.conditionIdx);
            
            if p.Results.validBasesOnly
                mask = pset.basisValid(basisIdx);
                basisIdx = basisIdx(mask);
            end
            
            data = cellvec(nAlign);
            for iAlign = 1:numel(alignIdx)
                idxAlign = alignIdx(iAlign);
                switch p.Results.type
                    case 'mean'
                        data{iAlign} = pset.dataMean{idxAlign}(basisIdx, conditionIdx, :);
                    case 'sem'
                        data{iAlign} = pset.dataSem{idxAlign}(basisIdx, conditionIdx, :);
                    case 'meanRandom'
                        assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');
                        assert(p.Results.dataRandomIndex >= 1 && p.Results.dataRandomIndex <= pset.nRandomSamples, ...
                            'dataRandomIndex must be in range 1:%d (nRandomSamples)', pset.nRandomSamples);
                        data{iAlign} = pset.dataMeanRandomized{idxAlign}(basisIdx, conditionIdx, :, p.Results.dataRandomIndex);
                    case 'semRandom'
                        assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');
                        assert(p.Results.dataRandomIndex >= 1 && p.Results.dataRandomIndex <= pset.nRandomSamples, ...
                            'dataRandomIndex must be in range 1:%d (nRandomSamples)', pset.nRandomSamples);
                        data{iAlign} = pset.dataSemRandomized{idxAlign}(basisIdx, conditionIdx, :, p.Results.dataRandomIndex);
                    otherwise
                        error('Unknown data type %s', p.Results.type);
                end
            end 
            
            [NbyCbyTA, avecRaw] = TensorUtils.catWhich(3, data{:});
            avec = makecol(alignIdx(avecRaw));
            tvec = cat(1, pset.tvecDataMean{alignIdx});
        end

        function [CTAbyN, cvec, tvec, avec, nvec] = buildCTAbyN(pset, varargin)
            % out is C*T*A x N concatenated matrix
            p = inputParser();
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('type', 'mean', @ischar);
            p.addParameter('dataRandomIndex', 1, @isscalar);
            p.parse(varargin{:});
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            conditionIdx = TensorUtils.vectorMaskToIndices(p.Results.conditionIdx);

            [NbyCbyTA, tvec, avec] = pset.buildNbyCbyTA(p.Results);
            nvec = basisIdx;
            cvec = conditionIdx;
            labels = {nvec, cvec, [tvec, avec]};
            [CTAbyN, labelsOut] = TensorUtils.reshapeByConcatenatingDims(NbyCbyTA, {[3 2], 1}, labels);

            cvec = labelsOut{1}(:, 1);
            tvec = labelsOut{1}(:, 2);
            avec = labelsOut{1}(:, 3);
            nvec = labelsOut{2};
        end

         function [NbyTAbyAttr, tvec, avec] = buildNbyTAbyConditionAttributes(pset, varargin)
             % build a tensor of N by T by nValsAttr1 by nValsAttr2 x ...
             % this tensor is the format used by dpca_covs
             p = inputParser();
             p.addParameter('type', 'mean', @ischar);
             p.KeepUnmatched = true;
             p.parse(varargin{:});
             
             [NbyCbyTA, tvec, avec] = pset.buildNbyCbyTA('type', p.Results.type, p.Unmatched);
             NbyTAbyC = permute(NbyCbyTA, [1 3 2]);
             N = size(NbyTAbyC, 1);
             TA = size(NbyTAbyC, 2);
             condSize = pset.conditionDescriptor.conditionsSize;
             NbyTAbyAttr = reshape(NbyTAbyC, [N TA makerow(condSize)]);
         end
         
         % added primarily for dpca-type noise floor
         function [NbyTAbyCbyR, nTrials_NbyC, tvec, avec, whichTrials] = buildNbyTAbyCbyTrials(pset, varargin)
            p = inputParser();
            p.addParameter('maxTrials', Inf, @isscalar);
            p.addParameter('chooseRandom', false, @islogical);
            p.addParameter('minimizeMissingSamples', false, @islogical);
            p.addParameter('randomSeed', 0, @isscalar);
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('message', 'Extracting individual trial data matrix by basis', @ischar);
            p.parse(varargin{:});
            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            nAlign = numel(alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            conditionIdx = makecol(p.Results.conditionIdx);
            nConditions = numel(conditionIdx);
            
            if p.Results.validBasesOnly
                mask = pset.basisValid(basisIdx);
                basisIdx = basisIdx(mask);
            end
            
            nBases = numel(basisIdx);
            
            % figure out how many trials to allocate for
            dataNTrials = pset.dataNTrials(alignIdx, basisIdx, conditionIdx);
            maxTrials = min(p.Results.maxTrials, max(dataNTrials(:)));
            
            if p.Results.chooseRandom
                s = RandStream('mt19937ar', 'Seed', p.Results.randomSeed);
            end
            
            % ultimately want N x TA x C x R, so we build up an N x T x C x R
            % tensor for each align and then cat along 2
            dataByAlign = cellvec(nAlign);
            nTrials_NbyC = nan(nBases, nConditions);
            whichTrials = cell(nBases, nConditions);
            
            % pre-allocate per-align
            for iAlignIdx = 1:nAlign
                iAlign = alignIdx(iAlignIdx);
                dataByAlign{iAlignIdx} = nan(nBases, pset.nTimeDataMean(iAlign), nConditions, maxTrials);
            end
                
            prog = ProgressBar(nBases, p.Results.message);
            for iBasisIdx = 1:nBases
                prog.update(iBasisIdx);
                iBasis = basisIdx(iBasisIdx);
                
                if ~pset.basisValid(iBasis)
                    continue;
                end
                
                fullListByCondition = pset.trialLists(iBasis, conditionIdx)';
                nTrialsByCondition = cellfun(@numel, fullListByCondition);
                nTrialsByCondition = min(nTrialsByCondition, maxTrials);

                % pick trials to include for this basis
                if p.Results.minimizeMissingSamples
                    % try to pick trials that look similar to the
                    % trial average in terms of where missing
                    % samples are located, over all aligns
                    
                    % collect the single trial data matrix and data mean
                    % across aligns for this basis
                    byTrialOverAligns = cellvec(nAlign);
                    dataMeanOverAligns = cellvec(nAlign);
                    for iAlignIdx = 1:nAlign
                        iAlign = alignIdx(iAlignIdx);
                        tvecThis = pset.tvecDataByTrial{iBasis,iAlign};
                        tMaskValidThis = tvecThis >= pset.tMinForDataMean(iAlign) & tvecThis <= pset.tMaxForDataMean(iAlign);
                        % grab the valid time portion of the nTrials x
                        % nTime data matrix
                        byTrialOverAligns{iAlign} = pset.dataByTrial{iBasis, iAlign}(:, tMaskValidThis);
                        dataMeanOverAligns{iAlign} = TensorUtils.squeezeDims(pset.dataMean{iAlign}(iBasis, :, :), 1);
                    end
                    
                    listByCondition = cellvec(nConditions);
                    for iC = 1:nConditions
                        % count number of valid timepoints in dataMean that are
                        % present in each trial, this will be the ranking order
                        % for the trials, we prefer this number to be higher.
                        matchingSamplesOverAligns = cellfun(@(byTrial, dataMean) sum(bsxfun(@and, ...
                            ~isnan(byTrial(fullListByCondition{iC}, :)), ~isnan(dataMean(iC, :))), 2), byTrialOverAligns, dataMeanOverAligns, ...
                            'UniformOutput', false);
                        totalMatchingSamples = sum(cat(2, matchingSamplesOverAligns{:}), 2);
                        
                        [~, sortIdx] = sort(totalMatchingSamples, 1, 'descend');
                        listByCondition{iC} = fullListByCondition{iC}(sortIdx(1:nTrialsByCondition(iC)))';
                    end

                elseif p.Results.chooseRandom
                    listByCondition = cellfun(@(list, n) randsample(s, list, n), ...
                        fullListByCondition, num2cell(nTrialsByCondition), ...
                        'UniformOutput', false);
                else
                    listByCondition = cellfun(@(list, n) list(1:n), ...
                        fullListByCondition, num2cell(nTrialsByCondition), ...
                        'UniformOutput', false);
                end    
                
                % store list for this by basis in matrix
                whichTrials(iBasisIdx, :) = listByCondition(:);
                    
                for iAlignIdx = 1:nAlign
                    iAlign = alignIdx(iAlignIdx);
                    
                    % compute time mask to select from this data by trial
                    % to match dataMean's tvec
                    tvecThis = pset.tvecDataByTrial{iBasis,iAlign};
                    tMaskValid = tvecThis >= pset.tMinForDataMean(iAlign) & tvecThis <= pset.tMaxForDataMean(iAlign);
                    
                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    byTrialValid = pset.dataByTrial{iBasis, iAlign}(:, tMaskValid);

                    % and reuse this list on subsequent aligns
                    listByCondition = whichTrials(iBasisIdx, :);
                    
                    for iConditionIdx = 1:nConditions
                        dataByAlign{iAlignIdx}(iBasisIdx, :, iConditionIdx, 1:nTrialsByCondition(iConditionIdx)) = ...
                            permute(byTrialValid(listByCondition{iConditionIdx}, :), [3 2 4 1]);
                        if iAlignIdx == 1
                            nTrials_NbyC(iBasisIdx, iConditionIdx) = nTrialsByCondition(iConditionIdx);
                        end
                    end
                 end
            end
            prog.finish();
            
            [NbyTAbyCbyR, avecRaw] = TensorUtils.catWhich(2, dataByAlign{:});
            avec = makecol(alignIdx(avecRaw));
            tvec = cat(1, pset.tvecDataMean{alignIdx});
         end
         
         % added primarily for dpca-type noise floor
         function [NbyTAbyAttrbyR, nTrials_NbyAttr, tvec, avec, whichTrials_NbyAttr] = buildNbyTAbyConditionAttributesbyTrials(pset, varargin)
             [NbyTAbyCbyR, nTrialsTensor, tvec, avec, whichTrials] = pset.buildNbyTAbyCbyTrials(varargin{:});
             N = size(NbyTAbyCbyR, 1);
             TA = size(NbyTAbyCbyR, 2);
             R = size(NbyTAbyCbyR, 4);
             condSize = makerow(pset.conditionDescriptor.conditionsSize);
             if pset.conditionDescriptor.nAxes == 1
                 condSize = condSize(1);
             end
             
             assert(prod(condSize) == size(NbyTAbyCbyR, 3), 'Sub-selecting conditionIdx not supported');
             
             NbyTAbyAttrbyR = reshape(NbyTAbyCbyR, [N TA condSize R]);
             
             % N x C --> N by condSize
             nTrials_NbyAttr = reshape(nTrialsTensor, [N condSize]);
             whichTrials_NbyAttr = reshape(whichTrials, [N condSize]);
         end
         
         function [NbyTAbyRbyAttr, nTrials_NbyAttr, tvec, avec, whichTrials_NbyAttr] = buildNbyTAbyTrialsbyConditionAttributes(pset, varargin)
             [NbyTAbyCbyR, nTrialsTensor, tvec, avec, whichTrials] = pset.buildNbyTAbyCbyTrials(varargin{:});
             N = size(NbyTAbyCbyR, 1);
             TA = size(NbyTAbyCbyR, 2);
             R = size(NbyTAbyCbyR, 4);
             NbyTAbyRbyC = permute(NbyTAbyCbyR, [1 2 4 3]);
             condSize = pset.conditionDescriptor.conditionsSize;
             
             assert(prod(condSize) == size(NbyTAbyCbyR, 3), 'Sub-selecting conditionIdx not supported');
             
             NbyTAbyRbyAttr = reshape(NbyTAbyRbyC, [N TA R makerow(condSize)]);
             
             % N x C --> N by condSize
             nTrials_NbyAttr = reshape(nTrialsTensor, [N makerow(condSize)]);
             whichTrials_NbyAttr = reshape(whichTrials, [N makerow(condSize)]);
         end

         function nTrials_NbyC = buildTrialCountsNbyC(pset)
             nTrials_NbyC = permute(min(pset.dataNTrials, [], 1), [2 3 1]);
         end
         
         function nTrials_NbyAttr = buildTrialCountsNbyConditionAttr(pset)
             nTrials_NbyC = pset.buildTrialCountsNbyC();
             condSize = pset.conditionDescriptor.conditionsSize;
             nTrials_NbyAttr = reshape(nTrials_NbyC, [size(nTrials_NbyC, 1) makerow(condSize)]);
         end
         
%         function [out tvecByConditionAlign] = buildCTByNEachA(pset, varargin)
%             % out: A cell vector of CT x N matrices 
%             % timeVec: nAlign cell of time vectors common to alignment
%             
%             p = inputParser();
%             p.addParameter('timeValidAcrossConditions', false, @islogical)
%             p.parse(varargin{:});
% 
%             % data is N bases x nConditions x nAlign cell array with vectors of length t
%             % allow different time vectors per condition
%             if p.Results.timeValidAcrossConditions
%                 [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
%                 tvecByConditionAlign = repmat(timeVecByAlign', pset.nConditions, 1);
%             else
%                 [data, timeData, tvecByConditionAlign] = pset.getDataTimeWindowedValidAcrossBases();
%             end
% 
%             % convert to nAlign cell of nCondition*time x nBases 
%             out = cell(pset.nAlign, 1);
%             for iAlign = 1:pset.nAlign
%                 out{iAlign} = cell2mat(data(:,:,iAlign)');
%             end
%         end
% 
%         function [out tvecByConditionAlign] = buildTByNEachCA(pset, varargin)
%             % out: C x A cell of N x T matrices 
%             % tvecByConditionAlign is C x A cell of time vectors common to condition x alignment
% 
%             p = inputParser();
%             p.addParameter('timeValidAcrossConditions', false, @islogical)
%             p.parse(varargin{:});
%                         
%             % data is N bases x nConditions x nAlign cell array with vectors of length t
%             if p.Results.timeValidAcrossConditions
%                 [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
%                 tvecByConditionAlign = repmat(makerow(timeVecByAlign), pset.nConditions, 1);
%             else
%                 [data, timeData, tvecByConditionAlign] = pset.getDataTimeWindowedValidAcrossBases();
%             end
% 
%             out = cell(pset.nConditions, pset.nAlign);
%             for iAlign = 1:pset.nAlign
%                 for iCondition = 1:pset.nConditions
% 
%                     % nBases cell of time vectors
%                     dataThisCA = data(:,iCondition, iAlign);
% 
%                     % nAlign*nTime x nBases matrix for this condition
%                     out{iCondition, iAlign} = cell2mat(dataThisCA');
%                 end
%             end
%         end
%         
    end

    methods % Single trial
        function assertSingleTrialData(pset)
            assert(pset.simultaneous && ~pset.dataSourceManual, 'This method may only be called on simultaneously collected data sets, i.e. where there is only one data source');
        end
        
        function [NbyTAbyR, tvec, avec, cvec] = buildNbyTAbyTrials(pset, varargin)
            p = inputParser();
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('validTrialsOnly', false, @islogical);
            p.parse(varargin{:});
            
            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            if p.Results.validBasesOnly
                mask = pset.basisValid(basisIdx);
                basisIdx = basisIdx(mask);
            end
            conditionIdx =  TensorUtils.vectorMaskToIndices(p.Results.conditionIdx);
            
            % dataByTrial is nBases x nAlign, insides are R x T_a
            % concatenate time over alignments
            [tvec, avec] = TensorUtils.catWhich(1, pset.tvecDataByTrial{1, :});
            RbyTA_eachBasis = TensorUtils.catInnerDimOverOuterDim(pset.dataByTrial(basisIdx, alignIdx), 2, 2);
            RbyTAbyN = cat(3, RbyTA_eachBasis{:});
            NbyTAbyR = permute(RbyTAbyN, [3 2 1]);
            
            cvec = pset.dataSources{1}.conditionIdx;
            rmask = ismember(cvec, conditionIdx);
            
            if p.Results.validTrialsOnly
                rmask = rmask & pset.dataSources{1}.valid;
            end
            
            cvec = cvec(rmask);
            NbyTAbyR = NbyTAbyR(:, :, rmask);  
        end
    end
    
    methods % Comparative statistics
        function [distByFromAlign, timeVecByFromAlign] = getDistanceBetween(pset, cFromList, cToList, varargin)
            % dist: length(fromAlign) cell of T x length(cFromList) distance traces as columns
            % comparing each condition cFromList(i) to condition cToList(i)

            p = inputParser;
            % if true, search the entire from trajectory for the closest point to 
            % each point on the to trajectory looking across multiple alignments as well.
            % If false, compute distance at each timepoint separately
            p.addParameter('searchEntire', true, @islogical); 
            % calculate distances within this subset of bases
            p.addParameter('basisIdx', true(pset.nBases, 1), @isvector);

            % when searchEntire is true, calculate distances from trajectories (with condition cTo) along
            % each of the alignments indexed in fromAlign, to the closest point in
            % trajectories in ANY of the alignments indexed in  toAlign
            % to the closest point on condition cTo within these alignments
            p.addParameter('fromAlign', 1:pset.nAlign, @isvector);
            p.addParameter('toAlign', 1:pset.nAlign, @isvector);
            
            % leave empty to compute a distance vs time trajectory for the
            % entire from trajectory at each alignment. Populate with a
            % nAlign x 1 cell array of time points to compute the distance
            % only from a specific set of time points
            p.addParameter('timepointsByFromAlign', [], @(x) isempty(x) || iscell(x));

            p.addParameter('showPlot', true, @islogical);

            p.parse(varargin{:});
            searchEntire = p.Results.searchEntire;
            basisIdx = p.Results.basisIdx;
            fromAlign = p.Results.fromAlign;
            toAlign = p.Results.toAlign;
            showPlot = p.Results.showPlot;
            timepointsByFromAlign = p.Results.timepointsByFromAlign;
            
            % scalar expansion to match cFrom and cTo
            if isscalar(cFromList)
                cFromList = repmat(cFromList, size(cToList));
            elseif isscalar(cToList)
                cToList = repmat(cToList, size(cFromList));
            end

            assert(isequal(size(cFromList), size(cToList)), 'Sizes of cFrom and cTo must match');
            nPair = length(cFromList);

            % C x A cell array of T x N data points
            [tByNEachCA, timeVecByAlign] = pset.buildTByNEachCA('timeValidAcrossConditions', true);
            timeVecByFromAlign = makecol(timeVecByAlign(1, fromAlign));

            % loop over pieces of the cFrom trajectories from each alignment
            distByFromAlign = cell(length(fromAlign), 1);
            for iIdxFromAlign = 1:length(fromAlign)
                iFromAlign = fromAlign(iIdxFromAlign);
                tvecFromAlign = timeVecByFromAlign{iFromAlign};
                nTime = length(tvecFromAlign);

                % are we computing distances only from selected timepoints on the from trajectories?
                if ~isempty(timepointsByFromAlign)
                    timepoints = timepointsByFromAlign{iFromAlign};
                    nTime = length(timepoints);
                    timeMask = nan(nTime, 1); 
                    for iTime = 1:nTime
                        ind = find(floor(tvecFromAlign) == floor(timepoints(iTime)), 1, 'first');
                        assert(~isempty(ind), 'Could not find timepoint %g in alignment %s', ...
                            timepoints(iTime),pset.alignNames{iFromAlign});
                        timeMask(iTime) = ind;
                    end
                    timeVecByFromAlign{iFromAlign} = timeVecByFromAlign{iFromAlign}(timeMask);
                else
                    timeMask = true(size(tvecFromAlign)); 
                end

                distByFromAlign{iFromAlign} = nan(nTime, nPair);

                % loop over comparison condition pairs (maybe avoidable?)
                for iPair = 1:nPair
                    cFrom = cFromList(iPair);
                    cTo = cToList(iPair);

                    % dataFrom is T*A x N matrix of concatenated trajectories for condition cFrom
                    dataFrom = tByNEachCA{cFrom, iFromAlign};
                    dataFrom = dataFrom(timeMask,basisIdx);

                    if searchEntire
                        % we search over every point in the cTo trajectories in all
                        % alignments in toAlign

                        % dataTo is an T*A x N matrix of concatenated trajectories for condition cTo
                        dataTo = cell2mat(tByNEachCA(cTo, toAlign)');
                        dataTo = dataTo(:, basisIdx);

                        % dist will be T*A x 1 vector of distances from
                        % each point along dataFrom to ANY point along
                        % dataTo
                        dist = pdist2(dataTo, dataFrom, 'euclidean', 'Smallest', 1);
                    else
                        % dataTo is T*A x N matrix of concatenated trajectories for condition cTo
                        dataTo = tByNEachCA{cTo, iFromAlign};
                        dataTo = dataTo(timeMask, basisIdx);

                        % compute distances and summing over N bases
                        dist = sqrt(sum((dataFrom - dataTo).^2, 2));
                    end
                    
                    distByFromAlign{iIdxFromAlign}(:, iPair) =  dist;
                end
            end

            if showPlot
                % display a 1 x nAlign row of all condition pairs distance traces superimposed
                fig();
                p = panel();
                p.pack(1, length(fromAlign));

                legstr = cell(nPair, 1);

                for iIdxFromAlign = 1:length(fromAlign)
                    p(1, iIdxFromAlign).select();
                    tvec = timeVecByFromAlign{iIdxFromAlign};

                    for iPair = 1:nPair
                        % use the appearance for the from condition, arbitrarily
                        appear = pset.conditionDescriptor.appearances(cFromList(iPair));

                        plot(tvec, distByFromAlign{iIdxFromAlign}(:, iPair), ...
                            'k-', 'LineWidth', 2, ...
                            'Color', appear.color, 'LineWidth', appear.lineWidth);
                        hold on

                        if iIdxFromAlign == 1
                            legstr{iPair} = sprintf('%s to %s', pset.conditionNames{cFromList(iPair)}, ...
                                pset.conditionNames{cToList(iPair)});
                        end
                    end
                    box off
                    hold off
                    xlim([min(tvec) max(tvec)]);
                end

                figsize(6, 6*length(fromAlign));
                for iIdxFromAlign = 1:length(fromAlign)
                    p(1, iIdxFromAlign).select();

                    if iIdxFromAlign == 1
                        legend(legstr, 'Location', 'NorthEast', 'FontSize', 10);
                        legend boxoff;
                    end
                    pset.drawTimeAxisForAlign(fromAlign(iIdxFromAlign));
                end
            end
        end

        function [distByFromAlign, timeVecByFromAlign, cFromList, cToList] = getDistanceAlongComparisonAxis(pset, compareAcross, varargin)
            p = inputParser;
            % default is to compare
            % compare from condition 1 to condition 2 instead of 2 to 1 along the axis?
            p.addRequired('compareAcross', @ischar);
            p.addParameter('reverse', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(compareAcross, varargin{:});
            reverse = p.Results.reverse;

            idxCompare = pset.conditionDescriptor.compareAlong(compareAcross);
            
            nCompare = length(idxCompare); 
            [cFromList, cToList] = deal(zeros(nCompare, 1));
            for iCompare = 1:nCompare
                idxThisComparison = idxCompare{iCompare};
                assert(length(idxThisComparison) == 2, 'Comparison axis must span exactly two elements');
                if ~reverse
                    idxThisComparison = idxThisComparison([2 1]);
                end
                cFromList(iCompare) = idxThisComparison(1);
                cToList(iCompare) = idxThisComparison(2);
            end

            [distByFromAlign, timeVecByFromAlign] = ...
                pset.getDistanceBetween(cFromList, cToList, p.Unmatched);
        end
    end

end
