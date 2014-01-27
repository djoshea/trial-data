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
    end

    properties(SetAccess=protected, Hidden, Transient)
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
        
        % SpikeFilter instance to use when converting spiking units to
        % firing rate channels
        spikeFilter
        
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
        
        % When multiple alignDescriptors are used to align the data, some
        % trials may be valid only for some alignments and not others.
        % Setting this flag to true ensures that only trials which are
        % valid for ALL alignments will be considered 
        includeOnlyTrialsValidAllAlignments
        
        % number of randomized samples to draw when generating
        % dataMeanRandomized
        nRandomSamples
        
        % random seed to use as initial seed when generating random data
        % sets
        randomSeed

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

        % ConditionDescriptor instance describing condition information 
        conditionDescriptor
        
        % StateSpaceTranslationNormalization instance describing the
        % translation and normalization to apply to each basis.
        % This will be applied during buildDataByTrial for dataSourceManual
        % = false psets (thus being reflected in the trial-averages automatically)
        % or applied manually to dataByTrial (if non-empty) and dataMean
        % (if non-empty) for dataSourceManual = false
        translationNormalization
    end

    % These properties store raw data sources (TrialDataConditionAlign instances)
    % from which data is extracted, as well as track from where each basis
    % originates.
    properties(SetAccess=?PopulationTrajectorySetBuilder)
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
    end
    
    % Properties whose values are computed dynamically and persist within odc
    % or are specified manually and persist within manualData
    properties(Hidden, Dependent, Transient, SetAccess=?PopulationTrajectorySetBuilder)
        % Alignment summary statistics by basis by align
        
        % nAlignSummary x nAlign cell containing AlignSummary instances for each
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
 
        %% BELOW ARE FOR TRIAL-AVERAGED DATA WITHIN CONDITION
        
        % trial averaged data contains common time vectors across
        % conditions and bases, but which differ across aligns. Invalid
        % portions of the bases with smaller windows will be filled with
        % NaNs.
        
        % nAlign x 1 numeric vectors indicating the time window used for each
        % alignmnet for the dataMean and dataIntervalHigh/Low cells
        tMinForDataMean
        tMaxForDataMean
        
        % nAlign x 1 numeric vectors indicating the number of timepoints
        % used for each alignment in the dataMean cell
        nTimeDataMean
        
        % nAlign x 1 cell array of time vectors corresponding to the 3rd
        % dimension of dataMean{iAlign}
        tvecDataMean
        
        % trial-averaged data within each condition
        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign) 
        % numeric array of single trial-averaged traces
        dataMean

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign) 
        % numeric array of single trial-averaged trace standard error
        dataSem
        
        % size(data) scalar array indicating how many
        % trials contributed to data in each cell
        dataNTrials
        
        % nAlign x nBases x nCondition logical array indicating whether there is  
        % valid data in the corresponding cell, which is required when the conditions
        % matrix isn't complete (all conditions valid) or only some alignments are valid
        % for a given condition, etc.
        dataValid
        
        % nAlign x 1 cell with nBases x nCondition x nTime(iAlign) x nRandomSamples numeric
        % tensors containing randomly generated dataMean
        dataMeanRandomized

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the lower confidence interval value for
        % dataMean
        dataIntervalLow
        
        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the upper confidence interval value for
        % dataMean
        dataIntervalHigh
    end
    
    % Properties within *Manual properties store manually-specified values for each of the
    % above properties. These are used to store persistent copies of
    % the data when .dataSourceManual is true
    properties(Hidden, Access=protected) 
        basisNamesManual
        basisUnitsManual
        alignSummaryDataManual
        alignSummaryAggregatedManual
        basisAlignSummaryLookupManual
        dataByTrialManual
        tMinForDataByTrialManual
        tMaxForDataByTrialManual
        alignValidByTrialManual
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
        dataValidManual
        dataMeanRandomizedManual
    end
    
    % Dependent properties which we compute on the fly rather than cache
    properties(Dependent)
        % indicates how the elements in data are computed from individual trials
        % true means all trials correspond one-to-one with each other
        % false means trials were not simultaneous, implying that only
        % trial-averages should be compared
        simultaneous = false;
        
        % number of data sources used across all bases
        nDataSources
        
        nBases % number of bases (e.g. units, analog channels)

        nConditions % number of conditions in conditionDescriptor
        
        nAlign % number of alignments in alignDescriptorSet
        
        conditions % shortcut to conditionDescriptor.conditions
        
        conditionsSize % pass-thru to .conditionDescriptor
        
        alignNames % names pulled from the alignDescriptors

        conditionNames % condition names pulled from conditionDescriptor
        
        dataIntervalQuantilesAsString % dataIntervalQuantileLow/High summarized as a string
    end
   
    % Constructor, initialization, cache invalidation
    methods
        function pset = PopulationTrajectorySet()
            pset = pset.initialize();
        end
    end
    
    methods(Static)
        function pset = loadobj(pset)
            pset = builtin('loadobj', pset);
            pset = pset.initialize();
        end
    end
    
    methods
        function pset = saveobj(pset)

        end
        
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
            
            if isempty(pset.translationNormalization) || pset.translationNormalization.nBases ~= pset.nBases
                pset.translationNormalization = ...
                    StateSpaceTranslationNormalization.buildIdentityForPopulationTrajectorySet(pset);
            end
            
            if isempty(pset.timeDelta)
                pset.timeDelta = 1;
            end 
            
            if isempty(pset.spikeFilter)
                pset.spikeFilter = SpikeFilter.getDefaultFilter();
            end
            
            if isempty(pset.minTrialsForTrialAveraging)
                pset.minTrialsForTrialAveraging = 1;
            end
            
            if isempty(pset.minFractionTrialsForTrialAveraging)
                pset.minFractionTrialsForTrialAveraging = 0;
            end
            
            if isempty(pset.includeOnlyTrialsValidAllAlignments)
                pset.includeOnlyTrialsValidAllAlignments = false;
            end
            
            if isempty(pset.dataIntervalQuantileLow)
                pset.dataIntervalQuantileLow = 0.025;
            end
            
            if isempty(pset.dataIntervalQuantileHigh)
                pset.dataIntervalQuantileHigh = 0.975;
            end
            
            if isempty(pset.nRandomSamples)
                pset.nRandomSamples = 100;
            end
            
            if isempty(pset.randomSeed)
                pset.randomSeed = 0;
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
            
            pset = pset.invalidateCache();
        end
            
        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function pset = invalidateCache(pset)
            pset.warnIfNoArgOut(nargout);

            % copy before writing to odc!
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flush();
            end
        end
        
        function pset = invalidateTrialAveragedData(pset)
            pset.warnIfNoArgOut(nargout);
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushTrialAveragedData();
            end
        end

        function pset = invalidateRandomizedTrialAveragedData(pset)
            pset.warnIfNoArgOut(nargout);
            if ~isempty(pset.odc)
                pset.odc = pset.odc.copy();
                pset.odc.flushRandomizedTrialAveragedData();
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
            tcprintf('inline', '{yellow}%s: {bright white}%d bases, %d conditions, %d alignments, %s\n', ...
                class(pset), pset.nBases, pset.nConditions, pset.nAlign, dataSourceStr);
            tcprintf('inline', '{yellow}Dataset: {none}%s\n\n', pset.datasetName);
            
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
            builtin('disp', pset);
        end 
    end
    
    % these methods are setters for property values which change the
    % behavior of the PTS. They automatically invalidate downstream cached
    % values that depend on the value of the property being set.
    methods
        function pset = set.spikeFilter(pset, v)
            % spikeFilter invalidates everything
            pset.spikeFilter = v;
            pset = pset.invalidateCache();
        end
        
        function pset = set.timeDelta(pset, v)
            % timeDelta invalidates everything
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
        
        function pset = set.includeOnlyTrialsValidAllAlignments(pset, v)
            % only affects trial averaging
            pset.includeOnlyTrialsValidAllAlignments = v;
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
    end
    
    % General utility methods
    methods(Access=protected)
        function warnIfNoArgOut(obj, nargOut)
        % call using obj.warnIfNoArgOut(nargout);
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('WARNING: %s is not a handle class. If the instance returned by this method is not stored, this call has no effect.\\n', ...
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
            if ~pset.dataSourceManual
                v = pset.odc.alignSummaryAggregated;
                if isempty(v)
                    pset.buildAlignSummaryData();
                    v = pset.odc.alignSummaryAggregated;
                end
            else
                v = pset.alignSummaryAggregatedManual;
            end
        end 
        
        function pset = set.alignSummaryAggregated(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.alignSummaryAggregated = v;
            else
                pset.alignSummaryAggregatedManual = v;
            end
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
        
        function v = get.tvecDataMean(pset)
            v = pset.odc.tvecDataMean;
            if isempty(v)
                % compute on-the-fly and cache
                tvecDataMean = cellvec(pset.nAlign);
                for iAlign = 1:pset.nAlign
                    tvecDataMean{iAlign} = makecol(pset.tMinForDataMean(iAlign):pset.timeDelta:pset.tMaxForDataMean(iAlign));
                end
                pset.odc.tvecDataMean = tvecDataMean;
                v = pset.odc.tvecDataMean;
            end
        end 
        
        function v = get.nTimeDataMean(pset)
            v = pset.odc.nTimeDataMean;
            if isempty(v)
                % compute on-the-fly and cache
                pset.odc.nTimeDataMean = makecol(cellfun(@numel, pset.tvecDataMean));
                v = pset.odc.nTimeDataMean;
            end
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
                    pset.buildDataMean();
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
                    pset.buildDataMean();
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
       
        function v = get.dataMeanRandomized(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataMeanRandomized;
                if isempty(v)
                    pset.buildDataRandomized();
                    v = pset.odc.dataMeanRandomized;
                end
            else
                v = pset.dataMeanRandomizedManual;
            end
        end
        
        function pset = set.dataMeanRandomized(pset, v)
            if ~pset.dataSourceManual
                pset.odc = pset.odc.copy();
                pset.odc.dataMeanRandomized = v;
            else
                pset.dataMeanRandomizedManual = v;
            end
        end
        
        function v = get.dataIntervalLow(pset)
            if ~pset.dataSourceManual
                v = pset.odc.dataIntervalLow;
                if isempty(v)
                    pset.buildDataRandomized();
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
                if isempty(v)
                    pset.buildDataRandomized();
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
    end
  
    % methods which set and apply the conditionDescriptor,
    % alignDescriptorSet, and translationNormalization applied to this pset
    methods 
        function pset = setConditionDescriptor(pset, cd)
            pset.warnIfNoArgOut(nargout);
            assert(isequal(class(cd), 'ConditionDescriptor'), ...
                'Must be ConditionDescriptor instance');
            
            assert(cd.allAxisValueListsManual, 'ConditionDescriptor must have manual axis value lists. Use .setAxisValueList or .fixValueListsByApplyingToTrialData');
            
            pset.conditionDescriptor = cd.getConditionDescriptor();
            
            pset = pset.applyConditionDescriptor();
        end
        
        function pset = applyConditionDescriptor(pset)
            pset.warnIfNoArgOut(nargout);
            
            % apply the condition descriptor to each trial data now and
            % store the resulting tdca
            prog = ProgressBar(pset.nDataSources, 'Grouping trials in data sources');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.setConditionDescriptor(pset.conditionDescriptor);
            end
            prog.finish();
            
            % only trial averaging needs to be done again
            pset = pset.invalidateTrialAveragedData();
        end
        
        function pset = setAlignDescriptorSet(pset, adSet)
            pset.warnIfNoArgOut(nargout);
            
            if ~iscell(adSet)
                adSet = {adSet};
            end
            assert(all(cellfun(@(x) isequal(class(x), 'AlignDescriptor'), ...
                adSet)), 'Must be AlignDescriptor instance');
            
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
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.align(pset.alignDescriptorSet{:});
            end
            prog.finish();
            
            % changing alignments invalidates everything
            pset = pset.invalidateCache();
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
            
            if pset.dataSourceManual
                % for manual data, we apply this now to all data stored
                % manually since it will not be regenerated later
                if ~isempty(pset.dataByTrial)
                    pset.dataByTrial = pset.dataByTrial;
                end
                if ~isempty(pset.dataMean)
                    pset.dataMean = cellApplyToDataFn(pset.dataMean);
                end
                if ~isempty(pset.dataIntervalHigh)
                    pset.dataIntervalHigh = cellApplyToDataFn(pset.dataIntervalHigh);
                end
                if ~isempty(pset.dataIntervalLow)
                    pset.dataIntervalLow = cellApplyToDataFn(pset.dataIntervalLow);
                end
            else
                % for auto-computed data, we can check whether the data has
                % been already computed (by checking the odc properties).
                % if so we can do the transformation to save time.
                % otherwise, we can defer until the buildData* methods
                % apply the translation / normalization
                if ~isempty(pset.odc.dataByTrial)
                    pset.dataByTrial = pset.dataByTrial;
                end
                if ~isempty(pset.odc.dataMean)
                    pset.dataMean = cellApplyToDataFn(pset.dataMean);
                end
                if ~isempty(pset.odc.dataIntervalHigh)
                    pset.dataIntervalHigh = cellApplyToDataFn(pset.dataIntervalHigh);
                end
                if ~isempty(pset.odc.dataIntervalLow)
                    pset.dataIntervalLow = cellApplyToDataFn(pset.dataIntervalLow);
                end
                
                % not necessary to invalidate if we're careful about making
                % updates to all derived quantities here.
                
                % and invalidate the trial averaged data to reflect this
                % normalization
                %pset = pset.invalidateTrialAveragedData();
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
                if ~isempty(pset.dataIntervalHigh)
                    pset.dataIntervalHigh = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalHigh, 'UniformOutput', false);
                end
                if ~isempty(pset.dataIntervalLow)
                    pset.dataIntervalLow = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalLow, 'UniformOutput', false);
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
                if ~isempty(pset.odc.dataIntervalHigh)
                    pset.dataIntervalHigh = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalHigh, 'UniformOutput', false);
                end
                if ~isempty(pset.odc.dataIntervalLow)
                    pset.dataIntervalLow = cellfun(@trNorm.undoTranslationNormalizationToData, pset.dataIntervalLow, 'UniformOutput', false);
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
        
        function pset = meanSubtractBases(pset)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translate(-pset.computeMeanByBasis(), 'translationDescription', 'mean-subtracted');
        end
        
        function pset = normalizeBasesByStd(pset)
            pset.warnIfNoArgOut(nargout);
            pset = pset.normalize(pset.computeStdByBasis(), 'normalizationDescription', 'std-normalized');
        end
        
        function pset = zscoreByBasis(pset)
            pset.warnIfNoArgOut(nargout);
            pset = pset.translateNormalize(-pset.computeMeanByBasis(), pset.computeStdByBasis(), ...
                'translationDescription', 'mean-subtracted', ...
                'normalizationDescription', 'std-normalized');
        end
    end
    
    % build methods for the odc properties, each must store results 
    % directly into the odc (without copying first, which allows the
    % results to persist)
    methods 
        function src = getDataSourceForBasis(pset, iBasis)
            % return the data source which provided that data for basis
            % iBasis
            src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};         
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
                    basisNames = chName;
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
        
        function buildDataByTrial(pset)
            % stores dataByTrial, alignValidByTrial, tMinByTrial, and tMaxByTrial in odc
            % This method fetches the aligned data for EVERY trial in 
            % each data source, for each alignment. It does not consider
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
            
            prog = ProgressBar(pset.nBases, 'Extracting data by basis');
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                
                % request the specified aligned analog channel from the
                % specified data source.
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                chName = pset.basisDataSourceChannelNames{iBasis};
                
                % unapply the condition descriptor so that we can grab
                % all trials in this call, even the ones that would be
                % marked invalid by this condition info. Manually
                % invalid trials will still not be considered.
                src = src.resetConditionInfo();
                
                for iAlign = 1:pset.nAlign
                    % mark this align as active
                    src = src.useAlign(iAlign);

                    % currently will request either analog trials or
                    % filtered spike rates channels
                    if src.hasAnalogChannel(chName)
                        [mat, tvec] = src.getAnalogAsMatrix(chName, 'timeDelta', pset.timeDelta);
                    elseif src.hasSpikeChannelOrUnit(chName)
                        src = src.padForSpikeFilter(pset.spikeFilter);
                        [mat, tvec] = src.getSpikeRateFilteredAsMatrix(chName, ...
                            'spikeFilter', pset.spikeFilter, 'timeDelta', pset.timeDelta);
                    else
                        error('Unknown channel type');
                    end
                    
                    dataByTrial{iBasis, iAlign} = mat;
                    
                    % store the time limits used in the time vector for 
                    % the dataByTrial{..} matrix
                    tMinForDataByTrial(iBasis, iAlign) = min(tvec);
                    tMaxForDataByTrial(iBasis, iAlign) = max(tvec);
                    
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
                prog.finish();
            end
            
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
        
        function buildDataMean(pset)
            % computes and stores dataMean, dataIntervalHigh/Low, and dataNTrials into odc
            % this function computes summary statistics across trials,
            % especially dataMean. It selects time windows which are
            % consistent across all bases and conditions for a specific
            % align to simplify subsequent operations.
                       
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
            
            if pset.dataSourceManual
                return;
            end
            
            % first, compute the all-inclusive time window for each basis,
            % for each alignment, using only condition and align valid trials
            [tMinValidByBasisAlign, tMaxValidByBasisAlign] = deal(nan(pset.nBases, pset.nAlign));
            prog = ProgressBar(pset.nBases, 'Computing trial-averaged time windows by basis');
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                for iAlign = 1:pset.nAlign
                    % note, this src will not be aligned to this iAlign,
                    % but this isn't necessary since we've already
                    % extracted the aligned data
                    src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                    tMinByTrial = pset.tMinByTrial{iBasis, iAlign};
                    tMaxByTrial = pset.tMaxByTrial{iBasis, iAlign};
                    
                    % filter for trials which are valid for this alignment
                    % and for the condition grouping.
                    
                    if pset.includeOnlyTrialsValidAllAlignments
                        % look at validity across all alignments
                        alignValidByTrial = all(cell2mat(pset.alignValidByTrial{iBasis, :}));
                    else
                        % precomputed in buildDataByTrial to avoid re-aligning
                        % this encompasses manual invalids and alignInvalids
                        alignValidByTrial = pset.alignValidByTrial{iBasis, iAlign};
                    end
                    
                    % critical to use computedValid here because the src
                    % is aligned to alignDescriptorSet{1} and .valid will
                    % be affected by the alignInfo as well
                    conditionValidByTrial = src.conditionInfo.computedValid;
                    
                    % also ignore trials manually marked as invalid in the
                    % trial data itself, which won't be reflected in the
                    % .computedValid properties of AlignInfo or
                    % ConditionInfo
                    manualValidByTrial = src.getManualValid();
                    
                    validByTrial = alignValidByTrial & conditionValidByTrial & manualValidByTrial;
                  
                    % for each basis, align, take the largest window across
                    % all VALID trials for this condition, align
                    tMinValidByBasisAlign(iBasis, iAlign) = nanmin(tMinByTrial(validByTrial));
                    tMaxValidByBasisAlign(iBasis, iAlign) = nanmax(tMaxByTrial(validByTrial));
                end
            end
            prog.finish();
            
            % then, compute the largest window that is valid for ALL bases,
            % for each align. these will be nAlign x 1
            tMinForDataMean = makecol(nanmax(tMinValidByBasisAlign, [], 1));
            tMaxForDataMean = makecol(nanmin(tMaxValidByBasisAlign, [], 1));
            
            % number of time points for each alignment
            nTimeByAlign = arrayfun(@(mn, mx) numel(mn:pset.timeDelta:mx), ...
                tMinForDataMean, tMaxForDataMean);
            
            assert(all(tMinForDataMean <= tMaxForDataMean), 'No time window is valid across all bases');
            
            % now that we've determined the time window, we can compute the
            % trial average using data from these windows
            
            [dataMean, dataSem] = deal(cellvec(pset.nAlign));
            for iAlign = 1:pset.nAlign
                [dataMean{iAlign}, dataSem{iAlign}] = ...
                    deal(nan(pset.nBases, pset.nConditions, nTimeByAlign(iAlign)));
            end
            [dataValid, dataNTrials] = deal(nan(pset.nAlign, pset.nBases, pset.nConditions));
 
            prog = ProgressBar(pset.nBases, 'Computing trial-averaged data by basis');
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                for iAlign = 1:pset.nAlign
                    src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                    % pull the by-trial data from .dataByTrial, and use the
                    % conditioned tdca source to group the trials by
                    % condition
                    byTrial = pset.dataByTrial{iBasis, iAlign};
                    
                    % lookup the time limits which describe the byTrial
                    % matrix
                    tMinAll = pset.tMinForDataByTrial(iBasis, iAlign);
                    tMaxAll = pset.tMaxForDataByTrial(iBasis, iAlign);
                    tvecAll = tMinAll:pset.timeDelta:tMaxAll;
                    
                    % lookup the new time limits which we'll compute the
                    % trial average within
                    tMinValid = tMinForDataMean(iAlign);
                    tMaxValid = tMaxForDataMean(iAlign);
                    tMaskValid = tvecAll >= tMinValid & tvecAll <= tMaxValid;
                    
                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    byTrialValid = byTrial(:, tMaskValid);
                    byCondition = src.groupElements(byTrialValid);
                    
                    for iCondition = 1:numel(byCondition)
                        mat = byCondition{iCondition};
                        nTrials = size(mat, 1);
                        
                        % minimum trial count at each time point needed to
                        % compute an average, otherwise NaN
                        minTrials = max(pset.minTrialsForTrialAveraging, ...
                            ceil(nTrials * pset.minFractionTrialsForTrialAveraging));
                        
                        nTrialsByTime = sum(~isnan(mat), 1);
                        
                        % compute mean, sem, and trial count
                        m = nanmean(mat, 1)';
                        m(nTrialsByTime < minTrials) = NaN;
                        
                        se = nansem(mat, 1)';
                        se(nTrialsByTime < minTrials) = NaN;

                        dataMean{iAlign}(iBasis, iCondition, :) = m;
                        dataSem{iAlign}(iBasis, iCondition, :) = se;
                        dataNTrials(iAlign, iBasis, iCondition) = nTrials;
                        dataValid(iAlign, iBasis, iCondition) = size(mat, 1) > 0;
                    end
                end
            end
            prog.finish();
            
            % no need to apply translation / normalization here since it is
            % already applied to dataByTrial!
            
            % store in odc without copying
            c = pset.odc;
            c.dataMean = dataMean;
            c.dataSem = dataSem;
            c.dataNTrials = dataNTrials;
            c.dataValid = dataValid;
            c.nTimeDataMean = nTimeByAlign;
            c.tMinForDataMean = tMinForDataMean;
            c.tMaxForDataMean = tMaxForDataMean;
        end
        
        function buildDataRandomized(pset)
            % build dataMeanRandomized and store in odc without copying
            prog = ProgressBar(pset.nBases, 'Computing randomized trial-averaged data by basis');
            seed = pset.randomSeed;
            
            dataMeanRandomized = cellvec(pset.nAlign); 
            for iAlign = 1:pset.nAlign
                dataMeanRandomized{iAlign} = nan(pset.nBases, pset.nConditions, ...
                    pset.nTimeDataMean(iAlign), pset.nRandomSamples);
            end
            
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                
                 listByConditionSamples = src.conditionInfo.generateMultipleRandomizedListByCondition(...
                        pset.nRandomSamples, 'initialSeed', seed, 'showProgress', false);
                 seed = seed + src.nTrials;

                for iAlign = 1:pset.nAlign
                    % pull the by-trial data from .dataByTrial
                    byTrial = pset.dataByTrial{iBasis, iAlign};
                    
                    % lookup the time limits which describe the byTrial
                    % matrix
                    tMinAll = pset.tMinForDataByTrial(iBasis, iAlign);
                    tMaxAll = pset.tMaxForDataByTrial(iBasis, iAlign);
                    tvecAll = tMinAll:pset.timeDelta:tMaxAll;
                    
                    % lookup the new time limits which we'll compute the
                    % trial average within
                    tMinValid = pset.tMinForDataMean(iAlign);
                    tMaxValid = pset.tMaxForDataMean(iAlign);
                    tMaskValid = tvecAll >= tMinValid & tvecAll <= tMaxValid;
              
                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    byTrialValid = byTrial(:, tMaskValid);
                    
                    byCondition = src.groupElements(byTrialValid);
                    minTrialsByCondition = nanvec(numel(byCondition));
                    for iCondition = 1:numel(byCondition)
                        mat = byCondition{iCondition};
                        nTrials = size(mat, 1);

                        % minimum trial count at each time point needed to
                        % compute an average, otherwise NaN
                        minTrialsByCondition(iCondition) = max(pset.minTrialsForTrialAveraging, ...
                            ceil(nTrials * pset.minFractionTrialsForTrialAveraging));
                    end
                    
                    for iSample = 1:pset.nRandomSamples
                        byCondition = cellfun(@(idx) byTrialValid(idx,:), ...
                            listByConditionSamples{iSample}, 'UniformOutput', false);
                    
                        for iCondition = 1:numel(byCondition)
                            mat = byCondition{iCondition};
                            nTrialsByTime = sum(~isnan(mat), 1);
                            
                            % compute mean
                            m = nanmean(mat, 1)';
                            m(nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;
                            
                            dataMeanRandomized{iAlign}(iBasis, iCondition, :, iSample) = m;
                        end
                    end
                end
            end
            prog.finish();
                
            % compute the quantiles to use as intervals
            [dataIntervalHigh, dataIntervalLow] = deal(cellvec(pset.nAlign));
            qLow = pset.dataIntervalQuantileLow;
            qHigh = pset.dataIntervalQuantileHigh;
            for iAlign = 1:pset.nAlign
                dataIntervals = quantile(dataMeanRandomized{iAlign}, [qLow qHigh], 4);
                dataIntervalLow{iAlign} = dataIntervals(:, :, :, 1);
                dataIntervalHigh{iAlign} = dataIntervals(:, :, :, 2);
            end

            c = pset.odc;
            c.dataIntervalHigh = dataIntervalHigh;
            c.dataIntervalLow = dataIntervalLow;
            c.dataMeanRandomized = dataMeanRandomized;
        end
        
        function buildAlignSummaryData(pset)
            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
            basisAlignSummaryLookup = pset.basisDataSourceIdx;
            alignSummaryData = cell(pset.nDataSources, pset.nAlign);
            
            % copy the align summary data from each data source, which may
            % take time since it is typically computed on the fly
            prog = ProgressBar(pset.nDataSources, 'Computing alignment summary statistics by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                alignSummaryData(iSrc, :) = pset.dataSources{iSrc}.alignSummarySet;
            end
            prog.finish();
            
            % and build the aggregated data across all bases too, for each
            % alignment
            alignSummaryAggregated = cell(pset.nAlign, 1);
            prog = ProgressBar(pset.nAlign, 'Computing aggregate alignment summary statistics');
            for iAlign = 1:pset.nAlign
                prog.update(iAlign);
                alignSummaryAggregated{iAlign} = AlignSummary.buildByAggregation(alignSummaryData(:, iAlign));
            end 
            prog.finish();
            
            c = pset.odc;
            c.basisAlignSummaryLookup = basisAlignSummaryLookup;
            c.alignSummaryData = alignSummaryData;
            c.alignSummaryAggregated = alignSummaryAggregated;
        end
    end
    
    methods % Filtering bases, conditions (NOT WORKING)
%         function filterAlign(pset, idx)
%             p = inputParser;
%             p.addRequired('alignIdx', @isvector);
%             p.parse(idx);
%             alignIdx = makecol(p.Results.alignIdx);
% 
%             pset.data = pset.data(:, :, alignIdx, :);
%             pset.dataSem = pset.dataSem(:, :, alignIdx, :);
%             pset.timeData = pset.timeData(:, :, alignIdx, :);
%             pset.alignTimeInfoData = pset.alignTimeInfoData(:, :, alignIdx, :);
%             pset.alignDescriptorSet = pset.alignDescriptorSet(alignIdx);
%             pset.dataValid = pset.dataValid(:, :, alignIdx);
%             pset.dataNTrials = pset.dataNTrials(:, :, alignIdx);
%             pset.tMinDataManual = pset.tMinDataManual(:, :, alignIdx);
%             pset.tMinDataManual = pset.tMaxDataManual(:, :, alignIdx);
%         end
%         
%         function filterBases(pset, idx)
%             % keep only bases listed in or masked by idx
%             p = inputParser;
%             p.addRequired('basisIdx', @isvector);
%             p.parse(idx);
%             basisIdx = makecol(p.Results.basisIdx);
% 
%             pset.data = pset.data(basisIdx, :, :);
%             pset.timeData = pset.timeData(basisIdx, :, :);
%             pset.alignTimeInfoData = pset.alignTimeInfoData(basisIdx, :, :);
%             pset.dataValid = pset.dataValid(basisIdx, :, :);
%             pset.dataNTrials = pset.dataNTrials(basisIdx, :, :);
%             
%             pset.tMinDataManual = pset.tMinDataManual(basisIdx, :, :);
%             pset.tMaxDataManual = pset.tMaxDataManual(basisIdx, :, :);
% 
%             pset.basisNames = makecol(pset.basisNames(basisIdx));
%             pset.basisMeta = makecol(pset.basisMeta(basisIdx));
% 
%             if pset.storeDataSources
%                 pset.dataSources = pset.dataSources(basisIdx, :);
%                 pset.dataSourcesOrig = pset.dataSourcesOrig(basisIdx, :);
%             end
%         end

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

    methods % Compute on-the-fly Dependent properties
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

    methods % Simple statistics     
        function meanByBasis = computeMeanByBasis(pset)
            CTAbyN = pset.buildCTAbyN();
            meanByBasis = nanmean(CTAbyN, 1)';
        end
        
        function varByBasis = computeVarByBasis(pset)
            varByBasis = nanvar(pset.buildCTAbyN(), 0, 1)';
        end
        
        function stdByBasis = computeStdByBasis(pset)
            stdByBasis = nanstd(pset.buildCTAbyN(), 0, 1)';
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
%             p.addParamValue('basisIdx', [1:6], @(x) isvector(x) && ...
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

        function plotBases(pset, varargin)
            % plot bases one above the next 
            p = inputParser;
            p.addParamValue('basisIdx', 1:pset.nBases, @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            p.addParamValue('conditionIdx', 1:pset.nConditions, @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nConditions])));
            % plot each basis at it's original scale or normalized to fit
            % the bands
            p.addParamValue('normalize', true, @islogical);
            % each basis starts at y = 1, 2, 3 and data is scaled to fit a
            % band with y-height scaling
            p.addParamValue('scaling', 0.8, @isscalar);
            
            % usually plot first basis at top, last basis at bottom
            p.addParamValue('reverse', false, @islogical);
            
            p.addParamValue('alignGapFraction', 0.02, @isscalar);
            p.addParamValue('xOffset', 0, @isscalar);
            p.addParamValue('yOffset', 0, @isscalar);
            p.addParamValue('plotArgs', {}, @iscell)
            p.parse(varargin{:});

            basisIdx = p.Results.basisIdx;
            nBasesPlot = numel(basisIdx);
            conditionIdx = p.Results.conditionIdx;
            nConditionsPlot = numel(conditionIdx);
            
            scaling = p.Results.scaling;
            normalize = p.Results.normalize;
            reverse = p.Results.reverse;
            alignGapFraction = p.Results.alignGapFraction;
            xOffset = p.Results.xOffset;
            yOffset = p.Results.yOffset;
            plotArgs = p.Results.plotArgs;
            
            timeWidthByAlign = pset.nTimeDataMean*pset.timeDelta;
            
            % compute absolute x-gap between alignments
            alignGap = alignGapFraction*sum(timeWidthByAlign) / (1 - alignGapFraction*pset.nAlign);
            
            for iAlign = 1:pset.nAlign
                tvec = pset.tvecDataMean{iAlign};
                tvec = tvec - min(tvec) + xOffset;
                
                data = pset.dataMean{iAlign}(basisIdx, conditionIdx, :);
                
                if normalize
                    % all data will range from 0 to 1 by row
                    mat = reshape(data, nBasesPlot, nConditionsPlot * numel(tvec));
                    mat = bsxfun(@minus, mat, nanmin(mat, [], 2));
                    mat = bsxfun(@rdivide, mat, nanmax(mat, [], 2));
                    data = reshape(mat, nBasesPlot, nConditionsPlot, numel(tvec));
                else
                    % all data will range from 0 to 1, but keep the same
                    % relative scaling
                    data = data - nanmin(data(:));
                    data = data / nanmax(data(:));
                end
                
                data = data * scaling + yOffset;
                data = bsxfun(@plus, data, (nBasesPlot:-1:1)');
                
                if reverse
                    data = flipud(data);
                end
                
                for iCondition = 1:nConditionsPlot
                    dataC = squeeze(data(:, iCondition, :));
                    plotArgsC = pset.conditionDescriptor.appearances(conditionIdx(iCondition)).getPlotArgs();
                    plot(tvec, dataC, '-', plotArgsC{:}, plotArgs{:})
                    hold on
                end
                
                xOffset = xOffset + timeWidthByAlign(iAlign) + alignGap;
            end
            
            box off;
        end

        function plotStateSpace(pset, varargin)
            % plot a 2d or 3d basis1 x basis2 x basis3 trajectory plot
            p = inputParser;
            p.addParamValue('basisIdx', 1:min(pset.nBases, 3), @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            p.addParamValue('plotArgs', {}, @iscell)
            p.parse(varargin{:});
            
            basisIdx = p.Results.basisIdx;
            nBasesPlot = numel(basisIdx);
            plotArgs = p.Results.plotArgs;

            switch(nBasesPlot)
                case 2
                    use3d = false;
                case 3
                    use3d = true;
                otherwise
                    error('Number of bases must be 2 or 3');
            end

            for iAlign = 1:pset.nAlign
                %tvec = pset.tvecDataMean{iAlign};
                data = pset.dataMean{iAlign};

                for iCondition = 1:pset.nConditions
                    appear = pset.conditionDescriptor.appearances(iCondition);
                    plotArgsC = appear.getPlotArgs();
                    
                    dataMat = squeeze(data(basisIdx, iCondition, :));
                    
                    if use3d
                        plot3(dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                            plotArgsC{:}, plotArgs{:});
                    else
                        dataMat = [dataVec1 dataVec2];
                        plot(dataMat(1, :), dataMat(2, :), ...
                            plotArgsC{:}, plotArgs{:});
                    end

                    hold on
                end
            
                % annotate data with marks / intervals
                for iBasis = 1:pset.nBases
                    %as = pset.alignSummaryData{pset.basisAlignSummaryLookup(iBasis), iAlign};
                    
                    % data is nBases x C x T
                    % drawOnTimeseriesByConditions needs T x nBasesPlot x C
                    %dataForDraw = permute(data(basisIdx, :, :), [3 1 2]);
                    %as.drawOnTimeseriesByCondition(dataForDraw); 
                end 
            end

            box off
            xlabel(pset.basisNames{basisIdx(1)});
            ylabel(pset.basisNames{basisIdx(2)});

             if use3d
                zlabel(pset.basisNames{basisIdx(3)});
                view([-40 20]);
            end

            axis tight
            axis square
            axis vis3d
        end
        
    end

    methods % Build data matrices
        % Notation:
        % C is nConditions
        % T is nTimepoints for a given alignment, such that T*A really
        %     means T(align1) + T(align2) + T(align3) + ...
        % A is nAlign
        % N is nBases
        
        function [NbyCbyTA, tvec, avec] = buildNbyCbyTA(pset, varargin)
            [NbyCbyTA, avec] = TensorUtils.catWhich(3, pset.dataMean{:});
            tvec = cat(1, pset.tvecDataMean{:});
        end

        function [CTAbyN, cvec, tvec, avec, nvec] = buildCTAbyN(pset, varargin)
            % out is C*T*A x N concatenated matrix

            [NbyCbyTA, tvec, avec] = pset.buildNbyCbyTA();
            nvec = (1:pset.nBases)';
            cvec = (1:pset.nConditions)';
            labels = {nvec, cvec, [tvec, avec]};
            [CTAbyN, labelsOut] = TensorUtils.reshapeByConcatenatingDims(NbyCbyTA, {[3 2], 1}, labels);

            cvec = labelsOut{1}(:, 1);
            tvec = labelsOut{1}(:, 2);
            avec = labelsOut{1}(:, 3);
            nvec = labelsOut{2};
        end

         function NbyTAbyAttr = buildNbyTAbyConditionAttributes(pset, varargin)
             % build a tensor of N by T by nValsAttr1 by nValsAttr2 x ...
             % this tensor is the format used by dpca_covs
             NbyTAbyC = permute(pset.buildNbyCbyTA(), [1 3 2]);
             N = pset.nBases;
             TA = size(NbyTAbyC, 2);
             condSize = pset.conditionDescriptor.conditionsSize;
             NbyTAbyAttr = reshape(NbyTAbyC, [N TA makerow(condSize)]);
         end

%         function [out tvecByConditionAlign] = buildCTByNEachA(pset, varargin)
%             % out: A cell vector of CT x N matrices 
%             % timeVec: nAlign cell of time vectors common to alignment
%             
%             p = inputParser();
%             p.addParamValue('timeValidAcrossConditions', false, @islogical)
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
%             p.addParamValue('timeValidAcrossConditions', false, @islogical)
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

    methods % Comparative statistics
        function [distByFromAlign, timeVecByFromAlign] = getDistanceBetween(pset, cFromList, cToList, varargin)
            % dist: length(fromAlign) cell of T x length(cFromList) distance traces as columns
            % comparing each condition cFromList(i) to condition cToList(i)

            p = inputParser;
            % if true, search the entire from trajectory for the closest point to 
            % each point on the to trajectory looking across multiple alignments as well.
            % If false, compute distance at each timepoint separately
            p.addParamValue('searchEntire', true, @islogical); 
            % calculate distances within this subset of bases
            p.addParamValue('basisIdx', true(pset.nBases, 1), @isvector);

            % when searchEntire is true, calculate distances from trajectories (with condition cTo) along
            % each of the alignments indexed in fromAlign, to the closest point in
            % trajectories in ANY of the alignments indexed in  toAlign
            % to the closest point on condition cTo within these alignments
            p.addParamValue('fromAlign', 1:pset.nAlign, @isvector);
            p.addParamValue('toAlign', 1:pset.nAlign, @isvector);
            
            % leave empty to compute a distance vs time trajectory for the
            % entire from trajectory at each alignment. Populate with a
            % nAlign x 1 cell array of time points to compute the distance
            % only from a specific set of time points
            p.addParamValue('timepointsByFromAlign', [], @(x) isempty(x) || iscell(x));

            p.addParamValue('showPlot', true, @islogical);

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

                    title(sprintf('Align %s', pset.alignNames{fromAlign(iIdxFromAlign)}));
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
            p.addParamValue('reverse', false, @islogical);
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
