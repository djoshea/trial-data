classdef(ConstructOnLoad) PopulationTrajectorySet
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

    properties(SetAccess=protected, Hidden, Transient)
        % odc is an instance of PopulationTrajectorySetOnDemandCache
        % which is a handle class. Properties which derive from 
        odc
    end
    
    % properties which control the behavior of the pset
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
        minTrialsForTrialAveraging = 1;
        
        % The minimum fraction of trials in a given condition over which to
        % compute a trial average, relative to the the total number of trials
        % in that condition. This parameter determines the valid time
        % windows for trial-averaged data (e.g. dataMean)
        minFractionTrialsForTrialAveraging = 0;
        
        % When multiple alignDescriptors are used to align the data, some
        % trials may be valid only for some alignments and not others.
        % Setting this flag to true ensures that only trials which are
        % valid for ALL alignments will be considered 
        includeOnlyTrialsValidAllAlignments = false;
        
        % number of randomized samples to draw when generating
        % dataMeanRandomized
        nRandomSamples = 100;
        
        % random seed to use as initial seed when generating random data
        % sets
        randomSeed = 0;
    end
    
    % alignDescriptors, conditionDescriptor, and translationNormalization
    % which align, group, and translate/normalize the data generated from
    % the dataSources. These may be set using eponymous methods prefixed with set*
    properties(SetAccess=protected)
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
        dataSourceManual = false;
        
        % TrialData data sources which source all data for the trajectories
        % this may be a single trial data object or many. If there is only
        % one, all bases are considered simultaneous.
        dataSources = {};
        
        % nBases x 1 index into dataSourceSet.
        basisDataSourceIdx
        
        % nBases x 1 cellstr indicating which channel name to extract data
        % from
        basisDataSourceChannelNames
    end
    
    % Properties whose values are computed dynamically and persist within odc
    % or are specified manually and persist within manualData
    properties(Hidden, Dependent, Transient, SetAccess=protected)
        % Alignment summary statistics by basis by align
        
        % nAlignSummary x nAlign cell containing AlignSummary instances for each
        % basis. built by buildAlignSummary
        alignSummaryData
        
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
        % alignmnet for the dataMean and dataSem cells
        tMinForDataMean
        tMaxForDataMean
        
        % trial-averaged data within each condition
        % nBases x nConditions x nAlign cell array of single trial-averaged traces
        dataMean
   
        % nBases x nConditions x nAlign cell array to firing rate traces
        dataSem
                        
        % size(data) scalar array indicating how many
        % trials contributed to data in each cell
        dataNTrials
        
        % logical array with size(data) indicating whether there is  
        % valid data in the corresponding cell, which is required when the conditions
        % matrix isn't complete (all conditions valid) or only some alignments are valid
        % for a given condition, etc.
        dataValid
        
        % nRandomSamples x 1 cell vector containing cell arrays like
        % dataMean which were generated by randomizing the 
        dataMeanRandomized
    end
    
    % Properties within *Manual properties store manually-specified values for each of the
    % above properties. These are used to store persistent copies of
    % the data when .dataSourceManual is true
    properties(Hidden, Access=protected) 
        basisNamesManual
        basisUnitsManual
        alignSummaryDataManual
        basisAlignSummaryLookupManual
        dataByTrialManual
        tMinForDataByTrialManual
        tMaxForDataByTrialManual
        alignValidByTrialManual
        tMinByTrialManual
        tMaxByTrialManual
        tMinForDataMeanManual
        tMaxForDataMeanManual
        dataMeanManual
        dataSemManual
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
        
        conditionsSize % pass-thru to .conditionDescriptor
        
        alignNames % names pulled from the alignDescriptors

        conditionNames % condition names pulled from conditionDescriptor
    end
   
    % Constructor, initialization, cache invalidation
    methods 
        function pset = PopulationTrajectorySet()
            pset = pset.initialize();
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
                pset = pset.applyFirstAlignDescriptor();
            end
            
            if isempty(pset.conditionDescriptor)
                pset.conditionDescriptor = ConditionDescriptor();
                pset = pset.applyConditionDescriptor();
            end
            
            if isempty(pset.spikeFilter)
                pset.spikeFilter = SpikeFilter.getDefaultFilter();
            end
            
            if isempty(pset.timeDelta)
                pset.timeDelta = 1;
            end 
            
            pset = pset.invalidateCache();
        end
            
        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function pts = invalidateCache(pts)
            pts.warnIfNoArgOut(nargout);

            % copy before writing to odc!
            pts.odc = pts.odc.copy();
            pts.odc.flush();
        end
        
        function pset = invalidateTrialAveragedData(pset)
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();
            pset.odc.flushTrialAveragedData();
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
            v = pset.odc.tMaxByTrial;
            if isempty(v)
                pset.buildDataByTrial();
                v = pset.odc.tMaxByTrialManual;
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
        
        function v = get.tMinForDataMean(pset)
            if ~pset.dataSourceManual
                v = pset.odc.tMinForDataMean;
                if isempty(v)
                    pset.buildDataMean();
                    v = pset.odc.tMinByTrial;
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
    end
  
    % methods which set and apply the conditionDescriptor,
    % alignDescriptorSet, and translationNormalization applied to this pset
    methods 
        function pset = setConditionDescriptor(pset, cd)
            pset.warnIfNoArgOut(nargout);
            assert(isequal(class(cd), 'ConditionDescriptor'), ...
                'Must be ConditionDescriptor instance');
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
            
            pset.alignDescriptorSet = adSet;
            
            pset = pset.applyFirstAlignDescriptor();
        end
        
        function pset = applyFirstAlignDescriptor(pset)
            pset.warnIfNoArgOut(nargout);
            
            % align all data sources to the FIRST alignDescriptor
            % since we have to hold onto one anyway
            prog = ProgressBar(pset.nDataSources, 'Aligning data sources to first alignDescriptor');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                pset.dataSources{iSrc} = pset.dataSources{iSrc}.align(pset.alignDescriptorSet{1});
            end
            prog.finish();
            
            % changing alignments invalidates everything
            pset = pset.invalidateCache();
        end
        
        function pset = applyTranslationNormalization(pset, trNorm)
            pset.warnIfNoArgOut(nargout);
            
            % only allow a single translation / normalization for
            % simplicity.
            if ~isempty(pset.translationNormalization)
                error('Translation / normalization already applied');
            end
            
            pset.translationNormalization = trNorm;
            
            if pset.dataSourceManual
                % for manual data, we apply this now to all data stored
                % manually since it will not be regenerated later
                if ~isempty(pset.dataByTrial)
                    pset.dataByTrial = trNorm.applyTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.dataMean)
                    pset.dataMean = trNorm.applyTranslationNormalizationToData(pset.dataMean);
                end
                if ~isempty(pset.dataSem)
                    pset.dataSem = trNorm.applyTranslationNormalizationToData(pset.dataSem);
                end
            else
                % for auto-computed data, we can check whether the data has
                % been already computed (by checking the odc properties).
                % if so we can do the transformation to save time.
                % otherwise, we can defer until the buildData* methods
                % apply the translation / normalization
                if ~isempty(pset.odc.dataByTrial)
                    pset.dataByTrial = trNorm.applyTranslationNormalizationToData(pset.dataByTrial);
                end
                
                % and invalidate the trial averaged data to reflect this
                % normalization
                pset = pset.invalidateTrialAveragedData();
            end
        end
        
        function pset = clearTranslationNormalization(pset, trNorm)
             pset.warnIfNoArgOut(nargout);
            
            if isempty(pset.translationNormalization)
                error('No translation / normalization applied');
            end
            
            if pset.dataSourceManual
                % for manual data, we reverse this now to all data stored
                % manually since it will not be regenerated later
                if ~isempty(pset.dataByTrial)
                    pset.dataByTrial = trNorm.undoTranslationNormalizationToData(pset.dataByTrial);
                end
                if ~isempty(pset.dataMean)
                    pset.dataMean = trNorm.undoTranslationNormalizationToData(pset.dataMean);
                end
                if ~isempty(pset.dataSem)
                    pset.dataSem = trNorm.undoTranslationNormalizationToData(pset.dataSem);
                end
            else
                % for auto-computed data, we can check whether the data has
                % been already computed (by checking the odc properties).
                % if so we can do the transformation to save time.
                % otherwise, we can defer until the buildData* methods
                % apply the translation / normalization
                if ~isempty(pset.odc.dataByTrial)
                    pset.dataByTrial = trNorm.undoTranslationNormalizationToData(pset.dataByTrial);
                end
                pset = pset.invalidateTrialAveragedData();
            end
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
            
            [dataByTrial, tMinByTrial, tMaxByTrial, alignValidByTrial] = ...
                deal(cell(pset.nBases, pset.nAlign));
            
            [tMinForDataByTrial, tMaxForDataByTrial] = deal(nan(pset.nBases, pset.nAlign));
           
            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
%             basisAlignSummaryLookup = pset.basisDataSourceIdx;
%             alignSummaryData = cell(pset.nDataSources, pset.nAlign);
            
            for iAlign = 1:pset.nAlign
                if iAlign == 1
                    % all data sources already aligned to the first
                    % alignDescriptor
                    aligned = pset.dataSources;
                else
                    % apply this alignment to each data source but hold on
                    % to the result separately (don't store it).
                    % also extract the alignSummary from each source now
                    aligned = cellvec(pset.nDataSources);
                    
                    ad = pset.alignDescriptorSet{iAlign};
                    for iSrc = 1:pset.nDataSources
                        aligned{iSrc} = pset.dataSources{iSrc}.align(ad);
                    end
                end
                
                prog = ProgressBar(pset.nBases, 'Extracting data by basis for align %d', iAlign);
                for iBasis = 1:pset.nBases
                    prog.update(iBasis);
                    
                    % request the specified aligned analog channel from the
                    % specified data source.
                    src = aligned{pset.basisDataSourceIdx(iBasis)};
                    chName = pset.basisDataSourceChannelNames{iBasis};
                    
                    % unapply the condition descriptor so that we can grab
                    % all trials in this call, even the ones that would be
                    % marked invalid by this condition info. Manually
                    % invalid trials will still not be considered.
                    src = src.ungroup();
                    
                    % currently will request either analog trials or
                    % filtered spike rates channels
                    if src.hasAnalogChannel(chName)
                        [mat, tvec] = src.getAnalogMatrix(chName, 'timeDelta', pset.timeDelta);
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
                    alignValidByTrial{iBasis, iAlign} = src.alignInfo.computedValid;
                    
                    % also store the precise time starts and stops for EACH
                    % trial that comprises that matrix. Essential that all
                    % padding be done to src before this call to ensure
                    % that the tvec returned above matches these numbers
                    [tMinByTrial{iBasis, iAlign}, ...
                     tMaxByTrial{iBasis, iAlign}] = src.getTimeStartStopEachTrial();
                 
                    assert(min(tvec) == nanmin(tMinByTrial{iBasis, iAlign}) && ...
                           max(tvec) == nanmax(tMaxByTrial{iBasis, iAlign}), ...
                        'Time vector returned by TrialData has invalid limits');
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
            % computes and stores dataMean, dataSem, and dataNTrials into odc
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
                    
                    % also ignore trials manually marked as invalid
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
            
            assert(all(tMinForDataMean <= tMaxForDataMean), 'No time window is valid across all bases');
            
            % now that we've determined the time window, we can compute the
            % trial average using data from these windows
            [dataMean, dataSem] = deal(cell(pset.nBases, pset.nConditions, pset.nAlign));
            [dataValid, dataNTrials] = deal(nan(pset.nBases, pset.nConditions, pset.nAlign));
 
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
                        
                        dataMean{iBasis, iCondition, iAlign} = m;
                        dataSem{iBasis, iCondition, iAlign} = se;
                        dataNTrials(iBasis, iCondition, iAlign) = nTrials;
                        dataValid(iBasis, iCondition, iAlign) = size(mat, 1) > 0;
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
            c.tMinForDataMean = tMinForDataMean;
            c.tMaxForDataMean = tMaxForDataMean;
        end
        
        function buildDataRandomized(pset)
            % build dataMeanRandomized and store in odc without copying
            prog = ProgressBar(pset.nBases, 'Computing randomized trial-averaged data by basis');
            seed = pset.randomSeed;
            
            dataMeanRandomized = cell(pset.nBases, pset.nConditions, pset.nAlign, pset.nRandomSamples);
            
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
                            
                            dataMeanRandomized{iBasis, iCondition, iAlign, iSample} = m;
                        end
                    end
                end
            end
            prog.finish();
            
            c = pset.odc;
            c.dataMeanRandomized = dataMeanRandomized;
        end
        
        function buildAlignSummaryData(pset)
            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
            basisAlignSummaryLookup = pset.basisDataSourceIdx;
            alignSummaryData = cell(pset.nDataSources, pset.nAlign);
            
            for iAlign = 1:pset.nAlign
                ad = pset.alignDescriptorSet{iAlign};
                prog = ProgressBar(pset.nBases, 'Computing alignment summary statistics for align %d', iAlign);
                for iSrc = pset.nDataSources
                    prog.update(iSrc);
                    if iAlign == 1
                        % already aligned
                        alignSummaryData{iSrc, iAlign} = pset.dataSources{iSrc}.alignSummary;
                    else
                        % need to realign the pset temporarily
                        alignSummaryData{iSrc, iAlign} = pset.dataSources{iSrc}.align(ad).alignSummary;
                    end
                end
                prog.finish();
            end
            
            c = pset.odc;
            c.basisAlignSummaryLookup = basisAlignSummaryLookup;
            c.alignSummaryData = alignSummaryData;
        end
    end
    
    methods % Filtering bases, conditions (NOT WORKING)
        function filterAlign(pset, idx)
            p = inputParser;
            p.addRequired('alignIdx', @isvector);
            p.parse(idx);
            alignIdx = makecol(p.Results.alignIdx);

            pset.data = pset.data(:, :, alignIdx, :);
            pset.dataSem = pset.dataSem(:, :, alignIdx, :);
            pset.timeData = pset.timeData(:, :, alignIdx, :);
            pset.alignTimeInfoData = pset.alignTimeInfoData(:, :, alignIdx, :);
            pset.alignDescriptorSet = pset.alignDescriptorSet(alignIdx);
            pset.dataValid = pset.dataValid(:, :, alignIdx);
            pset.dataNTrials = pset.dataNTrials(:, :, alignIdx);
            pset.tMinDataManual = pset.tMinDataManual(:, :, alignIdx);
            pset.tMinDataManual = pset.tMaxDataManual(:, :, alignIdx);
        end
        
        function filterBases(pset, idx)
            % keep only bases listed in or masked by idx
            p = inputParser;
            p.addRequired('basisIdx', @isvector);
            p.parse(idx);
            basisIdx = makecol(p.Results.basisIdx);

            pset.data = pset.data(basisIdx, :, :);
            pset.timeData = pset.timeData(basisIdx, :, :);
            pset.alignTimeInfoData = pset.alignTimeInfoData(basisIdx, :, :);
            pset.dataValid = pset.dataValid(basisIdx, :, :);
            pset.dataNTrials = pset.dataNTrials(basisIdx, :, :);
            
            pset.tMinDataManual = pset.tMinDataManual(basisIdx, :, :);
            pset.tMaxDataManual = pset.tMaxDataManual(basisIdx, :, :);

            pset.basisNames = makecol(pset.basisNames(basisIdx));
            pset.basisMeta = makecol(pset.basisMeta(basisIdx));

            if pset.storeDataSources
                pset.dataSources = pset.dataSources(basisIdx, :);
                pset.dataSourcesOrig = pset.dataSourcesOrig(basisIdx, :);
            end
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

    methods % Dependent properties
        function n = get.nDataSources(pset)
            n = numel(pset.dataSources);
        end
        
        function n = get.nBases(pset)
            n = numel(pset.basisDataSourceIdx);
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
    end
    
    methods % Resampling, shuffling : TODO for now this only works for SpikeRaster sources
        function shuffleSourcesAlong(pset, compareAlong)
            assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
            pset.dataSources = cellfun(@(sr) sr.buildShuffledAlong(compareAlong), pset.dataSourcesOrig, 'UniformOutput', false); 
            pset.updateFromDataSources();
        end

        function resampleSources(pset)
            assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
            pset.dataSources = cellfun(@(sr) sr.buildResampled(), pset.dataSourcesOrig, 'UniformOutput', false); 
            pset.updateFromDataSources();
        end

        function resampleSourcesFromSingleAttributeValue(pset, attr, value)
            assert(pset.storeDataSources, 'PopulationTrajectorySet must be configured with .storeDataSources==true in order to accomplish this');
            pset.dataSources = cellfun(@(sr) sr.buildResampledFromSingleAttributeValue(attr, value), pset.dataSourcesOrig, 'UniformOutput', false); 
            pset.updateFromDataSources();
        end
    end

    methods % Utility / internal mapping
        function varargout = alignBasisConditionDataFun(pset, fn, varargin)
            % runs fn(time, data[, iBasis, iCondition, iAlign)) 
            % on each alignment on each basis by condition, returns
            % results in nAlign cell array of {nBases x nCondition cell array}
            %
            % if multiple output arguments are requested, fn will be called
            % so as to request the same number of outputs
            %
            % asMatrix is a nargout x 1 boolean vector indicating whether to call
            % cell2mat on that output before returning so that you end up with
            % a nAlign cell array of nBases x nCondition matrices instead of cell array

            nOut = max(nargout, 1);

            p = inputParser;
            p.addParamValue('asMatrix', false(nOut, 1), @islogical);
            p.addParamValue('idxAlign', 1:pset.nAlign, @isvector);
            p.addParamValue('checkInvalid', true, @islogical);
            p.addParamValue('returnForInvalid', NaN, @(x) true); % can be scalar value or function handle
            p.parse(varargin{:});
            
            idxAlign = p.Results.idxAlign;
            asMatrix = p.Results.asMatrix;
            returnForInvalid = p.Results.returnForInvalid;
            checkInvalid = p.Results.checkInvalid;

            nAlign = pset.nAlign;
            nBases = pset.nBases;
            nConditions = pset.nConditions;

            for i = 1:nOut
                if asMatrix(i)
                    varargout{i} = nan(nBases, nConditions, nAlign);
                else
                    varargout{i} = cell(nBases, nConditions, nAlign);
                end
            end

            outCollect = cell(nOut, 1);
            nargs = nargin(fn);

            for iBasis = 1:nBases
                for iCondition = 1:nConditions
                    for iAlignIdx = 1:length(idxAlign)
                        iAlign = idxAlign(iAlignIdx);
                        
                        if checkInvalid && ~pset.dataValid(iBasis, iCondition, iAlign)
                            if isa(returnForInvalid, 'function_handle')
                                % call errorHandler on this data
                                [outCollect{1:nOut}] = returnForInvalid(time, data, iBasis, iCondition, iAlign);
                            else
                                [outCollect{1:nOut}] = deal(returnForInvalid);
                            end
                            continue;
                        end

                        time = pset.timeData{iBasis, iCondition, iAlign};
                        data = pset.data{iBasis, iCondition, iAlign};

                        if nargs == 5 || nargs < 0
                            [outCollect{1:nOut}] = fn(time, data, iBasis, iCondition, iAlign);
                        else
                            [outCollect{1:nOut}] = fn(time, data);
                        end

                        for iOut = 1:nOut
                            if asMatrix(iOut)
                                varargout{iOut}(iBasis, iCondition, iAlign) = outCollect{iOut};
                            else
                                varargout{iOut}{iBasis, iCondition, iAlign} = outCollect{iOut};
                            end
                        end
                    end
                end
            end

        end
    end

    methods % Time windowing
        function filterTimeWindowByAlign(pset, tWindowByAlign, varargin)
            % Window within time window tWindow = {[tMin tMax], [tMinAlign2 tMaxAlign2]} 
            % by setting tMinActive and tMaxActive to tMin and tMax
            
            if ~iscell(tWindowByAlign)
                tWindowByAlign = {tWindowByAlign};
            end
            
            assert(numel(tWindowByAlign) == pset.nAlign, 'Number of time windows must match nAlign');
            tMinByAlign = cellfun(@(x) x(1), tWindowByAlign);
            tMaxByAlign = cellfun(@(x) x(2), tWindowByAlign);
            
            tMinByTrial = TensorUtils.repmatSliceAlongDims(tMinByAlign, pset.dataSize, pset.DIM_ALIGN);
            tMaxByTrial = TensorUtils.repmatSliceAlongDims(tMaxByAlign, pset.dataSize, pset.DIM_ALIGN);
            
            pset.constrainTimeWindowManual(tMinByTrial, tMaxByTrial)
        end
        
        function constrainTimeWindowManual(pset, tMinByTrial, tMaxByTrial)
            % further constrain the existing .tMinDataManual and .tMaxDataManual
            pset.tMinDataManual = nanmin(pset.tMinDataManual, tMinByTrial);
            pset.tMaxDataManual = nanmin(pset.tMaxDataManual, tMaxByTrial);
        end
        
        function removeTimeWindowManual(pset)
            [pset.tMinDataManual, pset.tMaxDataManual] = deal(nan(pset.dataSize));
        end
        
        function applyTimeWindowData(pset, varargin)
            % permanently truncate data outside of tMinByTrial(...) and
            % tMaxByTrial(...)
            [pset.data, pset.timeData]= pset.getDataTimeWindowed(varargin{:});
        end
        
        function [tMinByTrial tMaxByTrial] = getTimeValidData(pset, varargin)
            % return the largest time window which excludes nan values in psthData
            % by basis x condition x align
            tMinByTrial = pset.alignBasisConditionDataFun(@getTimeValidStart, 'asMatrix', true);
            tMaxByTrial = pset.alignBasisConditionDataFun(@getTimeValidEnd, 'asMatrix', true);
        
            % constrain to lie within manual time windows
            tMinByTrial = nanmax(tMinByTrial, pset.tMinDataManual);
            tMaxByTrial = nanmin(tMaxByTrial, pset.tMaxDataManual);
            
            function t = getTimeValidStart(time, data, iBasis, iCondition, iAlign, varargin)
                validInd = find(~isnan(data), 1, 'first');
                if isempty(validInd)
                    t = NaN;
                else
                    t = time(validInd);
                end
            end

            function t = getTimeValidEnd(time, data, varargin)
                validInd = find(~isnan(data), 1, 'last');
                if isempty(validInd)
                    t = NaN;
                else
                    t = time(validInd);
                end
            end
        end

        function [tMinByTrial tMaxByTrial] = getTimeValidAcrossBasesConditions(pset, varargin)
            % finds the global time valid window for each alignment, 
            % ignoring bases that are entirely nan for the interval
            
            [tMinByTrial tMaxByTrial] = pset.getTimeValidData();

            for iAlign = 1:pset.nAlign 
                tMin = tMinByTrial(:, :, iAlign);
                tMinByTrial(:, :, iAlign) = nanmax(tMin(:));  % nanmax ignores all nan bases
                tMax = tMaxByTrial(:, :, iAlign);
                tMaxByTrial(:, :, iAlign) = nanmin(tMax(:));
            end
        end

        function [tMinByTrial tMaxByTrial] = getTimeValidAcrossConditions(pset, varagin)
            % tMin, tMax are N x C x A
            [tMinByTrial tMaxByTrial] = pset.getTimeValidData();

            % determine the largest valid window that works across all conditions (2)
            tMinAcrossBases = nanmax(tMinByTrial, [], 2);
            tMaxAcrossBases = nanmin(tMaxByTrial, [], 2);

            % re-expand the time windows
            tMinByTrial = repmat(tMinAcrossBases, [1 pset.nConditions 1]);
            tMaxByTrial = repmat(tMaxAcrossBases, [1 pset.nConditions 1]);
            
        end

        function [tMinByTrial tMaxByTrial] = getTimeValidAcrossBases(pset, varargin)
            % tMin, tMax are N x C x A
            [tMinByTrial tMaxByTrial] = pset.getTimeValidData();

            % determine the largest valid window that works across all bases (dim 1)
            tMinAcrossBases = nanmax(tMinByTrial, [], 1);
            tMaxAcrossBases = nanmin(tMaxByTrial, [], 1);

            % re-expand the time windows
            tMinByTrial = repmat(tMinAcrossBases, [pset.nBases 1 1]);
            tMaxByTrial = repmat(tMaxAcrossBases, [pset.nBases 1 1]);
        end

        function [data timeData] = getDataTimeWindowed(pset, varargin)
            % getDataTimeWindowed(tMinByTrial, tMaxByTrial)
            % window the data within tMinByTrial(iBasis, iCondition, iAlign) : tMaxByTrial(b,c,a)
            % by default tMinByTrial and tMaxByTrial are determined by getTimeValidData()
            % if scalar, tM??Data will be expanded to nBases x nConditions x nAlign
            %
            % tMinByTrial and tMaxByTrial are arrays of nBasis x nCondition x nAlign
            % data and timeData are cell arrays of nBasis x nCondition x nAlign
            p = inputParser;
            p.addOptional('tMinByTrial', [], @(x) isnumeric(x));
            p.addOptional('tMaxByTrial', [], @(x) isnumeric(x));
            
            % if tMin / tMax are outside the bounds or data is invalid, 
            % fill with appropriately sized NaNs to populate?
            p.addOptional('fillWithNaN', false, @islogical); 
            
            p.parse(varargin{:});
            tMinByTrial = p.Results.tMinByTrial;
            tMaxByTrial = p.Results.tMaxByTrial;
            fillWithNaN = p.Results.fillWithNaN;

            % by default use the valid window on each basis x condition x align
            if isempty(tMinByTrial) || isempty(tMaxByTrial)
                [tMinTemp tMaxTemp] = getTimeValidData(pset);
                if isempty(tMinByTrial)
                    tMinByTrial = tMinTemp;
                end
                if isempty(tMaxByTrial)
                    tMaxByTrial = tMaxTemp;
                end
            end

            % expand to size in case scalar
            if size(tMinByTrial) == 1 
                tMinByTrial = repmat(tMinByTrial, [pset.nBases, pset.nConditions, pset.nAlign]);
            end
            if size(tMaxByTrial) == 1 
                tMaxByTrial = repmat(tMaxByTrial, [pset.nBases, pset.nConditions, pset.nAlign]);
            end

            [data timeData]= pset.alignBasisConditionDataFun(@getPSTHTimeWindowed, 'checkInvalid', false);
            
            return;

            function [d t] = getPSTHTimeWindowed(time, data, iBasis, iCondition, iAlign)
                tMin = tMinByTrial(iBasis, iCondition, iAlign);
                tMax = tMaxByTrial(iBasis, iCondition, iAlign);
                if isnan(tMin) || isnan(tMax)
                    if fillWithNaN
                        t = NaN;
                        d = NaN;
                    else
                        t = [];
                        d = [];
                    end
                    return;
                end
                ind1 = find(floor(time) == floor(tMin), 1);
                ind2 = find(floor(time) == floor(tMax), 1);
                if isempty(ind1) || isempty(ind2)
                    if fillWithNaN
                        t = makecol(tMin:tMax);
                        d = nan(size(t));
                    else
                        error('Could not find timepoint. Check sampling is consistent and time window is valid');
                    end
                else
                    t = makecol(tMin:tMax);
                    d = makecol(data(ind1:ind2));
                end
            end
        end

        function [data timeData tvecByBasisByAlign] = getDataTimeWindowedValidAcrossConditions(pset, varargin)
            [tMinByTrial tMaxByTrial] = pset.getTimeValidAcrossConditions();

            [data timeData] = pset.getDataTimeWindowed(tMinByTrial, tMaxByTrial);

            % all tvecs in timeData for iBasis, iAlign are the same
            tvecByBasisByAlign = (timeData(:,1,:));
        end

        function [data timeData tvecByAlign] = getDataTimeWindowedValidAcrossBasesConditions(pset, varargin)
            % select the data within the maximum common valid time window for each alignment
            % i.e. ensuring that each basis x condition data has the same time limits within
            % each alignment
            %
            % data and timeData are cell arrays of nBasis x nCondition x nAlign
            % tvecByAlign is a nAlign cell vector with the time vector for each alignment
            
            [tMinByTrial tMaxByTrial] = pset.getTimeValidAcrossBasesConditions();

            [data timeData] = pset.getDataTimeWindowed(tMinByTrial, tMaxByTrial);

            % all tvecs in timeData for iAlign are the same
            tvecByAlign = squeezedim(timeData(1,1,:), [1 2]);
        end

        function [data, timeData, tvecByConditionAlign] = getDataTimeWindowedValidAcrossBases(pset, varargin)
            % window the data such that all bases share a common time window within each condition x align
            % data and timeData are cell arrays of nBasis x nCondition x nAlign
            % tvecByConditionAlign is a nCondition x nAlign cell vector with the time vector for each condition x align

            [tMinByTrial, tMaxByTrial] = getTimeValidAcrossBases(pset, varargin);
            [data timeData] = pset.getDataTimeWindowed(tMinByTrial, tMaxByTrial); 

            % all tvecs across bases are the same
            tvecByConditionAlign = squeezedim(timeData(1,:,:), 1);
        end
    end

    methods % Simple statistics
        function [maxValues maxTimes] = findMaximum(pset, varargin)
            % return nBases x nConditions x nAlign matrices with the maximum value and
            % time of maximum value for each basis x condition x align

            [maxValues maxTimes] = pset.alignBasisConditionDataFun(@findMaxFn, 'asMatrix', [true true]);

            function [valMax timeMax] = findMaxFn(time,data)
                [valMax i] = max(data);
                timeMax = time(i);
            end
        end

        function [minValues minTimes] = findMinimum(pset, varargin)
            % return nBases x nConditions x nAlign matrices with the minimum value and
            % time of minimum value for each basis x condition x align

            [minValues minTimes] = pset.alignBasisConditionDataFun(@findMinFn, 'asMatrix', [true true]);

            function [valMin timeMin] = findMinFn(time,data)
                [valMin i] = min(data);
                timeMin = time(i);
            end
        end

        function maxByBasis = findMaximumAcrossConditions(pset)
            % return a nBases x nAlign vector of maximum values across conditions
            
            maxValues = pset.findMaximum();
            maxByBasis = max(maxValues, [], 2);
        end

        function minByBasis = findMinimumAcrossConditions(pset)
            % return a nBases x nAlign vector of minimum values across conditions
            
            minValues = pset.findMinimum();
            minByBasis = min(minValues, [], 2);
        end
        
        function varData = computeVarianceOverTime(pset)
            % varData is nBases x nConditions x nAlign
            
            varData = pset.alignBasisConditionDataFun(@var, 'asMatrix', true);
        end

        function [ccvByBasisByAlign tVecByBasisByAlign] = crossConditionVariance(pset)
            % compute variance across conditions per condition over time for each basis x align
            % ccvByBasisByAlign is a cell array of nBases x nAlign with a T x 1 vector of ccv over time
            % tVecByAlign has the same size and carries the time vector associated with each ccv vector

            [data timeData tVecByBasisByAlign] = pset.getDataTimeWindowedValidAcrossConditions();

            ccvByBasisByAlign = TensorUtils.mapToSizeFromSubs([pset.nBases pset.nAlign], ...
                'contentsFn', @getCCV, 'asCell', true);

            function ccv = getCCV(iBasis, iAlign)
                % time by conditions matrix for this basis x align
                tByC = cell2mat(squeezedim(data(iBasis, :, iAlign), [1 3]));

                % t by 1 vector of variance across conditions
                ccv = makecol(var(tByC, [], 2));
            end
        end
        
        function [ccvByBasisByAlign tVecByBasisByAlign] = crossConditionStd(pset)
            % compute variance across conditions per condition over time for each basis x align
            % ccvByBasisByAlign is a cell array of nBases x nAlign with a T x 1 vector of ccv over time
            % tVecByAlign has the same size and carries the time vector associated with each ccv vector

            [data timeData tVecByBasisByAlign] = pset.getDataTimeWindowedValidAcrossConditions();

            ccvByBasisByAlign = TensorUtils.mapToSizeFromSubs([pset.nBases pset.nAlign], ...
                'contentsFn', @getCCS, 'asCell', true);

            function ccv = getCCS(iBasis, iAlign)
                % time by conditions matrix for this basis x align
                tByC = cell2mat(squeezedim(data(iBasis, :, iAlign), [1 3]));

                % t by 1 vector of variance across conditions
                ccv = makecol(std(tByC, [], 2));
            end
        end
    end

    methods % Visualization and plotting utilities
        function alignTimeInfo = drawTimeAxisForAlign(pset, iAlign, alignTimeInfo)
            ad = pset.alignDescriptorSet{iAlign};
            if ~exist('alignTimeInfo', 'var')
                alignTimeInfo = pset.alignTimeInfoData(:, :, iAlign);
                emptyMask = cellfun(@isempty, alignTimeInfo);
                alignTimeInfo = alignTimeInfo(~emptyMask);
                alignTimeInfo = cell2mat(makecol(alignTimeInfo(:)));
            end
            %ad.drawTimeAxis(alignTimeInfo);
        end
        
        function drawTimeAxisForConditionAlign(pset, iCondition, iAlign)
            ad = pset.alignDescriptorSet{iAlign};
            alignTimeInfo = pset.alignTimeInfoData(:, iCondition, iAlign);
            alignTimeInfo = cell2mat(makecol(alignTimeInfo(:)));
            ad.drawTimeAxis(alignTimeInfo);
        end

        function plotConditionPanels(pset, varargin)
            % draw all bases superimposed in figure panels by condition
            p = inputParser;
            p.addParamValue('windowed', false, @islogical);
            p.parse(varargin{:});
            windowed = p.Results.windowed;

            if windowed
                [data, time] = pset.getDataTimeWindowed();
            else
                data = pset.data;
                time = pset.timeData;
            end

            % determine how to layout the conditions, use a row if 1-d conditions
            % 2-d if necessary
            nDims = min(2, pset.conditionDescriptor.nAttributes);

            conditionInds = pset.conditionDescriptor.conditionsAsLinearInds;
            conditionAlignsValid = pset.conditionAlignsValidAllBases;
            
            for iAlign = 1:pset.nAlign
                
                % nConditions x 1 logical vector
                conditionsValidThisAlign = conditionAlignsValid(:, iAlign);
                
                if nDims == 1
                    % if there is only one attribute, plot it along one row
                    nRow = 1;

                    % we need a column for every valid condition on this align
                    nCol = nnz(conditionsValidThisAlign);
                    
                    conditionByRowCol = makerow(conditionInds(conditionsValidThisAlign));
                else
                    % one row for each value of the first attribute
                    nRow = pset.conditionDescriptor.nValuesByAttributeGroupBy(1);
                    
                    % reshape conditionsValid into nRow rows and the rest
                    % as columns
                    conditionIndsReshaped = reshape(conditionInds(:), nRow, []);
                    conditionsValidReshaped = reshape(conditionsValidThisAlign, nRow, []);
                    
                    % pare down the rows and columns that have at least one
                    % valid condition in them
                    rowMask = any(conditionsValidReshaped, 2);
                    colMask = any(conditionsValidReshaped, 1);
                    conditionByRowCol = conditionIndsReshaped(rowMask, colMask);
                    
                    nRow = nnz(rowMask);
                    nCol = nnz(colMask);
                end
                
                %fig();
                clf
                p = panel();
                p.pack(nRow, nCol);
                p.margin = 10;

                % build a nice colormap
                cmap = cbrewer('qual', 'Set1', pset.nBases);
                cmap = jet(pset.nBases);
                
                % loop over row and column panels
                yl = nan(pset.nConditions, 2);
                panelHasData = false(nRow, nCol);
                for iCol = 1:nCol
                    for iRow = 1:nRow
                        
                        iCondition = conditionByRowCol(iRow, iCol);
                        if ~conditionsValidThisAlign(iCondition)
                            continue;
                        end
                        
                        h(iCondition) = p(iRow, iCol).select();

                        [tMin, tMax, yMin, yMax] = deal(NaN);
                       
                        for iBasis = 1:pset.nBases
                            dataVec = data{iBasis, iCondition, iAlign};
                            tvec = time{iBasis, iCondition, iAlign};
                            
                            if ~isempty(dataVec) && any(~isnan(dataVec))
                                panelHasData(iRow, iCol) = true;
                            end
                            
                            tMin = min([tMin min(tvec)]);
                            tMax = max([tMax max(tvec)]);
                            yMin = min([yMin min(dataVec)]);
                            yMax = max([yMax max(dataVec)]);
                            
                            plot(tvec, dataVec, '-', 'LineWidth', 2, 'Color', cmap(iBasis,:));
                            hold on
                            
                            yl(iCondition, :) = get(gca, 'YLim');
                        end
                        
                        xlim([tMin tMax]);
                        ylim([yMin yMax]);
                        
                        title(sprintf('%s (%s)', pset.conditionNames{iCondition}, pset.alignNames{iAlign}));
                    end
                end
                
                % draw time axes
                for iCol = 1:nCol
                    for iRow = 1:nRow
                        p(iRow, iCol).select();
                        iCondition = conditionByRowCol(iRow, iCol);
                        if ~panelHasData(iRow, iCol)
                            axis off;
                        else
                            pset.drawTimeAxisForConditionAlign(iCondition, iAlign);
                        end
                    end
                end
                %whitebg(gcf, [0 0 0]);
                p.refresh();
            end
        end

        function plotBasisPanels(pset, varargin)
            p = inputParser;
            p.addParamValue('basisIdx', [1:6], @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            p.parse(varargin{:});

            basisIdx = intersect(p.Results.basisIdx, 1:pset.nBases);
            nBasesPlot = length(basisIdx);

            [data, time] = pset.getDataTimeWindowed();

            for iAlign = 1:pset.nAlign
                fig();
                clf;
                p = panel();
                p.pack(nBasesPlot,1);

                for iBasisIdx = 1:nBasesPlot
                    p(iBasisIdx,1).select();
                    iBasis = basisIdx(iBasisIdx);

                    for iCondition = 1:pset.nConditions
                        timeVec = time{iBasis, iCondition, iAlign};
                        dataVec = data{iBasis, iCondition, iAlign};

                        appear = pset.conditionDescriptor.appearances(iCondition);

                        plot(timeVec, dataVec, ...
                            'Color', appear.color, 'LineWidth', appear.lineWidth);
                        hold on
                    end
                    hold off
                    box off
                    title(sprintf('%s (%s)', pset.basisNames{iBasis}, pset.alignNames{iAlign}));
                    pset.drawTimeAxisForAlign(iAlign);
                    drawnow;
                end
                
                p.margin = 10;
            end
        end

        function plotStateSpace(pset, varargin)
            % plot a 2d or 3d basis1 x basis2 x basis3 trajectory plot
            p = inputParser;
            p.addParamValue('basisIdx', 1:min(pset.nBases, 3), @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            % plot alignments in separate state spaces
            p.addParamValue('separateAlign', false, @islogical);
            p.addParamValue('tMin', [], @isscalar);
            p.addParamValue('tMax', [], @isscalar);
            p.parse(varargin{:});

            separateAlign = p.Results.separateAlign();
            basisIdx = p.Results.basisIdx;
            if length(basisIdx) == 2
                use3d = false;
            elseif length(basisIdx) == 3;
                use3d = true;
            else
                error('Number of bases must be 2 or 3');
            end

            [data, time, tvecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();

            for iAlign = 1:pset.nAlign
                if separateAlign %|| iAlign == 1
                    fig();
                end
                timeVec = tvecByAlign{iAlign};

                for iCondition = 1:pset.nConditions
                    dataVec1 = data{basisIdx(1), iCondition, iAlign};
                    dataVec2 = data{basisIdx(2), iCondition, iAlign};

                    appear = pset.conditionDescriptor.appearances(iCondition);

                    if use3d
                        dataVec3 = data{basisIdx(3), iCondition, iAlign};
                        dataMat = [dataVec1 dataVec2 dataVec3];
                        plot3(dataVec1, dataVec2, dataVec3, ...
                            'Color', appear.color, 'LineWidth', appear.lineWidth);
                    else
                        dataMat = [dataVec1 dataVec2];
                        plot(dataVec1, dataVec2, ...
                            'Color', appear.color, 'LineWidth', appear.lineWidth);
                    end

                    hold on
                    ti = cell2mat(pset.alignTimeInfoData(:, iCondition, iAlign));
                    
                    pset.alignDescriptorSet{iAlign}.drawOnData({ti}, {timeVec}, {dataMat}, ...
                        'drawLegend', iCondition == 1);
                end

                if separateAlign || iAlign == pset.nAlign
                    hold off
                    box off
                    xlabel(pset.basisNames{basisIdx(1)});
                    ylabel(pset.basisNames{basisIdx(2)});
                   
                    if separateAlign
                        title(pset.alignNames{iAlign});
                    end
                     if use3d
                        zlabel(pset.basisNames{basisIdx(3)});
                        view([-40 20]);
                    end
                    
                    axis tight
                    axis square
                    axis vis3d
                end
            end
        end
    end

    methods % Build data matrices
        function [out tvecByConditionAlign] = buildCTAByN(pset, varargin)
            % out is C*T*A x N concatenated matrices for each alignment
            % timeVec: nAlign cell of time vectors common to alignment

            p = inputParser();
            p.addParamValue('timeValidAcrossConditions', false, @islogical)
            p.parse(varargin{:});


            % data is N bases x nConditions x nAlign cell array with vectors of length t
           if p.Results.timeValidAcrossConditions
                [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
                tvecByConditionAlign = repmat(timeVecByAlign', pset.nConditions, 1);
            else
                [data, timeData, tvecByConditionAlign] = pset.getDataTimeWindowedValidAcrossBases();
            end
            
            % convert to (nConditions*nAlign) x nBases cell of column vectors over time
            data = reshape(data, [pset.nBases, pset.nConditions * pset.nAlign])';

            % convert to nAlign*nConditions*T x Nbases matrix`
            out = cell2mat(data);
        end

        function [out tvecByConditionAlign] = buildCTByNEachA(pset, varargin)
            % out: A cell vector of CT x N matrices 
            % timeVec: nAlign cell of time vectors common to alignment
            
            p = inputParser();
            p.addParamValue('timeValidAcrossConditions', false, @islogical)
            p.parse(varargin{:});

            % data is N bases x nConditions x nAlign cell array with vectors of length t
            % allow different time vectors per condition
            if p.Results.timeValidAcrossConditions
                [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
                tvecByConditionAlign = repmat(timeVecByAlign', pset.nConditions, 1);
            else
                [data, timeData, tvecByConditionAlign] = pset.getDataTimeWindowedValidAcrossBases();
            end

            % convert to nAlign cell of nCondition*time x nBases 
            out = cell(pset.nAlign, 1);
            for iAlign = 1:pset.nAlign
                out{iAlign} = cell2mat(data(:,:,iAlign)');
            end
        end

        function [out tvecByConditionAlign] = buildTByNEachCA(pset, varargin)
            % out: C x A cell of N x T matrices 
            % tvecByConditionAlign is C x A cell of time vectors common to condition x alignment

            p = inputParser();
            p.addParamValue('timeValidAcrossConditions', false, @islogical)
            p.parse(varargin{:});
                        
            % data is N bases x nConditions x nAlign cell array with vectors of length t
            if p.Results.timeValidAcrossConditions
                [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
                tvecByConditionAlign = repmat(makerow(timeVecByAlign), pset.nConditions, 1);
            else
                [data, timeData, tvecByConditionAlign] = pset.getDataTimeWindowedValidAcrossBases();
            end

            out = cell(pset.nConditions, pset.nAlign);
            for iAlign = 1:pset.nAlign
                for iCondition = 1:pset.nConditions

                    % nBases cell of time vectors
                    dataThisCA = data(:,iCondition, iAlign);

                    % nAlign*nTime x nBases matrix for this condition
                    out{iCondition, iAlign} = cell2mat(dataThisCA');
                end
            end
        end
        
        function [out] = buildTAByCByN(pset, varargin)
            [data timeData timeVecByAlign] = pset.getDataTimeWindowedValidAcrossBasesConditions();
            
            % data is nBases x nConditions x nAlign cell of nTime(iAlign) x 1
            % reshape it to be nAlign x nConditions x nBases
            data = permute(data, [3 2 1]);
            % then flatten into cell
            out = cell2mat(data);
        end
        
        function out = buildNByTAByAttributeTensor(pset, varargin)
            % build a tensor of N by T by nValsAttr1 by nValsAttr2 x ...
            % this tensor is the format used by dpca_covs
            taByCByN = pset.buildTAByCByN();
            nTime = size(taByCByN, 1);
            condSize = pset.conditionDescriptor.conditionsSize;
            nByTAByC = permute(taByCByN, [3 1 2]);
            out = reshape(nByTAByC, [pset.nBases nTime makerow(condSize)]);
        end
    end

    methods % Comparative statistics
        function [distByFromAlign timeVecByFromAlign] = getDistanceBetween(pset, cFromList, cToList, varargin)
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
            [tByNEachCA timeVecByAlign] = pset.buildTByNEachCA('timeValidAcrossConditions', true);
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

        function [distByFromAlign timeVecByFromAlign cFromList cToList] = getDistanceAlongComparisonAxis(pset, compareAcross, varargin)
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
            [cFromList cToList] = deal(zeros(nCompare, 1));
            for iCompare = 1:nCompare
                idxThisComparison = idxCompare{iCompare};
                assert(length(idxThisComparison) == 2, 'Comparison axis must span exactly two elements');
                if ~reverse
                    idxThisComparison = idxThisComparison([2 1]);
                end
                cFromList(iCompare) = idxThisComparison(1);
                cToList(iCompare) = idxThisComparison(2);
            end

            [distByFromAlign timeVecByFromAlign] = ...
                pset.getDistanceBetween(cFromList, cToList, p.Unmatched);
        end
    end

end
