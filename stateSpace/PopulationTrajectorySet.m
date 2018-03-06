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

    % metadata properties that won't affect anything once set
    properties
        datasetName

        % Set this to true to keep all comute-on-demand properties when
        % using .saveFast or caching via CacheManger
        % Setting it to false will save space on disk
        keepComputedOnSaveFast = false;
    end
    
    % Properties which determine whether specific other properties are computed or stored manually
    properties(SetAccess=?PopulationTrajectorySetBuilder)
        % There are currently four modes in which a pset can be operating:
        % 1. auto: fully automatic
        % 2. dataByTrialManualmanually specified data by trial, but with data sources. For example, this would be used after splicing.
        %      It enables new condition groupings to be used but not new alignments.
        % 3. manually specified data by trial, manual trialLists, without data sources. This allows dataMean to be
        %      computed on demand, as well as randomizations to be done. But no new condition groupings nor new alignments.
        % 4. manual data mean, without data by trial.
        %
        % For each of thse modes, the properties below are set as follows:

        % Mode | dataByTrialManual | dataInfoManual | dataMeanManual
        % -----|-------------------|------------------|---------------
        % auto    | false             | false            | false
        % autoWithManualDataByTrial    | true              | false            | false
        % auto   | true              | true             | false
        % 4    | true             | true             | true
        %
        % Note that dataMeanManual is really just ~hasDataByTrial.

        % Are data sources included with this pset? If false, values are computed dynamically
        % from the dataSources and stored in the .odc. If true, values are derived from *Manual properties
        dataInfoManual = false

        % is dataByTrial and associated properties being extracted directly from the dataSources (false)
        % or stored manually in the corresponding dataByTrialManual properties.
        % this is useful if the data of interest was computed manually by the user on a single trial basis in correspondence
        % with trial data objects, so that subsequent grouping operations can take
        dataByTrialManual = false
        
        % is dataMean computed from dataByTrial? This must be true if dataByTrial is not present or stored with the pset. This
        % would be the case if the means are computed in some way, e.g. via a projection or concatenation operations
        dataMeanManual = false;
    end

    methods % ensure consistency for the data manual properties above
        function str = describeDataManualMode(pset)
            % generate a short description that encapuslates the settings of the manual mode properties above
            if pset.dataInfoManual
                if ~pset.dataByTrialManual
                    % mode 4
                    str = 'manual trial-averaged, no single trial';
                else
                    % mode 3
                    str = 'manual single trial, manual trial groupings';
                end
            else
                if pset.dataByTrialManual
                    % mode 2
                    str = 'manual single trial, trial groupings via data sources';
                else
                    % mode 1
                    str = 'single trial from data sources';
                end
            end
        end
    end

    % properties which control the behavior of the pset and will invalidate
    % computed values when they are set by one of the corresponding setProperty methods
    properties(SetAccess=?PopulationTrajectorySetBuilder)

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

        % Ignore all zero spike trials when building trial averages
        ignoreAllZeroSpikeTrials = false;

        % Ignore zero spike trials at the beginnning and end of the list of trials
        % when building trial averages (i.e. but not in between the other
        % trials)
        ignoreLeadingTrailingZeroSpikeTrials = false;

        % alignDescriptors, conditionDescriptor, and translationNormalization
        % which align, group, and translate/normalize the data generated from
        % the dataSources. These may be set using eponymous methods prefixed with set*

        % nAlign x 1 cell of alignDescriptors
        alignDescriptorSet = {};

        % ConditionDescriptor instance describing condition information
        conditionDescriptor

        % StateSpaceTranslationNormalization instance describing the
        % translation and normalization to apply to each basis.
        % This will be applied during buildDataByTrial for dataInfoManual
        % = false psets (thus being reflected in the trial-averages automatically)
        % or applied manually to dataByTrial (if non-empty) and dataMean
        % (if non-empty) for dataInfoManual = false
        translationNormalization

        % SpikeFilter instance to use when converting spiking units to
        % firing rate channels
        spikeFilter

        % These properties store raw data sources (TrialDataConditionAlign instances)
        % from which data is extracted, as well as track from where each basis
        % originates. Along with other misc settings
        timeUnitName % string name of common time units

        timeUnitsPerSecond % scalar conversion factor

        % TrialData data sources which source all data for the trajectories
        % this may be a single trial data object or many. If there is only
        % one, all bases are considered simultaneous.
        dataSources = {}

        % nBases x 1 index into dataSourceSet.
        basisDataSourceIdx

        % nBases x 1 cellstr indicating which channel name to extract data
        % from
        basisDataSourceChannelNames

        % random seed used as initial seed when generating random data
        % sets
        randomSeed

        % these have get methods to populate with default vectors if the property is empty
        % nBases x 1 logical
        basisValidPermanent % stores permanent invalid data

        % nBases x 1 cellstrvec
        basisInvalidCausePermanent % causes of permanent invalid

        % nBases x 1 cellstrvec
        basisValidTemporary % temporary invalid data

        % nBases x 1 cellstrvec
        basisInvalidCauseTemporary

        % holds randomized data (will be walked by propMeta)
        randomized
        
        % holds stored data (will be walked by propMeta)
        stored

        % stores manual values of properties that are persistent (as opposed to ODC)
        manual

            % some things that may be stored within stored
            % FOR INDIVIDUAL TRIAL DATA, BUT LINKED TO THE CONDITIONS SINCE ONLY VALID TRIALS ARE INCLUDED
            % these are computed from dataByTrial if dataInfoManual is false,
            % else stored in the corresponding *Manual property

            % nBases x nAlign cell each containing nTrials x nTime analog data
            % for that bases for tha alignment in order by trial number. These
            % will be guaranteed to share the same time vector across bases,
            % unlike dataByTrial. This time vector will be given by pset.tvecDataMean
            % Unlike dataByTrial, this will depend on the ConditionDescriptor
            % as well for the grouping
%             dataByTrialCommonTimeGrouped

            % nAlign x 1 cell of N x T x C x R
            % for a pset with single trial data, this is of the same form that
            % is returned by arrangeNbyTAbyConditionAttributesbyTrials. The
            % rationale for being able to cache it here is that when projecting
            % we can compute pseudo-single trials through the projection so that
            % arrangeNbyTAbyConditionAttributesbyTrials can sample this tensor if
            % no actual single trial data is available
            % dataCachedSampledTrials

            % N x C number of trials in dataCachedSampledTrials that
            % contain usable trials
            % dataCachedSampledTrialCounts

            % nAlign x 1 cell of N x T x C x nTrials
            % the mean for each basis when the trial sampled in
            % dataCachedSampledTrials is excluded from the mean
            % this is mainly used for cross-validation
            % dataCachedMeanExcludingSampledTrials

            % ConditionDescriptor for dataMeanRandomized describing the
            % randomization method used
            %         conditionDescriptorRandomized

        % nAlign x 1 cell with nBases x nCondition x nTime(iAlign) x nRandomSamples numeric
        % tensors containing randomly generated dataMean
%         dataMeanRandomized
%         dataSemRandomized
%         dataNumTrialsRawRandomized % same as dataNumTrialsRaw

        % nAlign x 1 cell of nBases x T x nConditions x nRandomSamples numeric array
        % containing randomized samples of differences of pairs of trials,
        % scaled by 1/sqrt(2*nTrials). These serve as a sample from the
        % centered distribution of the estimation noise on the dataMean
        % traces. These are used by StateSpaceProjectionStatistics in
        % estimating a noise floor (buildDifferenceOfTrials)
        % same as dataDifferenceOfTrialsScaledNoiseEstimate except for
%         dataDifferenceOfTrialsScaledNoiseEstimateRandomized
    end

    % ON DEMAND PROPERTIES
    % Properties whose values are computed on-demand and persist within odc or .manual
    % depending on dataInfoManual or dataByTrialManual. They have special get and set
    % methods that use the PropMeta metadata to be built on the fly and invalidated when appropriate
    % They require write-access by PopulationTrajectorySetBuilder
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
        basisNames

        % nBases x 1 cellstr: units for each basis
        basisUnits

        %% BELOW ARE FOR INDIVIDUAL TRIAL DATA
        % these are computed from dataSources if dataByTrialManual is false,
        % else stored in the corresponding *Manual property

        % data by trial cells contain common time vectors across trials,
        % which differ across bases and aligns

        % nBases x nAlign cell containing ordered data by trial
        % each cell contains nTrials x nTime analog data for that basis,
        % for that alignment, in order by trial number. The time vector
        % along the columns is given by
        % nanmin(tMinByTrial{..}) : nanmax(tMaxByTrial{..})
        % raw implies does not reflect .basisValid
        dataByTrialRaw

        % nBases x nAlign numeric matrix containing the start and stop
        % times for the time vector indicating the time along the columns
        % of the corresponding cell of dataByTrial
        % raw implies does not reflect .basisValid
        tMinForDataByTrialRaw
        tMaxForDataByTrialRaw

        % nBasis x nAlign cell arrays indicating the start and stop
        % timepoints for each trial in dataByTrial
        % raw implies does not reflect .basisValid
        tMinValidByTrialRaw
        tMaxValidByTrialRaw

        %% BELOW ARE FOR TRIAL-AVERAGED DATA WITHIN CONDITION
        % these are computed from dataByTrial if dataInfoManual is false,
        % else stored in the corresponding *Manual property

        % nBases x nConditions cell of trial (into dataByTrial) for each
        % condition (buildTrialLists). This will not respect basisValid, as it will
        % be used to determine basisValid
        trialListsRaw

        % trial averaged data contains common time vectors across
        % conditions and bases, but which differ across aligns. Invalid
        % portions of the bases with smaller windows will be filled with
        % NaNs.

        % containing the start and stop times for
        % which sufficient trials exist to compute a trial-average for each
        % align, basis, and conditions
        tMinValidByAlignBasisConditionRaw
        tMaxValidByAlignBasisConditionRaw

        % nBases x 1 logical: which bases are considered included in the
        % present analyses and time-vector computations. updating this will
        % invalidate data means. it is added to make basis filtering
        % operations simpler (since all psets can maintain the same logical
        % vector
        basisValid
        basisInvalidCause

        % properties below will respect basisValid:

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

        % nBases x nCondition scalar array indicating how many
        % trials contributed to data in each cell  (computed on demand in get method)
        % this has to be an ODC property as sometimes it depends on trialListsRaw, but it can
        % also be manual when transformations of data are computed
        dataNumTrials

        % nBases x nCondition logical array indicating whether there is
        % valid data in the corresponding cell
        dataMeanValid

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the lower confidence interval value for
        % dataMean (buildDataRandomizedIntervals)
%         dataIntervalLow

        % nAlign x 1 cell with nBases x nConditions x nTimeDataMean(iAlign)
        % numeric array containing the upper confidence interval value for
        % dataMean (buildDataRandomizedIntervals)
%         dataIntervalHigh

        % nAlign x 1 cell of nBases x T x nConditions x nSamples numeric array
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
        % Common units across basisUnits or ''
        dataUnitsCommon

        % numeric scalar indicating spacing between successive time points
        % all channels timeseries will be interpolated to a time vector
        % with this spacing
        timeDelta

        % nAlign x 1 numeric vectors indicating the number of timepoints
        % used for each alignment in the dataMean cell
        nTimeDataMean

        % nAlign x 1 cell array of time vectors corresponding to the 3rd
        % dimension of dataMean{iAlign}
        tvecDataMean

        % nAlign x 1 cell array of time vectors which are automatically
        % offset from tvecDataMean for plotting
        tvecDataMeanForPlotting

        % nAlign x 1 vector of offsets to where zero of align i will be
        % plotted (the horizontal offset for each align's time zero)
        offsetsTimeDataMeanForPlotting

        % nAlign x nBases cell of time vectors for each .dataByTrial matrix
        tvecDataByTrialRaw
        tvecDataByTrial

        % nAlign x nBases x nConditions numeric of number of non-nan
        % timepoints in each data mean
        nTimeValidByAlignBasisCondition

        % is dataByTrial empty?
        hasDataByTrial

        % same as dataByTrial except masked by basisValid along first
        % dimension
%         dataByTrialValidOnly

        % nBases cell of nTrials x 1 logical masks
        trialHasSpikesMaskByBasis

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

        % nAlign x nBases x nConditions logical indicating whether that
        % basis provides a valid trial average
        hasValidTrialAverageByBasisCondition

        % nConditions x 1 logical vector indicating whether the condition
        % has trial average data for EVERY basis on ALL aligns
        conditionHasValidTrialAverageAllAlignsBases

        % nConditions x 1 logical: conditions not listed valid will be
        % ignored (and thus NaN) in .dataMean when it is computed
        % this can be used if some bases but not all bases have data for
        % that condition, and the condition is messing up the time vector
        % limits
        conditionIncludeMask

        % nAlign x nConditions matrix containing the start and stop times
        % for which sufficient trials exist to compute a trial-average for
        % ALL valid bases, on each align and condition. only Conditions for
        % which conditionHasValidTrialAverageAllAlignsBases is true are
        % considered. Only bases which are valid are considered
        tMinValidAllBasesByAlignCondition
        tMaxValidAllBasesByAlignCondition

        % nBases x 1 cell of conditions associated with each trial, as in
        % dataByTrial
        conditionIdxByTrial

        % indicates how the elements in data are computed from individual trials
        % true means all trials correspond one-to-one with each other
        % false means trials were not simultaneous, implying that only
        % trial-averages should be compared
        simultaneous

        nTrials % only valid when simultaneous is true

        nTrialsValid % only valid when simultaneous is true

        nTrialsByCondition % only valid when simultaneous is true

        % number of data sources used across all bases
        nDataSources

        nBases % number of bases (e.g. units, analog channels)

        nBasesValid % number of bases currently marked valid

        nBasesValidPermanent
        nBasesPermanentlyInvalid

        nBasesTemporarilyInvalid

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

        % nCondition x 1 logical indicating which condition/alignments
        % have any trials for any bases
        conditionsWithTrials

        % nAlign x nBases x nConditions logical indicating which bases have valid
        % trial averages on each condition
        %         alignBasisConditionsWithValidTrialAverage

        % nAlign x nConditions logical indicating which align/conditions
        % have trial averages on at least one basis
        alignConditionsWithTrialAverage

        % nBases x nConditions logical indicating which bases have valid
        % trial averages on each condition on ALL aligns
        basisConditionsWithValidTrialAverage

        % nBases x 1 logical array indicating bases which have no valid trials on some condition for which other bases have trials
        basesMissingTrialAverageForNonEmptyConditionAligns

        basisValidWithTrialAverageAllNonEmptyConditionAligns

        % nBases x 1 logical indicating which valid basess have at least one valid
        % trial average on some condition, on some alignment
        basesNonEmpty

        % nConditions x 1 logical indicating which conditions are present
        % in all nonEmptyBases
        conditionsWithValidTrialAverageOnNonEmptyBases
    end


    properties(SetAccess=protected, Hidden)
        % odc is an instance of OnDemandCache
        % which is a handle class. Properties which derive from
        odc
    end

    properties(SetAccess=protected, Hidden)
        propMeta = PopulationTrajectorySet.buildPropertyMeta();
    end

    properties(Hidden)
        debugPropMeta = true;
    end

    % Constructor
    methods
        function pset = PopulationTrajectorySet()
            pset = pset.initialize();
        end
    end

    methods(Static, Hidden)
        function meta = buildPropertyMeta()
            meta = struct();

            % settings
            meta.datasetName = PropertyShapeMeta({}, 'char', 'group', 'settings');
            meta.timeUnitName = PropertyShapeMeta({}, 'char', 'group', 'settings');
            meta.timeUnitsPerSecond = PropertyShapeMeta({}, 'char', 'group', 'settings');
            meta.spikeFilter = PropertyShapeMeta({}, 'SpikeFilter', 'group', 'settings');
            meta.minTrialsForTrialAveraging = PropertyShapeMeta({}, 'int', 'group', 'settings');
            meta.minFractionTrialsForTrialAveraging = PropertyShapeMeta({}, 'float', 'group', 'settings');
            meta.ignoreAllZeroSpikeTrials = PropertyShapeMeta({}, 'logical', 'group', 'settings');
            meta.ignoreLeadingTrailingZeroSpikeTrials = PropertyShapeMeta({}, 'logical', 'group', 'settings');

            % descriptors
            meta.alignDescriptorSet = PropertyShapeMeta({'A'}, 'AlignDescriptor', 'group', 'descriptors');
            meta.conditionDescriptor = PropertyShapeMeta({}, 'AlignDescriptor', 'group', 'descriptors');
            meta.translationNormalization = PropertyShapeMeta({}, 'StateSpaceTranslationNormalization', ...
              'group', 'descriptors', 'emptyOkay', true, ...
              'customDims', {'N'}, 'customSelectAlongDimFn', @trNormSelectBasesFn);

           function trNorm = trNormSelectBasesFn(meta, trNorm, dimNames, masksByDim, varargin) %#ok<INUSL>
               assert(strcmp(dimNames{1}, 'N'));
               if ~isempty(trNorm)
                   trNorm = trNorm.selectBases(masksByDim{1});
               end
           end
          
            % basis info
            meta.basisNames = PropertyShapeMeta({'N'}, 'char', 'group', 'basisInfo', ...
              'odc', true, 'useManualProp', 'dataInfoManual', 'depends', {'dataSources'}, ...
              'buildFn', 'buildBasisNamesUnits');

            meta.basisUnits = PropertyShapeMeta({'N'}, 'char', 'group', 'basisInfo', ...
              'odc', true, 'useManualProp', 'dataInfoManual', 'depends', {'dataSources'}, ...
              'buildFn', 'buildBasisNamesUnits');

            meta.basisValidPermanent = PropertyShapeMeta({'N'}, 'logical', 'group', 'basisInfo');
            meta.basisInvalidCausePermanent = PropertyShapeMeta({'N'}, 'logical', 'group', 'basisInfo');
            meta.basisValidTemporary = PropertyShapeMeta({'N'}, 'logical', 'group', 'basisInfo');
            meta.basisInvalidCauseTemporary = PropertyShapeMeta({'N'}, 'logical', 'group', 'basisInfo');

            % data source info
            meta.dataSources = PropertyShapeMeta({'nDataSources'}, 'TrialDataConditionAlign', 'group', 'dataSourceInfo', ...
                'customDims', {'A'}, 'customSelectAlongDimFn', @dataSourcesSelectAlongDimFn);

           function dataSources = dataSourcesSelectAlongDimFn(meta, dataSources, dimNames, masksByDim, varargin) %#ok<INUSL>
               for iD = 1:numel(dimNames)
                   mask = masksByDim{iD};
                   switch dimNames{iD}
                       case 'A'
                           for iDS = 1:numel(dataSources)
                               dataSources{iDS} = dataSources{iDS}.selectAlign(mask);
                           end
                   end
               end
           end
           
            meta.basisDataSourceIdx = PropertyShapeMeta({'N'}, 'int', 'group', 'dataSourceInfo');
            meta.basisDataSourceChannelNames = PropertyShapeMeta({'N'}, 'char', 'group', 'dataSourceInfo');

            % single trial
            meta.dataByTrialRaw = PropertyShapeMeta({{'N', 'A'}, {'R', 'Tsingle'}}, 'float', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataByTrialManual', ...
              'depends', {'dataSources', 'alignDescriptorSet', 'conditionDescriptor'}, ...
              'buildFn', 'buildDataByTrial', 'translate', true, 'normalize', true);

            meta.tMinForDataByTrialRaw = PropertyShapeMeta({'N', 'A'}, 'float', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataByTrialManual', ...
              'depends', {'dataSources', 'alignDescriptorSet', 'conditionDescriptor'}, ...
              'buildFn', 'buildDataByTrial');

            meta.tMaxForDataByTrialRaw = PropertyShapeMeta({'N', 'A'}, 'float', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataByTrialManual', ...
              'depends', {'dataSources', 'alignDescriptorSet', 'conditionDescriptor'}, ...
              'buildFn', 'buildDataByTrial');

            meta.tMinValidByTrialRaw = PropertyShapeMeta({{'N', 'A'}, {'R'}}, 'float', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataByTrialManual', 'depends', 'dataByTrialRaw', ...
              'buildFn', 'buildDataByTrialPerTrialLimits');

            meta.tMaxValidByTrialRaw = PropertyShapeMeta({{'N', 'A'}, {'R'}}, 'float', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataByTrialManual', 'depends', 'dataByTrialRaw', ...
              'buildFn', 'buildDataByTrialPerTrialLimits');

            meta.trialListsRaw = PropertyShapeMeta({{'N', 'C'}, {'Rc'}}, 'int', 'group', 'singleTrial', ...
              'odc', true, 'useManualProp', 'dataInfoManual', ...
              'depends', {'dataSources', 'conditionDescriptor', 'trialHasSpikesMaskByBasis', 'ignoreAllZeroSpikeTrials', 'ignoreLeadingTrailingZeroSpikeTrials'}, ...
              'buildFn', 'buildTrialLists');

            % automatically computed always, since they are determined directly from tM**ValidByTrialRaw
            meta.tMinValidByAlignBasisConditionRaw = PropertyShapeMeta({'A', 'N', 'C'}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'depends', {'alignDescriptorSet', 'conditionDescriptor', 'tMinValidByTrialRaw', ...
              'tMaxValidByTrialRaw', 'minFractionTrialsForTrialAveraging', 'minTrialsForTrialAveraging'}, ...
              'buildFn', 'buildTimeWindowsByAlignBasisCondition');

            meta.tMaxValidByAlignBasisConditionRaw = PropertyShapeMeta({'A', 'N', 'C'}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'depends', {'alignDescriptorSet', 'conditionDescriptor', 'tMinValidByTrialRaw', ...
              'tMaxValidByTrialRaw', 'minFractionTrialsForTrialAveraging', 'minTrialsForTrialAveraging'}, ...
              'buildFn', 'buildTimeWindowsByAlignBasisCondition');

            % automatically computed from tM**ValidByAlignBasisConditionRaw and basisValidPermanent and basisValidTemporary
            meta.basisValid = PropertyShapeMeta({'N'}, 'logical', 'group', 'trialAverage', ...
              'odc', true, 'depends', {'basisValidPermanent', 'basisValidTemporary', 'tMinValidByAlignBasisConditionRaw', ...
              'tMaxValidByAlignBasisConditionRaw', 'conditionDescriptor', 'alignDescriptorSet'}, ...
              'buildFn', 'buildBasisValid');

            meta.basisInvalidCause = PropertyShapeMeta({'N'}, 'char', 'group', 'trialAverage', ...
              'odc', true, 'depends', {'basisValidPermanent', 'basisValidTemporary', 'basisInvalidCauseTemporary', ...
              'basisInvalidCausePermanent', 'tMinValidByAlignBasisConditionRaw', 'tMaxValidByAlignBasisConditionRaw', 'conditionDescriptor', 'alignDescriptorSet'}, ...
              'buildFn', 'buildBasisValid');

            meta.tMinForDataMean = PropertyShapeMeta({'A'}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataMeanManual', 'depends', {'alignDescriptorSet', 'conditionDescriptor', ...
              'basisValid', 'tMinValidByAlignBasisConditionRaw', 'tMaxValidByAlignBasisConditionRaw', 'spikeFilter'}, ...
              'buildFn', 'buildTvecDataMean');

            meta.tMaxForDataMean = PropertyShapeMeta({'A'}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataMeanManual', 'depends', {'alignDescriptorSet', 'conditionDescriptor', ...
              'basisValid', 'tMinValidByAlignBasisConditionRaw', 'tMaxValidByAlignBasisConditionRaw', 'spikeFilter'}, ...
              'buildFn', 'buildTvecDataMean');

            meta.dataMean = PropertyShapeMeta({{'A'}, {'N', 'C', 'Tmean'}}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataMeanManual', 'translate', true, 'normalize', true, 'buildFn', 'buildDataMean', ...
              'depends', {'tMinForDataMean', 'tMaxForDataMean', 'spikeFilter', 'basisValid'});

            meta.dataSem = PropertyShapeMeta({{'A'}, {'N', 'C', 'Tmean'}}, 'float', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataMeanManual', 'translate', false, 'normalize', true, 'buildFn', 'buildDataMean', ...
              'depends', {'tMinForDataMean', 'tMaxForDataMean', 'spikeFilter', 'basisValid'});

            meta.dataNumTrials = PropertyShapeMeta({'N', 'C'}, 'int', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataMeanManual', 'buildFn', 'buildDataNumTrials', ...
              'depends', {'trialListsRaw', 'basisValid'});

            % align summary aggregation
            meta.alignSummaryData = PropertyShapeMeta({'nDataSources', 'A'}, 'AlignSummary', 'group', 'trialAverage', ... % nDataSources == nAlignSummaryData
              'odc', true, 'useManualProp', 'dataInfoManual', 'depends', 'alignDescriptorSet', 'buildFn', 'buildAlignSummaryData', ...
              'customDims', {'C'});

            meta.alignSummaryAggregated = PropertyShapeMeta({'A'}, 'AlignSummary', 'group', 'trialAverage', ...
              'odc', true, 'depends', 'alignSummaryData', 'buildFn', 'buildAlignSummaryAggregated', ...
              'customDims', {'C'});

            meta.basisAlignSummaryLookup = PropertyShapeMeta({'N'}, 'int', 'group', 'trialAverage', ...
              'odc', true, 'useManualProp', 'dataInfoManual', 'depends', 'alignDescriptorSet', 'buildFn', 'buildAlignSummaryData');

            % noise estimate
            meta.dataDifferenceOfTrialsScaledNoiseEstimate = PropertyShapeMeta({{'A'}, {'N', 'Tmean', 'C', 'R'}}, 'float', ...
              'group', 'noise', 'emptyOkay', true, ...
              'useManualProp', 'dataMeanManual', ...
              'odc', true, 'translate', true, 'normalize', true, 'buildFn', 'buildDataNoiseEstimate', ...
              'depends', {'dataByTrial', 'trialLists'});
          
            meta.stored = struct();
            meta.randomized = struct();
        end

        function metaFiltered = getPropMetaInGroup(group)
            meta = PopulationTrajectorySet.buildPropertyMeta();
            if ~iscell(group), group = {group}; end
            metaFiltered = TrialDataUtilities.Struct.filterFields(meta, @(m, prop) ~isstruct(m) && ismember(m.getAttrWithDefault('group', ''), group));
        end

        function [lookup, vec, keepIdx] = filterUsedUpdateLookup(lookup, vec)
            % lookup is a vector of indices into vec
            % keep elements of vec that are mentioned by lookup
            % and then update lookup to reflect the new positions of the values in vec
            N = numel(vec);
            oldVecIdx = (1:N)';
            [keepIdx, lookup] = intersect(lookup, oldVecIdx);
            vec = vec(keepIdx);
        end
    end

    % Related to on demand computed and PropMeta properties (including .stored)
    methods(Hidden)
        function useManual = isPropValueStoredInManual(pset, prop)
            propMeta = pset.propMeta.(prop);
            if isfield(propMeta.attr, 'useManualProp')
                %  is determined by named property value
                useManual = pset.(propMeta.attr.useManualProp);
            elseif propMeta.getAttrWithDefault('useManualIfSet', false)
                % manual is used if not empty
                useManual = isfield(pset.manual, prop) && ~isempty(pset.manual.(prop));
            else
                % if its not odc, then it's not using .manual
                useManual = false;
            end
        end
        
        function isODC = isPropODC(pset, prop)
            propMeta = pset.propMeta.(prop);
            isODC = propMeta.getAttrWithDefault('odc', false);
        end

%         function tf = propCanBeSetCurrently(pset, prop)
%             % determines whether this property value can be specified externally,
%             % this is true if the property is not odc, or if its odc but using .manual
%             if ~pset.isPropODC(prop)
%                 % non-ODC properties can always be assigned
%                 tf = true;
%             else
%                 % ODC properties can be assigned when stored in manual
%                 tf = pset.isPropValueStoredInManual(prop);
%             end
%         end

        % retrieves the value either from .manual or from .odc, and if not present in .odc, runs buildFn to compute it
        function value = processPropGet(pset, prop)
            propMeta = pset.propMeta.(prop);
            if ~pset.isPropODC(prop)
                error('Get for property %s should not use processODCPropGet', prop);
            end
            useManual = pset.isPropValueStoredInManual(prop);
            if useManual
                value = pset.manual.(prop);
            else
                % check for value already pre-computed ODC
                if isfield(pset.odc.data, prop) && ~isempty(pset.odc.data.(prop))
                    % non-empty value found in ODC, return it
                    value = pset.odc.data.(prop);
                else
                    % call the build fn which is expected to set the value in odc, then return the result from odc
                    buildFn = propMeta.attr.buildFn;
                    if pset.debugPropMeta
                        debug('PropMeta: Computing %s on demand via %s\n', prop, buildFn);
                    end
                    pset.(buildFn)();
                    value = pset.odc.data.(prop);
                end
            end
        end

        function [value, storedManual, isODC] = processPropGetWithoutComputing(pset, prop)
            isODC = pset.isPropODC(prop);
            storedManual = pset.isPropValueStoredInManual(prop);

            if ~isODC
                value = pset.(prop);
            elseif storedManual
                value = pset.manual.(prop);
            else
                % check for value already pre-computed ODC
                if isfield(pset.odc.data, prop) && ~isempty(pset.odc.data.(prop))
                    % non-empty value found in ODC, return it
                    value = pset.odc.data.(prop);
                else
                    value = [];
                end
            end
        end

        % stores the value in .manual or in .odc. DOES NOT INVALIDATE DEPENDENT PROPERTIES.
        % this is important as it makes it easier to assign multiple properties at once without worrying about complex sequencing
        function pset = processODCPropSet(pset, prop, value)
            pset.warnIfNoArgOut(nargout);
            if ~pset.isPropODC(prop)
                error('Set for property %s should not use processODCPropSet', prop);
            end
            useManual = pset.isPropValueStoredInManual(prop);
            if useManual
                pset.manual.(prop) = value;
            else
                pset.odc.data.(prop) = value;
            end
        end

        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function pset = invalidateCache(pset, varargin)
            pset.warnIfNoArgOut(nargout);
            pset = pset.invalidateProperties(fieldnames(pset.propMeta));
        end

        function pset = invalidateProperties(pset, props)
            if ischar(props)
                props = {props};
            end
            for i = 1:numel(props)
                pset.(props{i}) = [];
            end
            pset = pset.invalidateDerivedProperties(props);
        end
        
        function varargout = walkPropMeta(pset, applyFn, numargout, varargin)
            % utility for applying a function over all propMeta properties 
            % [out1, out2, ..., outN] = applyFn(dataRoot, propName, propMeta, pathToContainer)
            % where out# will be a struct whos form recapitualates the location of the property within pset
            % pathToRoot will be something like '', 'stored', or 'random' and follow the path down to the prop, e.g. 'stored.structWithinStored'
            
            p = inputParser();
            p.addParameter('propMeta', pset.propMeta, @isstruct);
            p.parse(varargin{:});
            
            propMeta = p.Results.propMeta;
            
            if nargin < 3
                numargout = min(0, nargout(applyFn));
            end
            
            % walk properties in pset
            outputs = walkInner(pset, propMeta, '');
            
            function outputs = walkInner(dataRoot, propMetaRoot, pathToContainer)
                outputs = repmat({struct()}, numargout, 1);
                
                if isempty(pathToContainer)
                    containerStr = 'pset';
                else
                    containerStr = ['pset.' pathToContainer];
                end
 
                props = fieldnames(propMetaRoot);
                for iP = 1:numel(props)
                    prop = props{iP};
                    propMeta = propMetaRoot.(prop);
                    
                    if isstruct(propMeta)
                        % recurse on subfields
                        if isobject(dataRoot)
                            assert(isprop(dataRoot, prop), 'Field %s not found in class %s', prop, class(dataRoot));
                        else
                            assert(isfield(dataRoot, prop), 'Field %s not found in property %s', prop, containerStr);
                        end
                        outputsThis = walkInner(dataRoot.(prop), propMeta, [pathToContainer '.' prop]);
                    else
                        % applyFn to this property
                        outputsThis = cell(numargout, 1);
                        if isobject(dataRoot)
                            assert(isprop(dataRoot, prop), 'Property %s not found in class %s', prop, class(dataRoot));
                        else
                            assert(isfield(dataRoot, prop), 'Field %s not found in property %s', prop, containerStr);
                        end
                        [outputsThis{:}] = applyFn(dataRoot, prop, propMeta, pathToContainer);
                    end
                   
                    % expand the outputs into the output struct
                    for iO = 1:numargout
                        outputs{iO}.(prop) = outputsThis{iO};
                    end
                end
            end
            
            varargout = outputs;
        end
        
        function values = getValuesFromPropMeta(pset, varargin)
            p = inputParser();
            p.addOptional('propMeta', pset.propMeta, @isstruct);
            p.parse(varargin{:});
            
            values = pset.walkPropMeta(@(x, varargin) x, 1, 'propMeta', p.Results.propMeta);
        end
        
        function pset = walkPropMetaAssignValues(pset, newValues, propWasUpdated)
            % replace each propMeta property value in pset with the corresponding value from newValues 
            % if propWasUpdated is provided, it will contain nested booleans in the same form as newValues, and the 
            % value will only be assigned if true.
            
            if nargin < 3 || isempty(propWasUpdated)
                propWasUpdated = true;
            end
            pset = walkInner(pset, pset.propMeta, newValues, propWasUpdated);
            
            function dataRoot = walkInner(dataRoot, propMetaRoot, valueRoot, doAssignRoot)
                props = fieldnames(propMetaRoot);
                for iP = 1:numel(props)
                    prop = props{iP};
                    meta = propMetaRoot.(prop);
                    if ~isfield(valueRoot, prop)
                        continue;
                    end
                    value = valueRoot.(prop);
                    if islogical(doAssignRoot)
                        doAssign = doAssignRoot;
                    elseif isfield(doAssignRoot, prop)
                        doAssign = doAssignRoot.(prop);
                    else
                        doAssign = false;
                    end
                    
                    if isstruct(meta)
                        % recurse on subfields
                        dataRoot.(prop) = walkInner(dataRoot.(prop), meta, value, doAssign);
                    elseif doAssign
                        dataRoot.(prop) = value;
                    end
                end
            end
        end

        function list = convertBooleanWalkStructToPathStrings(pset, walkStruct) %#ok<INUSL>
            % takes a struct like that returned from walkPropMeta whose values are true or false and assembles a list of strings
            
            list = walkInner(walkStruct, '');
            
            function list = walkInner(walkStruct, path)
                props = fieldnames(walkStruct);
                list = cell(0, 1);
                for p = 1:numel(props)
                    val = walkStruct.(props{p});
                    if isempty(path)
                        pathThis = props{p};
                    else
                        pathThis = [path '.' props{p}];
                    end
                    if isstruct(val)
                        list = cat(1, list, walkInner(val, pathThis));
                    elseif islogical(val) && isscalar(val)
                        if val
                            list = cat(1, list, {pathThis});
                        end
                    else
                        error('All fields must contain scalar logical values');
                    end
                end
            end
        end
        
        function [pset, propsInvalidated] = invalidateDerivedProperties(pset, propsChanged, varargin)
            % recursively and efficiently clear out any properties that depend on properties in the set propsChanged
            % handles pset (.propMeta), stored (.storedpropMeta) and .randomized (.randomizedpropMeta)
            % updateFn can be empty, in which case the properties will be cleared to []
            % or it can have signature newValue = updateFn(oldValue, prop, propMeta, pathToContainer, varargin)
            % where propContainer will be 'pset', 'stored', or 'random'
            pset.warnIfNoArgOut(nargout);
            
            p = inputParser();
            p.addParameter('propsAlreadyUpdated', {}, @(x) ischar(x) || iscellstr(x));
            p.addParameter('updateFn', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('errorIfManualNotUpdated', true, @islogical);
            p.parse(varargin{:});

            if ischar(propsChanged)
                propsChanged = {propsChanged};
            end

            % keep a running list of properties that have already been handled
            propsAlreadyUpdated = p.Results.propsAlreadyUpdated;
            if ischar(propsAlreadyUpdated)
                propsAlreadyUpdated = {propsAlreadyUpdated};
            end
            updateFn = p.Results.updateFn;
           
            % apply updateFn or generate cleared values 
            [newData, propsInvalidated] = pset.walkPropMeta(@walkFn, 2);
            
            % assign the new values where invalidated
            pset = pset.walkPropMetaAssignValues(newData, propsInvalidated);
            
            % then recursively invalidate properties dependent on those values just cleared
            listInvalidated = pset.convertBooleanWalkStructToPathStrings(propsInvalidated);
            if ~isempty(listInvalidated)
                pset = pset.invalidateDerivedProperties(listInvalidated);
            end
            
            function [valueOut, thisPropInvalidated] = walkFn(dataRoot, prop, propMeta, pathToContainer)
                if ~isempty(pathToContainer)
                    fullPathProp = [pathToContainer '.' prop];
                    isRootPsetProp = false;
                else
                    fullPathProp = prop;
                    isRootPsetProp = true;
                end
                    
                if ismember(prop, propsChanged) || ismember(prop, propsAlreadyUpdated)
                    thisPropInvalidated = false;
                    valueOut = [];
                else
                    attr = propMeta.attr;

                    if isfield(attr, 'depends') && any(ismember(propsChanged, attr.depends))
                        % prop depends on propChanged, we should invalidate or update it.

                        if isRootPsetProp
                            % this is a property of pset
                            % get the old value without computing on-demand if not already in odc
                            [propValue, useManual, isODC] = pset.processPropGetWithoutComputing(prop);
                        else
                            propValue = dataRoot.(prop);
                            useManual = true;
                        end

                        if useManual || ~isODC
                            if isempty(updateFn)
                                if isRootPsetProp
                                    % internal logic error checking - this needs to be handled directly since we can't just invalidate manual property values
                                    if p.Results.errorIfManualNotUpdated
                                        error('Property %s is stored manually but has been invalidated by the dependency logic', fullPathProp);
                                    else
                                        warning('PropMeta: Property %s is stored manually but has been invalidated by the dependency logic, clearing value', fullPathProp);
                                    end
                                else
                                    if isempty(updateFn)
                                        debug('PropMeta: Clearing %s\n', fullPathProp);
                                    else
                                        debug('PropMeta: Updating %s via provided function\n', fullPathProp);
                                    end
                                end

                                valueOut = [];
                                thisPropInvalidated = true;
                            else
                                if pset.debugPropMeta, debug('PropMeta: Updating property %s via provided function\n', fullPathProp); end
                                valueOut = updateFn(propValue, prop, propMeta, pathToContainer);
                                thisPropInvalidated = true;
                            end
                        else
                            if isempty(propValue)
                                if pset.debugPropMeta, debug('PropMeta: Skipping compute-on-demand %s not yet computed\n', fullPathProp); end
                                valueOut = [];
                                thisPropInvalidated = false;   
                            elseif isempty(updateFn)
                                if pset.debugPropMeta, debug('PropMeta: Invalidating compute-on-demand %s\n', fullPathProp); end
                                valueOut = [];
                                thisPropInvalidated = true;
                            else
                                if pset.debugPropMeta, debug('PropMeta: Updating compute-on-demand %s via provided function\n', prop); end
                                valueOut = updateFn(value, prop, propMeta, pathToContainer);
                                thisPropInvalidated = true;
                            end
                        end % useManual
                        
                    else % this prop does not depend on an invalidated prop?
                        thisPropInvalidated = false;
                        valueOut = [];
                    end
                    
                end % was this one of the props already handled
            end % walkFn
        end

        function pset = transformInternalProperties(pset, fn)
            % fn should look like [data, wasUpdated] = fn(data, propName, propMeta, pathToContainer, varargin)
            % propContainer will be 'pset' or 'stored' or 'random'

            function [dataRoot, wasUpdated] = walkImpl(dataRoot, propName, propMeta, pathToContainer)
                if isa(dataRoot, 'PopulationTrajectorySet')
                    [value, storedManual, isODC] = pset.processPropGetWithoutComputing(propName);
                else
                    value = dataRoot.(propName);
                    storedManual = true;
                    isODC = false;
                end
                if isempty(value) && ~storedManual && isODC
                    if pset.debugPropMeta, debug('PropMeta: transformInternalProperties skipping not-yet-computed %s\n', propName); end
                    wasUpdated = false;
                else
                    [dataRoot, wasUpdated] = fn(value, propName, propMeta, pathToContainer);
                end    
            end
            
            [newValues, wasTransformed] = pset.walkPropMeta(@walkImpl, 2); 
            pset = pset.walkPropMetaAssignValues(newValues, wasTransformed);
        end

        function pset = transformInternalProperties_selectAlongDimension(pset, dimNames, masksByDim)
            pset.warnIfNoArgOut(nargout);

            function [data, wasUpdated] = transformFn(data, prop, propMeta, varargin) %#ok<INUSL>
                [data, wasUpdated] = propMeta.selectAlongDimByName(data, dimNames, masksByDim);
            end

            % apply slicing operation along those dims on all properties
            pset = pset.transformInternalProperties(@transformFn);
        end

        function pset = storeDataInStored(pset, name, data, propMeta, varargin)
            p = inputParser();
            p.addParameter('ignoreOverwrite', false, @islogical);
            p.parse(varargin{:});

            pset.warnIfNoArgOut(nargout);

            if isfield(pset.stored, name) && ~p.Results.ignoreOverwrite
                warning('Overwriting existing stored custom data field %s', name);
            end
            pset.stored.(name) = data;
            pset.storedMeta.(name) = propMeta;
        end

        function varargout = getDataInStored(pset, varargin)
            varargout = cell(numel(varargin), 1);
            for i = 1:numel(varargin)
                name = varargin{i};
                if isfield(pset.stored, name)
                    varargout{i} = pset.stored.(name);
                end
            end
        end

        function pset = clearDataInStored(pset, name)
            pset.warnIfNoArgOut(nargout);

            if isfield(pset.stored, name)
                pset.stored = rmfield(pset.stored, name);
            end
        end

        function tf = hasDataInStored(pset, name)
            tf = isfield(pset.stored, name);
        end
    end

    % get and set properties for ODC properties that utilize propMeta
    methods
        function v = get.basisNames(pset)
            v = pset.processPropGet('basisNames');
        end
        function pset = set.basisNames(pset, v)
            pset = pset.processODCPropSet('basisNames', v);
        end

        function v = get.basisUnits(pset)
            v = pset.processPropGet('basisUnits');
        end
        function pset = set.basisUnits(pset, v)
            pset = pset.processODCPropSet('basisUnits', v);
        end

        function v = get.alignSummaryData(pset)
            v = pset.processPropGet('alignSummaryData');
        end
        function pset = set.alignSummaryData(pset, v)
            pset = pset.processODCPropSet('alignSummaryData', v);
        end

        function v = get.basisAlignSummaryLookup(pset)
            v = pset.processPropGet('basisAlignSummaryLookup');
        end
        function pset = set.basisAlignSummaryLookup(pset, v)
            pset = pset.processODCPropSet('basisAlignSummaryLookup', v);
        end

        function v = get.alignSummaryAggregated(pset)
            v = pset.processPropGet('alignSummaryAggregated');
        end
        function pset = set.alignSummaryAggregated(pset, v)
            pset = pset.processODCPropSet('alignSummaryAggregated', v);
        end

        function v = get.dataByTrialRaw(pset)
            v = pset.processPropGet('dataByTrialRaw');
        end
        function pset = set.dataByTrialRaw(pset, v)
            pset = pset.processODCPropSet('dataByTrialRaw', v);
        end

        function v = get.tMinForDataByTrialRaw(pset)
            v = pset.processPropGet('tMinForDataByTrialRaw');
        end
        function pset = set.tMinForDataByTrialRaw(pset, v)
            pset = pset.processODCPropSet('tMinForDataByTrialRaw', v);
        end

        function v = get.tMaxForDataByTrialRaw(pset)
            v = pset.processPropGet('tMaxForDataByTrialRaw');
        end
        function pset = set.tMaxForDataByTrialRaw(pset, v)
            pset = pset.processODCPropSet('tMaxForDataByTrialRaw', v);
        end

        function v = get.tMinValidByTrialRaw(pset)
            v = pset.processPropGet('tMinValidByTrialRaw');
        end
        function pset = set.tMinValidByTrialRaw(pset, v)
            pset = pset.processODCPropSet('tMinValidByTrialRaw', v);
        end

        function v = get.tMaxValidByTrialRaw(pset)
            v = pset.processPropGet('tMaxValidByTrialRaw');
        end
        function pset = set.tMaxValidByTrialRaw(pset, v)
            pset = pset.processODCPropSet('tMaxValidByTrialRaw', v);
        end
        function v = get.basisValid(pset)
            v = pset.processPropGet('basisValid');
        end
        function pset = set.basisValid(pset, v)
            pset = pset.processODCPropSet('basisValid', v);
        end

        function v = get.basisInvalidCause(pset)
            v = pset.processPropGet('basisInvalidCause');
        end
        function pset = set.basisInvalidCause(pset, v)
            pset = pset.processODCPropSet('basisInvalidCause', v);
        end

        function v = get.trialListsRaw(pset)
            v = pset.processPropGet('trialListsRaw');
        end

        function pset = set.trialListsRaw(pset, v)
            pset = pset.processODCPropSet('trialListsRaw', v);
        end

        function v = get.tMinValidByAlignBasisConditionRaw(pset)
            v = pset.processPropGet('tMinValidByAlignBasisConditionRaw');
        end
        function pset = set.tMinValidByAlignBasisConditionRaw(pset, v)
            pset = pset.processODCPropSet('tMinValidByAlignBasisConditionRaw', v);
        end

        function v = get.tMaxValidByAlignBasisConditionRaw(pset)
            v = pset.processPropGet('tMaxValidByAlignBasisConditionRaw');
        end
        function pset = set.tMaxValidByAlignBasisConditionRaw(pset, v)
            pset = pset.processODCPropSet('tMaxValidByAlignBasisConditionRaw', v);
        end

        function v = get.tMinForDataMean(pset)
            v = pset.processPropGet('tMinForDataMean');
        end
        function pset = set.tMinForDataMean(pset, v)
            pset = pset.processODCPropSet('tMinForDataMean', v);
        end

        function v = get.tMaxForDataMean(pset)
            v = pset.processPropGet('tMaxForDataMean');
        end
        function pset = set.tMaxForDataMean(pset, v)
            pset = pset.processODCPropSet('tMaxForDataMean', v);
        end

        function v = get.dataMean(pset)
            v = pset.processPropGet('dataMean');
        end
        function pset = set.dataMean(pset, v)
            pset = pset.processODCPropSet('dataMean', v);
        end

        function v = get.dataSem(pset)
            v = pset.processPropGet('dataSem');
        end

        function pset = set.dataSem(pset, v)
            pset = pset.processODCPropSet('dataSem', v);
        end

        function v = get.dataNumTrials(pset)
            v = pset.processPropGet('dataNumTrials');
        end

        function pset = set.dataNumTrials(pset, v)
            pset = pset.processODCPropSet('dataNumTrials', v);
        end

        function v = get.dataDifferenceOfTrialsScaledNoiseEstimate(pset)
            v = pset.processPropGet('dataDifferenceOfTrialsScaledNoiseEstimate');
        end

        function pset = set.dataDifferenceOfTrialsScaledNoiseEstimate(pset, v)
            pset = pset.processODCPropSet('dataDifferenceOfTrialsScaledNoiseEstimate', v);
        end
    end

    % Non-dependent get methods that may modify values on retrieval to ensure consistency
    methods
       function sf = get.spikeFilter(pset)
            sf = pset.spikeFilter;
            if isempty(sf)
                sf = SpikeFilter.getDefaultFilter();
            end
       end

       function permValid = get.basisValidPermanent(pset)
            permValid = pset.basisValidPermanent;
            if isempty(permValid)
                permValid = truevec(pset.nBases);
            end
        end

        function permCause = get.basisInvalidCausePermanent(pset)
            permValid = pset.basisValidPermanent;
            permCause = pset.basisInvalidCausePermanent;
            if isempty(permCause)
                permCause = cellstrvec(pset.nBases);
            end
            emptyMask = cellfun(@isempty, permCause);
            permCause(~permValid & emptyMask) = {'marked invalid permanently'};
            permCause(permValid) = {''};
        end

        function tempValid = get.basisValidTemporary(pset)
            tempValid = pset.basisValidTemporary;
            if isempty(tempValid)
                tempValid = truevec(pset.nBases);
            end
        end

        function tempCause = get.basisInvalidCauseTemporary(pset)
            tempValid = pset.basisValidTemporary;
            tempCause = pset.basisInvalidCauseTemporary;
            if isempty(tempCause)
                tempCause = cellstrvec(pset.nBases);
            end
            emptyMask = cellfun(@isempty, tempCause);
            tempCause(~tempValid & emptyMask) = {'marked invalid temporarily'};
            tempCause(tempValid) = {''};
        end
    end

    % Dependent on-the-fly computation (not cached in ODC)
    methods
        function v = get.conditionIncludeMask(pset)
            v = pset.conditionDescriptor.conditionIncludeMask;
        end

        function v = get.dataUnitsCommon(pset)
            units = unique(pset.basisUnits);
            if numel(units) == 1
                v = units{1};
            else
                v = '';
            end
        end

        function d = get.timeDelta(pset)
            d = pset.spikeFilter.timeDelta;
        end

        function v = get.tvecDataMean(pset)
            % generate on the fly, no caching
            v = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                v{iAlign} = makecol(pset.tMinForDataMean(iAlign):pset.timeDelta:pset.tMaxForDataMean(iAlign));
            end
        end

        function v = get.tvecDataMeanForPlotting(pset)
            % generate on the fly, no caching
            v = cellvec(pset.nAlign);
            offsets = pset.getAlignPlottingTimeOffsets(pset.tvecDataMean);
            for iAlign = 1:pset.nAlign
                v{iAlign} = makecol(pset.tvecDataMean{iAlign} + offsets(iAlign));
            end
        end

        function off = get.offsetsTimeDataMeanForPlotting(pset)
            off = pset.getAlignPlottingTimeOffsets(pset.tvecDataMean);
        end

        function v = get.nTimeValidByAlignBasisCondition(pset)
            v = (pset.tMaxValidByAlignBasisConditionRaw - pset.tMinValidByAlignBasisConditionRaw) / pset.timeDelta + 1;
        end

        function v = get.nTimeDataMean(pset)
            v = makecol(cellfun(@numel, pset.tvecDataMean));
        end

        function tf = get.hasDataByTrial(pset)
            % we have dataByTrial if we have data sources, or if the dataByTrial is manually specified
            tf = ~pset.dataInfoManual || pset.dataByTrialManual;
        end

        % has one of the storeDataRandomized* methods been called to
        % populate dataMeanRandomized?
        function tf = get.hasDataRandomized(pset) %#ok<MANU>
            tf = false;
            % TODO
%             tf = ~isempty(pset.randomizedData) && numel(fieldnames(pset.randomizedData)) > 0;
        end

        function c = get.conditionsWithTrialsAllBasesAligns(pset)
            c = squeeze(all(pset.dataNumTrialsRaw(pset.basisValid, :), 1));
        end

        function v = get.tvecDataByTrialRaw(pset)
            % generate on the fly, no caching
            v = cell(pset.nBases,pset.nAlign);
            tMins = pset.tMinForDataByTrialRaw;
            tMaxs = pset.tMaxForDataByTrialRaw;
            delta = pset.timeDelta;
            for iAlign = 1:pset.nAlign
                for iBasis = 1:pset.nBases
                    v{iBasis, iAlign} = makecol(tMins(iBasis,iAlign):delta:tMaxs(iBasis,iAlign));
                end
            end
        end

        function v = get.tvecDataByTrial(pset)
            % generate on the fly, no caching
            v = cell(pset.nBases,pset.nAlign);
            tMins = pset.tMinForDataByTrial;
            tMaxs = pset.tMaxForDataByTrial;
            delta = pset.timeDelta;
            for iAlign = 1:pset.nAlign
                for iBasis = 1:pset.nBases
                    v{iBasis, iAlign} = makecol(tMins(iBasis,iAlign):delta:tMaxs(iBasis,iAlign));
                end
            end
        end

        function v = get.hasValidTrialAverageByBasisCondition(pset)
            v = permute(all(pset.tMinValidByAlignBasisConditionRaw <= pset.tMaxValidByAlignBasisConditionRaw, 1), [2 3 1]);
        end

        function v = get.conditionHasValidTrialAverageAllAlignsBases(pset)
            % here is where conditionIncludeMask is factored in
            hasAvg = pset.hasValidTrialAverageByBasisCondition;pset.manualSliceOrExpandTimeWindow
            v = makecol(squeeze(all(hasAvg(pset.basisValid, :), 1)));
            v(~pset.conditionIncludeMask) = false;
        end

        function v = get.tMinValidAllBasesByAlignCondition(pset)
            % generate on the fly, no caching
            % generate the minimum time over all bases that have a trial
            % average for that align condition

            % @djoshea changing this 20160728 since not all bases need to
            % be present. the window can be the tightest for all bases with
            % any valid window, ignoring the ones that don't have any data
            %  cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            cMask = pset.conditionIncludeMask(:);

            if ~any(cMask)
                v = nan(pset.nAlign, pset.nConditions);
            else
                tMinABC = pset.tMinValidByAlignBasisConditionRaw;
                % for each included condition, take the max tMin over all bases
                tMinABC(:, :, ~cMask) = NaN;
                v = TensorUtils.squeezeDims(nanmax(tMinABC, [], 2),  2);
            end
        end

        function v = get.tMaxValidAllBasesByAlignCondition(pset)
            % generate on the fly, no caching

            % see above comment 20160728
            %             cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            cMask = pset.conditionIncludeMask(:);

            if ~any(cMask)
                v = nan(pset.nAlign, pset.nConditions);
            else
                tMaxABC = pset.tMaxValidByAlignBasisConditionRaw;
                % for each included condition, take the min tMax over all bases
                tMaxABC(:, :, ~cMask) = NaN;
                v = TensorUtils.squeezeDims(nanmin(tMaxABC, [], 2),  2);
            end
        end

        function masks = get.trialHasSpikesMaskByBasis(pset)
            if ~pset.hasDataByTrial
                masks = {};
            else
                masks = cellvec(pset.nBases);
                %                 prog = ProgressBar(pset.nBases, 'Computing has spikes by trial for each basis');
                for iBasis = 1:pset.nBases

                    isNonNanZero = @(x) ~isnan(x) & x~=0;
                    masksByAlign = cellfun(@(x) isNonNanZero(nanmax(x, [], 2)), pset.dataByTrialRaw(iBasis, :), 'UniformOutput', false);
                    masksCat = cat(2, masksByAlign{:});
                    maskAny = any(masksCat, 2);
                    masks{iBasis} = maskAny;
                end
            end
        end
    end

    % internal build methods that write directly to ODC
    methods(Hidden)
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
            c.data.basisNames = basisNames;
            c.data.basisUnits = basisUnits;
        end

        function buildBasisValid(pset)
            permInvalid = ~pset.basisValidPermanent;
            permCause = pset.basisInvalidCausePermanent;
            tempInvalid = ~pset.basisValidTemporary & ~permInvalid; % for convenience in assigning causes below
            tempCause = pset.basisInvalidCauseTemporary;

            % also mark additional bases temporarily invalid that have no
            % valid trial averages for any condition or align
            cMask = pset.conditionIncludeMask;
            missingTrialAvg = ~permInvalid & ~tempInvalid & ~makecol(any(pset.hasValidTrialAverageByBasisCondition(:, cMask), 2));

            valid = ~permInvalid & ~tempInvalid & ~missingTrialAvg;

            % use permanent cause first
            cause = cellstrvec(pset.nBases);
            cause(permInvalid) = cellfun(@(x) ['Permanent: ' x], permCause(permInvalid), 'UniformOutput', false);

            % then temporary cause if not permanently invalid
            cause(tempInvalid) = cellfun(@(x) ['Temporary: ' x], tempCause(tempInvalid), 'UniformOutput', false);
            cause(missingTrialAvg) = {'Missing: basis has no valid trial averages for any align./ condition'};

            c = pset.odc; % this notation makes explicit that we're manipulating the internal handle
            c.data.basisValid = valid;
            c.data.basisInvalidCause = cause;
        end

        function buildDataByTrial(pset)
            % stores dataByTrial, tMinByTrial, and tMaxByTrial in odc
            % This method fetches the aligned data for EVERY trial in
            % each data source, for each alignment. It does NOT consider
            % the condition grouping, allowing this to be handled later.
            % It stores this data as a matrix, whose time vector is given
            % by the widest time vector along any trial. Missing samples in
            % this matrix are NaNs.

            if pset.dataInfoManual
                return;
            end

            dataByTrial = cell(pset.nBases, pset.nAlign);

            [tMinForDataByTrial, tMaxForDataByTrial] = deal(nan(pset.nBases, pset.nAlign));

            % alignSummary instances are built by dataSource, so the
            % basis to alignSummary lookup is the same as the basis to dataSource lookup
            %             basisAlignSummaryLookup = pset.basisDataSourceIdx;
            %             alignSummaryData = cell(pset.nDataSources, pset.nAlign);

            dataSourcesByBasis = pset.dataSources(pset.basisDataSourceIdx);
            basisDataSourceChannelNames = pset.basisDataSourceChannelNames; %#ok<*PROP>
            nAlign = pset.nAlign;
            nBases = pset.nBases;
            spikeFilter = pset.spikeFilter;

            isSpikeChannel = false(pset.nBases, 1);

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
                        % use the specs from the spikeFilter
                        [dataByTrial{iBasis, iAlign}, tvec] = src.getAnalogAsMatrix(chName, ...
                            'timeDelta', spikeFilter.timeDelta, 'resampleMethod', spikeFilter.resampleMethod, ...
                            'binAlignmentMode', spikeFilter.binAlignmentMode);
                        isSpikeChannel(iBasis) = false;

                    elseif src.hasSpikeChannel(chName)
                        [dataByTrial{iBasis, iAlign}, tvec] = src.getSpikeRateFilteredAsMatrix(chName, ...
                            'spikeFilter', spikeFilter);
                        isSpikeChannel(iBasis) = true;
                    else
                        error('Unknown channel type');
                    end

                    % store the time limits used in the time vector for
                    % the dataByTrial{..} matrix
                    if ~isempty(tvec)
                        tMinForDataByTrial(iBasis, iAlign) = min(tvec);
                        tMaxForDataByTrial(iBasis, iAlign) = max(tvec);
                    end
                end
            end
            prog.finish();

            % apply translation / normalization to data
            if ~isempty(pset.translationNormalization)
                dataByTrial = pset.translationNormalization.applyTranslationNormalizationToData(dataByTrial);
            end

            % store the results in the odc without copying
            c = pset.odc;
            c.data.dataByTrialRaw = dataByTrial;
            c.data.tMinForDataByTrialRaw = tMinForDataByTrial;
            c.data.tMaxForDataByTrialRaw = tMaxForDataByTrial;
        end

        function buildDataByTrialPerTrialLimits(pset)
            % stores tMinByTrial and tMaxByTrial in odc
            % Store the precise time starts and stops for EACH
            % trial that comprises that matrix.

            if ~pset.hasDataByTrial
                return;
            end

            [tMinByTrial, tMaxByTrial] = ...
                deal(cell(pset.nBases, pset.nAlign));

            nAlign = pset.nAlign;
            nBases = pset.nBases;
            dataByTrial = pset.dataByTrialRaw;
            tvecCell = pset.tvecDataByTrialRaw;
            for iBasis = 1:nBases
                for iAlign = 1:nAlign
                    [tMinByTrial{iBasis, iAlign}, ...
                        tMaxByTrial{iBasis, iAlign}] = ...
                        TrialDataUtilities.Data.getValidTimeExtents(...
                            tvecCell{iBasis, iAlign}, dataByTrial{iBasis, iAlign});
                end
            end

            % store the results in the odc without copying
            c = pset.odc;
            c.data.tMinValidByTrialRaw = tMinByTrial;
            c.data.tMaxValidByTrialRaw = tMaxByTrial;
        end

        function buildTrialLists(pset)
            if pset.dataInfoManual
                return;
            end
            % computes and stores dataNumTrialsRaw and dataValid into odc
            trialLists = cell(pset.nBases, pset.nConditions);

            hasSpikesByBasis = pset.trialHasSpikesMaskByBasis;
            prog = ProgressBar(pset.nBases, 'Computing trial-counts by condition');
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                % note, this src will not be aligned to this iAlign,
                % but this isn't necessary since we've already
                % extracted the aligned data
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};

                listsThis = src.conditionInfo.listByCondition(:);

                % filter out trials that don't have spikes, or lie outside
                % the band of trials with spikes
                if pset.ignoreAllZeroSpikeTrials || pset.ignoreLeadingTrailingZeroSpikeTrials
                    if pset.ignoreAllZeroSpikeTrials
                        keepTrials = find(hasSpikesByBasis{iBasis});
                    else
                        mask = hasSpikesByBasis{iBasis};
                        firstNonZero = find(mask, 1, 'first');
                        lastNonZero = find(mask, 1, 'last');
                        keepTrials = firstNonZero:lastNonZero;
                    end

                    for iC = 1:numel(listsThis)
                        listsThis{iC} = intersect(listsThis{iC}, keepTrials);
                    end
                end

                trialLists(iBasis, :) = listsThis;
            end
            prog.finish();

            c = pset.odc;
            c.data.trialListsRaw = trialLists;
        end

        function buildTimeWindowsByAlignBasisCondition(pset)
            % computes and stores tMin/MaxValidByBasisAlignCondition, the
            % time windows for each basis over which enough trials exist to
            if pset.dataInfoManual
                return;
            end

            % first, compute the all-inclusive time window for each basis,
            % for each alignment, using only condition and align valid trials
            [tMinValidByAlignBasisConditionRaw, tMaxValidByAlignBasisConditionRaw] = ...
                deal(nan(pset.nAlign, pset.nBases, pset.nConditions));

            % do this first to force computation of data by trial at the beginning,
            % rather than having it happen on the first loop iteration
            temp = pset.tMinValidByTrialRaw; %#ok<NASGU>
            dataNumTrialsRaw = cellfun(@numel, pset.trialListsRaw);

            prog = ProgressBar(pset.nBases, 'Computing trial-averaged time windows by basis/align/condition');

            for iBasis = 1:pset.nBases
                for iAlign = 1:pset.nAlign
                    prog.update(iBasis);

                    % note, this src will not be aligned to this iAlign,
                    % but this isn't necessary since we've already
                    % extracted the aligned data
                    src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                    tMinByTrial = pset.tMinValidByTrialRaw{iBasis, iAlign};
                    tMaxByTrial = pset.tMaxValidByTrialRaw{iBasis, iAlign};

                    % group the condition windows by conditionLists
                    [tMinByTrialGrouped, tMaxByTrialGrouped] = src.conditionInfo.groupElements(tMinByTrial, tMaxByTrial);

                    % for each basis, align, take the widest window we can that
                    % is valid for a sufficient number of trials on this basis
                    for iCondition = 1:pset.nConditions
                        nTrialsThis = dataNumTrialsRaw(iBasis, iCondition);
                        if nTrialsThis == 0, continue, end
                        trialCountThresh = max(pset.minTrialsForTrialAveraging, ...
                            ceil(pset.minFractionTrialsForTrialAveraging*nTrialsThis));

                        % it's possible for a trial to be valid, but for
                        % there to be no valid data for this trial because
                        % the particular channel being requested had no
                        % data here. So we decide here that
                        % minFractionTrialsForTrialAveraging is relative to
                        % the total number of trials included, which might
                        % mean there is no data for this window, but that's
                        % what the user has requested
                        nTrialsThisWithData = nnz(~isnan(tMinByTrialGrouped{iCondition}) & ~isnan(tMaxByTrialGrouped{iCondition}));
                        if nTrialsThisWithData == 0, continue, end
                        if nTrialsThisWithData >= trialCountThresh
                            tMinSorted = sort(removenan(tMinByTrialGrouped{iCondition}), 1, 'ascend');
                            tMinValidByAlignBasisConditionRaw(iAlign, iBasis, iCondition) = tMinSorted(trialCountThresh);
                            tMaxSorted = sort(removenan(tMaxByTrialGrouped{iCondition}), 1, 'descend');
                            tMaxValidByAlignBasisConditionRaw(iAlign, iBasis, iCondition) = tMaxSorted(trialCountThresh);
                        end
                    end
                end
            end
            prog.finish();

            c = pset.odc;
            c.data.tMinValidByAlignBasisConditionRaw = tMinValidByAlignBasisConditionRaw;
            c.data.tMaxValidByAlignBasisConditionRaw = tMaxValidByAlignBasisConditionRaw;
        end

        function buildTvecDataMean(pset)
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
            %             conditionsWithTrialsAllBasesAligns = squeeze(all(all(pset.dataNumTrialsRaw(:, basisMask, :), 1), 2));
            %             if ~any(conditionsWithTrialsAllBasesAligns)
            %                 error('No conditions have trial averages for all bases on all aligns. Try lowering minTrialsForTrialAveraging or minFractionTrialsForTrialAveraging?');
            %             end
            %
            %             tMinValidByAlignCondition = TensorUtils.squeezeDims(max(pset.tMinValidByAlignBasisConditionRaw(:, basisMask, conditionsWithTrialsAllBasesAligns), [], 2), 2);
            %             tMaxValidByAlignCondition = TensorUtils.squeezeDims(min(pset.tMaxValidByAlignBasisConditionRaw(:, basisMask, conditionsWithTrialsAllBasesAligns), [], 2), 2);

            % for each align, compute the widest window valid for ANY
            % considered condition for ALL valid bases. A condition is
            % considered if at least one basis has a valid trial average
            % for it on this align.

            % 20171221 - moving to its own function, then
            % changing these time windows to reflect the
            % condition include mask, since the others dont matter

            cMask = pset.conditionIncludeMask;
            tMinForDataMean = makecol(nanmin(pset.tMinValidAllBasesByAlignCondition(:, cMask), [], 2));
            tMaxForDataMean = makecol(nanmax(pset.tMaxValidAllBasesByAlignCondition(:, cMask), [], 2));

            alignInvalid = isnan(tMinForDataMean) | isnan(tMaxForDataMean);
            if any(alignInvalid)
                warning('No valid time window is valid across all bases across all conditions for alignment %s. Use .explain to advise', ...
                    TrialDataUtilities.String.strjoin(find(alignInvalid), ',')); %#ok<FNDSB>
            end

            for iAlign = 1:pset.nAlign
                if alignInvalid(iAlign)
                    continue;
                end
                % realign time vector to run through t=0 with increments of timeDelta
                tempVec = TrialDataUtilities.Data.linspaceIntercept(tMinForDataMean(iAlign), pset.timeDelta, tMaxForDataMean(iAlign), 0);
                tMinForDataMean(iAlign) = min(tempVec);
                tMaxForDataMean(iAlign) = max(tempVec);
            end

            % store in odc without copying
            c = pset.odc;
            c.data.tMinForDataMean = tMinForDataMean;
            c.data.tMaxForDataMean = tMaxForDataMean;
        end

        function buildDataMean(pset)
            % computes and stores dataMean, dataIntervalHigh/Low, and dataNumTrialsRaw into odc
            % this function computes summary statistics across trials,
            % especially dataMean. It selects time windows which are
            % consistent across all bases and conditions for a specific
            % align to simplify subsequent operations.

            % warn about valid bases missing trials in at least one condition for
            % which other bases have trials. This is because the existence
            % of these bases will invalidate ALL data on those conditions
            % since there is no time window with valid trial averaged data
            % across ALL bases.
            %
            % update: we do not warn anymore since the concept of temporary invalid
            %   already covers these bases.
            % pset.warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns();

            tMinForDataMean = pset.tMinForDataMean;
            tMaxForDataMean = pset.tMaxForDataMean;
            basisMask = pset.basisValid;
            nTimeByAlign = arrayfun(@(tmin, tmax) numel(tmin:pset.timeDelta:tmax), tMinForDataMean, tMaxForDataMean);

            % now that we've determined the time window, we can compute the
            % trial average using data from these windows
            [dataMean, dataSem] = deal(cellvec(pset.nAlign));
            for iAlign = 1:pset.nAlign
                [dataMean{iAlign}, dataSem{iAlign}] = ...
                    deal(nan(pset.nBases, pset.nConditions, nTimeByAlign(iAlign)));
            end

            if pset.nBasesValid == 0
                warning('No valid bases found to compute trial-averaged data. Check .basisValid and .basisInvalidCause');
            end

            % copy temp values for for slicing
            %             conditionHasValidTrialAverageAllAlignsBases = pset.conditionHasValidTrialAverageAllAlignsBases;
            dataByTrial = pset.dataByTrialRaw;
            tMinForDataByTrial = pset.tMinForDataByTrialRaw;
            tMaxForDataByTrial = pset.tMaxForDataByTrialRaw;
            timeDelta = pset.timeDelta;
            minTrialsForTrialAveraging = pset.minTrialsForTrialAveraging;
            minFractionTrialsForTrialAveraging = pset.minFractionTrialsForTrialAveraging;

            trialLists = pset.trialListsRaw;
            cMask = pset.conditionIncludeMask;

            for iAlign = 1:pset.nAlign
                prog = ProgressBar(pset.nBases, 'Computing trial-averaged data for align %d', iAlign);

                basisHasNoValidTimepoints = falsevec(pset.nBases);
                for iBasis = 1:pset.nBases
                    prog.update(iBasis);
                    % don't process invalid bases, leave these as NaNs
                    if ~basisMask(iBasis)
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
                        basisHasNoValidTimepoints(iBasis) = true;
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
                    % have all shorter trials

                    tMaskValid = tvecAll >= tMinValid & tvecAll <= tMaxValid;

                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    byTrialValid = byTrial(:, tMaskValid);

                    % IMPORTANT no need to remove zero spike trials since they are
                    % omitted in trialLists already

                    byCondition = cellfun(@(idx) byTrialValid(idx,:), trialLists(iBasis, :)', ...
                        'UniformOutput', false);

                    for iCondition = 1:pset.nConditions
                        %                         if ~conditionHasValidTrialAverageAllAlignsBases(iCondition), continue, end
                        if ~cMask(iCondition), continue; end
                        mat = byCondition{iCondition};
                        nTrials = size(mat, 1);
                        if nTrials == 0, continue, end

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

            %             if ~any(pset.conditionHasValidTrialAverageAllAlignsBases)
            %                 warning('No conditions have valid trial averages for all alignments and bases, so no dataMean will be computed. Use .findBasesMissingTrialAverageForNonEmptyConditionAligns to identify the bases lacking trial averages on each condition');
            %             end

            if any(basisHasNoValidTimepoints)
                warning('Valid trial averages were not computed (left NaN) for %d bases. Check .explain()', nnz(basisHasNoValidTimepoints));
            end

            % no need to apply translation / normalization here since it is
            % already applied to dataByTrial!

            % store in odc without copying
            c = pset.odc;
            c.data.dataMean = dataMean;
            c.data.dataSem = dataSem;
        end

%         function buildDataRandomizedIntervals(pset)
%             % compute the quantiles to use as intervals
%             [dataIntervalHigh, dataIntervalLow] = deal(cellvec(pset.nAlign));
%             qLow = pset.dataIntervalQuantileLow;
%             qHigh = pset.dataIntervalQuantileHigh;
%             for iAlign = 1:pset.nAlign
%                 dataIntervals = quantile(pset.dataMeanRandomized{iAlign}, [qLow qHigh], 4);
%                 dataIntervalLow{iAlign} = dataIntervals(:, :, :, 1);
%                 dataIntervalHigh{iAlign} = dataIntervals(:, :, :, 2);
%             end
%
%             if ~isempty(pset.translationNormalization)
%                 trNorm = pset.translationNormalization;
%                 dataIntervalHigh = cellfun(@trNorm.applyTranslationNormalizationToData, dataIntervalHigh, 'UniformOutput', false);
%                 dataIntervalLow = cellfun(@trNorm.applyTranslationNormalizationToData, dataIntervalLow, 'UniformOutput', false);
%             end
%
%             c = pset.odc;
%             c.data.dataIntervalHigh = dataIntervalHigh;
%             c.data.dataIntervalLow = dataIntervalLow;
%         end

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
            c.data.basisAlignSummaryLookup = basisAlignSummaryLookup;
            c.data.alignSummaryData = alignSummaryData;
        end

        function buildAlignSummaryAggregated(pset)
            % build the aggregated data across all bases too, for each alignment
            alignSummaryAggregated = cell(pset.nAlign, 1);
            alignSummaryData = pset.alignSummaryData;
            dsMask = unique(pset.basisAlignSummaryLookup(pset.basisValid));
            prog = ProgressBar(pset.nAlign, 'Computing aggregate alignment summary statistics');
            for iAlign = 1:pset.nAlign
                prog.update(iAlign);
                if ~any(dsMask)
                    alignSummaryAggregated{iAlign} = AlignSummary.buildEmptyFromConditionAlignDescriptor(pset.conditionDescriptor, pset.alignDescriptorSet{iAlign});
                else
                    alignSummaryAggregated{iAlign} = AlignSummary.buildByAggregation(alignSummaryData(dsMask, iAlign));
                end
            end
            prog.finish();

            c = pset.odc;
            c.data.alignSummaryAggregated = alignSummaryAggregated;
        end

        function buildDataNumTrials(pset)
            numTrials = cellfun(@numel, pset.trialListsRaw);
            numTrials(~pset.basisValid, :) = NaN;
            numTrials(:, ~pset.conditionIncludeMask) = NaN;

            valid = numTrials > 0;

            c = pset.odc;
            c.data.dataNumTrials = numTrials;
            c.data.dataMeanValid = valid;
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

            c = pset.odc;
            c.dataDifferenceOfTrialsScaledNoiseEstimate = computeDataNoiseEstimate(pset);

            function eachAlign_diff_NbyTbyC = computeDataNoiseEstimate(pset)
                if ~pset.hasDataByTrial
                    % TODO can implement based on cached trials here if needed
                    eachAlign_diff_NbyTbyC = [];
                    return;
                    %                 error('Cannot store noise estimate without single trial data. You can pre-compute this using storeDataNoiseEstimate while single trial data is still present');
                end

                % grab random pair of trials for each basis, each condition,
                % concatenate aligns in time
                debug( 'Building difference of trials noise estimate\n');
                eachAlign_NbyTbyCby2 = pset.arrangeEachAlign_NbyTbyCbyR('numTrials', 2, 'ignoreTrialsWithTooFewSamples', true, ...
                    'samplingMode', 'randomWithReplacement');

                % grab trial counts
                nTrials_NbyC = pset.dataNumTrials;

                % diff be N x T x C
                eachAlign_diff_NbyTAbyC = cellfun(@(t) diff(t, 1, 4), eachAlign_NbyTbyCby2, 'UniformOutput', false);

                % scale appropriately each timeseries (running along dim 2) taking into account trial counts
                nTrials_Nby1byC = permute(nTrials_NbyC, [1 3 2]);

                % N x TA x C
                eachAlign_diff_NbyTbyC = cellfun(@(tensor) bsxfun(@rdivide, tensor, sqrt(2*nTrials_Nby1byC)), eachAlign_diff_NbyTAbyC, 'UniformOutput', false);

                if ~isempty(pset.translationNormalization)
                    eachAlign_diff_NbyTbyC = cellfun(@(tensor) pset.translationNormalization.applyNormalizationToData(tensor), eachAlign_diff_NbyTAbyC, 'UniformOutput', false);
                end
            end
        end

%         function eachAlign_diff_NbyTbyC = getDataNoiseEstimate(pset)
%             % get if pre-computed or compute and store if not
%             if ~pset.hasDataInStored('dataDifferenceOfTrialsScaledNoiseEstimate')
%                 eachAlign_diff_NbyTbyC = pset.computeDataNoiseEstimate();
%             else
%                 eachAlign_diff_NbyTbyC = pset.getDataInStored('dataDifferenceOfTrialsScaledNoiseEstimate');
%             end
%         end

        %% DATA BY TRIAL COMMON TIME
        function dataByTrial = computeDataByTrialCommonTime(pset)
            % dataByTrial will nBases x nAlign cell containing ordered data by trial
            % each cell contains nTrials x nTime analog data for that
            % basis, where nTime is now set by .tvecDataMean instead of
            % .tvecDataByTrial

            assert(pset.hasDataByTrial, 'PopulationTrajectorySet must have dataByTrial');

            nBases = pset.nBases;
            tMinDataMean = pset.tMinForDataMean;
            tMaxDataMean = pset.tMaxForDataMean;
            tvecDataMean = pset.tvecDataMean;

            dataByTrial = cell(nBases, pset.nAlign);
            for iAlign = 1:pset.nAlign
                if nBases > 1
                    prog = ProgressBar(nBases, 'Computing dataByTrial for align %d', iAlign);
                end
                for iBasis = 1:nBases
                    if nBases > 1
                        prog.update(iBasis);
                    end

                    % don't process invalid bases, leave these as NaNs
                    if ~pset.basisValid(iBasis)
                        continue;
                    end

                    % lookup the time limits which describe the byTrial
                    % matrix
                    tMinAll = pset.tMinForDataByTrial(iBasis, iAlign);
                    tMaxAll = pset.tMaxForDataByTrial(iBasis, iAlign);
                    if isnan(tMinAll) || isnan(tMaxAll)
                        error('Basis %d has no valid timepoints', iBasis);
                    end
                    tvecAll = tMinAll:pset.timeDelta:tMaxAll;

                    % lookup the new time limits which we'll compute the
                    % trial average within
                    tMinValid = tMinDataMean(iAlign);
                    tMaxValid = tMaxDataMean(iAlign);
                    tMaskValid = tvecAll >= tMinValid & tvecAll <= tMaxValid;
                    takeFirst = find(tMaskValid, 1);
                    insertFirst = find(tvecDataMean{iAlign} == tvecAll(takeFirst), 1);

                    tMaskInsert = insertFirst + (0:nnz(tMaskValid)-1);

                    % grab the valid time portion of the nTrials x
                    % nTime data matrix
                    thisMat = pset.dataByTrial{iBasis, iAlign};
                    dataByTrial{iBasis, iAlign} = nan(size(thisMat, 1), numel(tvecDataMean{iAlign}));
                    dataByTrial{iBasis, iAlign}(:, tMaskInsert) = pset.dataByTrial{iBasis, iAlign}(:, tMaskValid);
                end
                prog.finish();
            end
        end

        function dataByTrialGrouped = computeDataByTrialCommonTimeGrouped(pset)
            dataByTrial = pset.computeDataByTrialCommonTime();
            trialLists = pset.trialLists;
            prog = ProgressBar(pset.nBases, 'Grouping dataByTrial into conditions');

            dataByTrialGrouped = cell(pset.nBases, pset.nAlign, pset.nConditions);
            nTrials = zeros(pset.nBases, pset.nAlign, pset.nConditions);
            for iBasis = 1:pset.nBases
                prog.update(iBasis);
                if ~pset.basisValid(iBasis), continue; end
                for iAlign= 1:size(dataByTrial, 2)
                    byCondition = cellfun(@(idx) dataByTrial{iBasis, iAlign}(idx,:), trialLists(iBasis, :)', ...
                        'UniformOutput', false);

                    dataByTrialGrouped(iBasis, iAlign, :) = byCondition;
                    nTrials(iBasis, iAlign, :) = cellfun(@(x) size(x, 1), byCondition);
                end
            end
            prog.finish();
        end

        function dataByTrialGrouped = getDataByTrialCommonTimeGrouped(pset)
            % get if pre-computed or compute and store if not
            if ~pset.hasDataInStored('dataByTrialCommonTimeGrouped')
                dataByTrialGrouped = pset.computeDataByTrialCommonTimeGrouped();
            else
                dataByTrialGrouped = pset.getDataInStored('dataByTrialCommonTimeGrouped');
            end
        end

        function pset = storeDataByTrialCommonTimeGrouped(pset)
            % inserts into .stored
            % dataByTrialGrouped will nBases x nAlign x nConditions cell containing ordered data by trial
            % each cell contains nTrialsThisCondition x nTime analog data for that
            % basis, where nTime is now set by .tvecDataMean instead of
            % .tvecDataByTrial.
            %
            % nTrials will be nBases x nAlign x nConditions

            % here we do the work of generating the dataByTrial for all
            % trials first, since non-grouped all trial version of this function may be a useful function in the
            % future to have. Then we group it here and store it in ODC as
            % a property worth storing.
            pset.warnIfNoArgOut(nargout);

            dataByTrialGrouped = pset.computeDataByTrialCommonTimeGrouped

            propMeta = PropertyShapeMeta({{'N', 'A'}, {'R', 'T'}}, 'odc', true, 'translate', true, 'normalize', true);
            pset = pset.storeDataInStored('dataByTrialCommonTimeGrouped', dataByTrialGrouped, propMeta);
        end

        function [sampledTrials, sampledTrialCounts] = computeSampledTrials(pset, varargin)
            % We suppress additional arguments like validBasesOnly
            % because we want the size and shape of this to be reliable so that it can be sized
            p = inputParser();
            p.addRequired('numTrials', @isscalar);
            p.KeepUnmatched = false;
            p.parse(varargin{:});

            [sampledTrials, indexInfo] = pset.arrangeEachAlign_RbyTbyNbyC('sampleMethod', 'randomWithReplacement', ...
                'numTrials', p.Results.numTrials, p.Unmatched);
            sampledTrialCounts = indexInfo.numTrialsActual; % N x C
        end

        % store cached single trials in dataSampledTrials as
        function pset = storeSampledTrials(pset, varargin)
            % takes same parameters as arrangeEachAlign_RbyTbyNbyC, in particular,
            % you must provide 'numTrials' argument. We suppress additional arguments like validBasesOnly
            % because we want the size and shape of this to be reliable
            p = inputParser();
            p.addRequired('numTrials', @isscalar);
            p.addParameter('prefix', 'sampledTrials', @ischar);
            p.KeepUnmatched = false;
            p.parse(varargin{:});

            pset.warnIfNoArgOut(nargout);

            [sampledTrials, sampledTrialCounts] = pset.computeSampledTrials('numTrials', p.Results.numTrials);

            propMeta = PropertyShapeMeta({{'A'}, {'R', 'T', 'N', 'C'}}, 'translate', true, 'normalize', true);
            pset = pset.storeDataInStored(sprintf('%s_data', p.Results.prefix), sampledTrials, propMeta, 'ignoreOverwrite', true);

            propMeta = PropertyShapeMeta({{'N', 'C'}}, 'translate', true, 'normalize', true);
            pset = pset.storeDataInStored(sprintf('%s_trialCounts', p.Results.prefix), sampledTrialCounts, propMeta, 'ignoreOverwrite', true);
        end

        function tf  = hasSampledTrials(pset, varargin)
            p = inputParser();
            p.addParameter('prefix', 'sampledTrials', @ischar);
            p.parse(varargin{:});
            sampledTrialsField = sprintf('%s_data', p.Results.prefix);
            sampledTrialCountsField = sprintf('%s_trialCounts', p.Results.prefix);

            tf = pset.hasDataInStored(sampledTrialsField) && ~pset.hasDataInStored(sampledTrialCountsField);
        end

        function [sampledTrials, sampledTrialCounts] = getSampledTrials(pset, varargin)
            % fetches the sampled trials from stored
            p = inputParser();
            p.addParameter('prefix', 'sampledTrials', @ischar);
            p.parse(varargin{:});

            sampledTrialsField = sprintf('%s_data', p.Results.prefix);
            sampledTrialCountsField = sprintf('%s_trialCounts', p.Results.prefix);
            if ~pset.hasDataInStored(sampledTrialsField) || ~pset.hasDataInStored(sampledTrialCountsField)
                error('%s not found. Call pset = pset.storeSampledTrials() to store it', p.Results.prefix);
            end

            [sampledTrials, sampledTrialCounts] = pset.getDataInStored(sampledTrialsField, sampledTrialCountsField);
        end

        function [meansExcludingTrials_A_NbyCbyTbyR, sampledTrials_A_NbyCbyTbyR, nTrials_NbyC] = computeDataMeanExcludingSampledTrials(pset, varargin)
            % Used primarily for cross-validation with DPCA. This function can use one of two sources for data.
            % If the pset has cached trials already, we'll use those to generate the means excluding those trials.
            % If not, we'll use hasDataByTrial to sample trials directly. If ~hasDataByTrial and no cached single trials are present,
            % an error is raised.
            %
            % We use a simplifying
            % assumption here that every trial contributes to the mean,
            % even though this is not necessarily true. Technically the
            % number of trials contributing to the mean varies over time
            % But we don't keep track of this information anywhere. So
            % instead we use dataNumTrialsRaw as the number of trials
            % contributing to the mean everywhere when computing the
            % average without this trial

            p = inputParser();
            p.addParameter('numTrialsDefault', 10, @isscalar); % used if there are no sampled trials already
            p.addParameter('prefix', 'sampledTrials', @ischar);
            p.addParameter('validBasesOnly', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if pset.hasSampledTrials('prefix', p.Results.prefix)
                % use the sampled trials
                [sampledTrials, sampledTrialCounts] = pset.getSampledTrials(pset, 'prefix', p.Results.prefix);
            else
                % sample trials now but dont bother storing them
                [sampledTrials, sampledTrialCounts] = pset.computeSampledTrials(pset, 'numTrials', p.Results.numTrialsDefault);
            end

            % sampledTrials is A x 1 of R x T x N x C x numTrials
            % sampledTrialCounts is N x C

            % dataMean is A x 1 of N x C x T
            % output is going to be A x 1 of N x C x T x numTrials

            [meansExcludingTrials_A_NbyCbyTbyR, sampledTrials_A_NbyCbyTbyR] = deal(cell(pset.nAlign));

            % N x C --> N x 1 x C
            nTrials_NbyC = sampledTrialCounts;
            nTrials_N1C = permute(sampledTrialCounts, [1 3 2]);

            dataMean = pset.dataMean;

            for iA = 1:pset.nAlign
                % R x T x N x C --> N x C x T x R
                sampledTrials_A_NbyCbyTbyR{iA} = permute(sampledTrials{iA}, [3 4 2 1]);

                % N x C x T
                dataSum = bsxfun(@times, dataMean{iA}, nTrials_N1C);

                meansExcludingTrials_A_NbyCbyTbyR{iA} = bsxfun(@rdivide, bsxfun(@minus, dataSum, sampledTrials_A_NbyCbyTbyR{iA}), nTrials_N1C-1);

                if p.Results.validBasesOnly
                    meansExcludingTrials_A_NbyCbyTbyR{iA} = meansExcludingTrials_A_NbyCbyTbyR{iA}(pset.basisValid, :, :, :, :, :);
                    sampledTrials_A_NbyCbyTbyR{iA} = sampledTrials_A_NbyCbyTbyR{iA}(pset.basisValid, :, :, :, :, :);
                end
            end
            if p.Results.validBasesOnly
                nTrials_NbyC = nTrials_NbyC(pset.basisValid, :);
            end
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
        function pset = precomputeProperties(pset)
            pset.dataMean;
            pset.alignSummaryData;
            pset.alignSummaryAggregated;
            pset.trialListsRaw;
            pset.dataDifferenceOfTrialsScaledNoiseEstimate;
        end

        function saveFast(pset, location, varargin)
            p = inputParser();
            p.addParameter('recursive', false, @islogical); % calls saveFast on each source too
            p.addParameter('keepComputed', false, @islogical); % drop everything that can later be computed to save space
            p.parse(varargin{:});

            sources = pset.dataSources;

            if ~pset.dataInfoManual
                pset.dataSources = 'saved separately, use PopulationTrajectorySet.loadFast to load';
            end

            if ~p.Results.keepComputed && ~pset.keepComputedOnSaveFast
                % comment this out to not clear this to save the alignment we've done
                pset.odc = [];
            end

            if ~exist('location', 'var') || isempty(location)
                location = pwd;
            end

            mkdirRecursive(location);
            savefast(fullfile(location, 'pset.mat'), 'pset');

            % save elements of sources
            if ~pset.dataInfoManual
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

            if ~exist('location', 'var') || isempty(location)
                location = pwd;
            end
            if ~exist(location, 'dir')
                error('Directory %s not found. Did you save with saveFast?', location);
            end
            loaded = load(fullfile(location, 'pset.mat'));
            pset = loaded.pset;

            % load elements of sources
            if ~pset.dataInfoManual
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
            tf = ~pset.dataInfoManual;
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
                pset.odc = OnDemandCache();
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

            if isempty(pset.dataInfoManual)
                pset.dataInfoManual = false;
            end

            if isempty(pset.dataByTrialManual)
                pset.dataByTrialManual = pset.dataInfoManual;
            end

            if isempty(pset.dataSources)
                pset.dataSources = {};
            end

            if isempty(pset.basisDataSourceChannelNames)
                pset.basisDataSourceChannelNames = {};
            end

            if isempty(pset.alignDescriptorSet)
                pset = pset.setAlignDescriptorSet({AlignDescriptor()});
            end

            if isempty(pset.conditionDescriptor)
                pset = pset.setConditionDescriptor(ConditionDescriptor());
            end
        end
    end

    % Display / description
    methods
        function printDescription(pset)
            if pset.dataInfoManual
                dataSourceStr = '';
            else
                dataSourceStr = sprintf(', %d data sources', pset.nDataSources);
            end

            % get.basisValid requires data to be extracted which is time
            % consuming, here we collect the data
            if isfield(pset.odc.data, 'basisValid')
                valid = pset.odc.data.basisValid;
            else
                valid = [];
            end
            if isempty(valid)
                tcprintf('inline', '{yellow}%s: {bright white}%d bases {red}(??? valid){none}, %d permanently invalid, {bright white}%d conditions, %d alignments%s\n', ...
                    class(pset), pset.nBases, pset.nBasesPermanentlyInvalid, pset.nConditions, pset.nAlign, dataSourceStr);
                tcprintf('inline', '  {darkGray}Note: basis validity will be determined after .trialLists is computed.\n');
            else
                tcprintf('inline', '{yellow}%s: {bright white}%d bases ({red}%d valid:{none} %d perm, %d temp invalid), {bright white}%d conditions, %d alignments%s\n', ...
                    class(pset), pset.nBases, pset.nBasesValid, pset.nBasesPermanentlyInvalid, pset.nBasesTemporarilyInvalid, pset.nConditions, pset.nAlign, dataSourceStr);
            end
            tcprintf('inline', '{yellow}Dataset: {none}%s\n', pset.datasetName);
            
            tcprintf('inline', '{yellow}Data Source Mode: {none}%s\n', pset.describeDataManualMode());

            if pset.simultaneous
                tcprintf('inline', '{yellow}Simultaneous: {none}%d trials {red}(%d valid)\n', pset.nTrials, pset.nTrialsValid);
            end
            tcprintf('inline', '{yellow}Trial Averaging: {none}at least %d and %g%% of trials for average\n', ...
                pset.minTrialsForTrialAveraging, pset.minFractionTrialsForTrialAveraging*100);
            tcprintf('inline', '{yellow}Time Delta: {none}%g ms, {yellow}Spike Filter: {none}%s\n', ...
                pset.timeDelta, pset.spikeFilter.getDescription);
            if pset.ignoreAllZeroSpikeTrials
                zeroSpikeMode = 'ignore all';
            elseif pset.ignoreLeadingTrailingZeroSpikeTrials
                zeroSpikeMode = 'ignore leading and trailing';
            else
                zeroSpikeMode = 'keep all';
            end
            tcprintf('inline', '{yellow}Zero Spike Trial Mode: {none}%s\n', zeroSpikeMode);
%             if pset.hasDataRandomized
%                 tcprintf('inline', '{yellow}Data Randomized: {none}%d random samples, {red}%s\n', ...
%                     pset.nRandomSamples, pset.conditionDescriptorRandomized.randomizationDescription);
%             end

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

        function explain(pset, varargin)
            p = inputParser();
            p.addParameter('dropFraction', 0.05, @isscalar);
            p.parse(varargin{:});

            temp = pset.tMinForDataMean; %#ok<NASGU>

            hcprintf('%d / %d bases are valid\n  {999}[.basisValid, .nBasesValid, .nBases, .basisInvalidCause]\n', pset.nBasesValid, pset.nBases);
            hcprintf('%d / %d bases are permanently invalid\n{999}[.basisValidPermanent, .nBasesValidPermanentlyInvalid]\n', pset.nBasesPermanentlyInvalid, pset.nBases);
            hcprintf('%d / %d bases are temporarily invalid\n{999}[.basisValidTemporary, .nBasesValidTemporarilyInvalid]\n', pset.nBasesTemporarilyInvalid, pset.nBases);
            hcprintf('%d / %d conditions marked for inclusion\n  {999}[.conditionIncludeMask, .nConditions]\n', pset.nConditions, nnz(pset.conditionIncludeMask));

            hcprintf('\nExistence of valid trial averages:\n');
            hcprintf('  %d / %d bases are marked temporarily invalid as they lack a trial average for all align/conditions.\n', ...
                nnz(pset.basesMissingTrialAverageForNonEmptyConditionAligns & pset.basisValidPermanent & ~pset.basisValid), pset.nBases);
            hcprintf('   %d / %d valid bases lack a trial average for some align/conditions where other bases have trial averages\n  .  {999}These will be ignored when computing .dataMean and trial averaged time windows\n    .basesMissingTrialAverageForNonEmptyConditionAligns, findBasesMissingTrialAverageForNonEmptyConditionAligns]\n', ...
                nnz(pset.basesMissingTrialAverageForNonEmptyConditionAligns & pset.basisValid), pset.nBases);
            if nnz(pset.basesMissingTrialAverageForNonEmptyConditionAligns) > 0
                hcprintf('    {FF9}Consider using {999}.markBasesPermanentlyInvalidMissingTrialAverages{FF9} to invalidate these bases.\n');
            end
            hcprintf('  Trial averages require at least %d trials and at least %d %% of trials at each timepoint\n    {999}[minTrialsForTrialAveraging, minFractionTrialsForTrialAveraging]\n', ...
                pset.minTrialsForTrialAveraging, pset.minFractionTrialsForTrialAveraging * 100);

            hcprintf('  %d / %d conditions have valid trial averages for all aligns, non-empty bases\n    {999}[.conditionHasValidTrialAverageAllAlignsBases]\n', ...
                nnz(pset.conditionHasValidTrialAverageAllAlignsBases), pset.nConditions);

            hcprintf('\nValid time windows for dataMean: hypotheticals for ''dropFraction'' = %g%% \n  {999}[.tMinForDataMean, .tMaxForDataMean]\n', ...
                p.Results.dropFraction * 100);
            listByAlign = pset.listBasesConstrainingTimeWindowValidByAlign('dropFraction', p.Results.dropFraction);
            listByAlignConditionMin = pset.listBasesConstrainingTimeWindowValidByAlignCondition('dropFraction', p.Results.dropFraction, 'mode', 'min');
            listByAlignConditionMax = pset.listBasesConstrainingTimeWindowValidByAlignCondition('dropFraction', p.Results.dropFraction, 'mode', 'max');

            tMinValidAC = pset.tMinValidAllBasesByAlignCondition;
            tMaxValidAC = pset.tMaxValidAllBasesByAlignCondition;
            tMinValidABC = pset.tMinValidByAlignBasisConditionRaw;
            tMaxValidABC = pset.tMaxValidByAlignBasisConditionRaw;

            for iA = 1:pset.nAlign
                if ~isempty(listByAlign{iA})
                    [newMinAC, newMaxAC] = pset.computeNewTimeWindowValidAfterInvalidatingBases(listByAlign{iA});
                    newMinA = nanmin(newMinAC, [], 2);
                    newMaxA = nanmax(newMaxAC, [], 2);
                    hcprintf('  Align %d : [%g - %g %s], would be [%g - %g %s] without %d bases\n', ...
                        iA, pset.tMinForDataMean(iA), pset.tMaxForDataMean(iA), pset.timeUnitName, ...
                        newMinA(iA), newMaxA(iA), pset.timeUnitName, numel(listByAlign{iA}));
                else
                    hcprintf('  Align %d : [%g - %g %s] not constrained by any bases\n', ...
                        iA, pset.tMinForDataMean(iA), pset.tMaxForDataMean(iA), pset.timeUnitName);
                end
                hcprintf('    {999}[.listBasesConstrainingTimeWindowValidByAlign, .computeNewTimeWindowValidAfterInvalidatingBases]\n');

                for iC = 1:pset.nConditions
                    if ~isempty(listByAlignConditionMin{iA, iC}) ||  ~isempty(listByAlignConditionMax{iA, iC})
                        listCombined = unique(cat(1, listByAlignConditionMin{iA, iC}, listByAlignConditionMax{iA, iC}));
                        [newMinAC, newMaxAC] = pset.computeNewTimeWindowValidAfterInvalidatingBases(listCombined);
                        hcprintf('    Condition %d : [%g - %g %s], would be [%g - %g %s] without %d bases\n', ...
                            iC, tMinValidAC(iA, iC), tMaxValidAC(iA, iC), pset.timeUnitName, ...
                            newMinAC(iA, iC), newMaxAC(iA, iC), pset.timeUnitName, ...
                            numel(listCombined));
                        if ~isempty(listByAlignConditionMin{iA, iC})
                            list = listByAlignConditionMin{iA, iC};
                            [uniqTimes, ~, whichTimeIdx] = unique(tMinValidABC(iA, list, iC));
                            fprintf('      min ');
                            for iT = 1:numel(uniqTimes)
                                hcprintf('{f99}b%s{999}(%d)', strjoin(list(whichTimeIdx==iT), ','), uniqTimes(iT));
                            end
                            fprintf('\n');
                        end
                        if ~isempty(listByAlignConditionMax{iA, iC})
                            list = listByAlignConditionMax{iA, iC};
                            [uniqTimes, ~, whichTimeIdx] = unique(tMaxValidABC(iA, list, iC));
                            fprintf('      max ');
                            for iT = numel(uniqTimes):-1:1 % sort in reverse order
                                hcprintf('{f99}b%s{999}(%d)', strjoin(list(whichTimeIdx==iT), ','), uniqTimes(iT));
                            end
                            fprintf('\n');
                        end
                    else
                        hcprintf('    Condition %d : [%g - %g %s] not constrained by any bases\n', ...
                            iC, tMinValidAC(iA, iC), tMaxValidAC(iA, iC), pset.timeUnitName);
                    end
                end
                hcprintf('    {999}[.tMinValidByAlignBasisConditionRaw, .tMinValidAllBasesByAlignCondition]\n    {999}[.listBasesConstrainingTimeWindowValidByAlignCondition, .computeNewTimeWindowValidAfterInvalidatingBases]\n');
            end


        end
    end

    % these methods are setters for property values which change the
    % behavior of the PTS. They automatically invalidate downstream cached
    % values that depend on the value of the property being set.
    methods
        function pset = setMinTrialsForTrialAveraging(pset, v)
            pset.warnIfNoArgOut(nargout);
            assert(isscalar(v));

            assert(~pset.dataInfoManual, 'minTrialsForTrialAveraging cannot be changed with dataInfoManual == true');

            if ~isequal(pset.minTrialsForTrialAveraging, v)
                pset.minTrialsForTrialAveraging = v;
                pset = pset.invalidateDerivedProperties('minTrialsForTrialAveraging');
            end
        end

        function pset = setMinFractionTrialsForTrialAveraging(pset, v)
            pset.warnIfNoArgOut(nargout);
            assert(isscalar(v));

            assert(~pset.dataInfoManual, 'minFractionTrialsForTrialAveraging cannot be changed with dataInfoManual == true');

            if ~isequal(pset.minFractionTrialsForTrialAveraging, v)
                pset.minFractionTrialsForTrialAveraging = v;
                pset = pset.invalidateDerivedProperties('minFractionTrialsForTrialAveraging');
            end
        end

        function pset = setIgnoreAllZeroSpikeTrials(pset, tf)
            pset.warnIfNoArgOut(nargout);
            assert(isscalar(tf) && islogical(tf));

            assert(~pset.dataInfoManual, 'ignoreAllZeroSpikeTrials cannot be changed with dataInfoManual == true');

            if ~isequal(pset.ignoreAllZeroSpikeTrials, tf)
                pset.ignoreAllZeroSpikeTrials = tf;
                pset = pset.invalidateDerivedProperties('ignoreAllZeroSpikeTrials');
            end
        end

        function pset = setIgnoreLeadingTrailingZeroSpikeTrials(pset, tf)
            pset.warnIfNoArgOut(nargout);
            assert(isscalar(tf) && islogical(tf));

            assert(~pset.dataInfoManual, 'ignoreLeadingTrailingZeroSpikeTrials cannot be changed with dataInfoManual == true');

            if ~isequal(pset.ignoreLeadingTrailingZeroSpikeTrials, tf)
                pset.ignoreLeadingTrailingZeroSpikeTrials = tf;
                pset = pset.invalidateDerivedProperties('ignoreLeadingTrailingZeroSpikeTrials');
            end
        end

        function pset = setSpikeFilter(pset, f)
            % changing spikeFilter invalidates everything and we reapply
            % all alignDescriptors to get the padding right
            pset.warnIfNoArgOut(nargout);
            assert(isa(f, 'SpikeFilter'));

            assert(~pset.dataInfoManual, 'spikeFilter cannot be changed with dataInfoManual == true');

            if ~isequal(pset.spikeFilter,  f)
                pset.spikeFilter = f;
                pset = pset.invalidateDerivedProperties('spikeFilter');

                % spikeFilters also set padding within alignDescriptorSet
                pset = pset.setAlignDescriptorSet(pset.alignDescriptorSet);
            end
        end

        function pset = setTimeDelta(pset, v)
            pset.warnIfNoArgOut(nargout);
            if pset.spikeFilter.timeDelta ~= v
                sf = pset.spikeFilter;
                sf.timeDelta = v;
                pset = pset.setSpikeFilter(sf);
            end
        end

        function pset = setConditionIncludeMask(pset, v)
            pset.warnIfNoArgOut(nargout);
            assert(islogical(v) && isvector(v) && numel(v) == pset.nConditions, ...
                'conditionIncludeMask must be logical vector with length nConditions');

            % TODO this can be changed to handle manual masking of the appropriate fields for dataInfoManual case
            assert(~pset.dataInfoManual, 'conditionIncludeMask cannot be changed with dataInfoManual == true');

            cd = pset.conditionDescriptor;
            cd = cd.setConditionIncludeMask(v);
            pset = pset.setConditionDescriptor(cd);
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

        function c = listBasesConstrainingTimeWindowValidByAlignCondition(pset, varargin)
            % if tMinByAlign and tMaxByAlign are not specified, they will
            % be picked according to the widest window available for any
            % basis on that align / condition
            p = inputParser();

            % specify either these manually
            p.addParameter('tMinByAlignCondition', nan(pset.nAlign, pset.nConditions), @(x) isequal(size(x), [pset.nAlign pset.nConditions]));
            p.addParameter('tMaxByAlignCondition', nan(pset.nAlign, pset.nConditions), @(x) isequal(size(x), [pset.nAlign pset.nConditions]));

            % proceed outward from both ends or just to widen the min or
            % max?
            p.addParameter('mode', 'both', @(x) ismember(x, {'both', 'min', 'max'}));

            % or use this to allow only a certain percent of the
            p.addParameter('dropFraction', 0.05, @isscalar);
            p.parse(varargin{:});

            % changing 20160728
            %             cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            cMask = pset.conditionIncludeMask;
            c = cell(pset.nAlign, pset.nConditions);

            tMinValidABC = pset.tMinValidByAlignBasisConditionRaw;
            tMaxValidABC = pset.tMaxValidByAlignBasisConditionRaw;
            basisValid = pset.basisValid; %#ok<*PROPLC>

            widenMin = ismember(p.Results.mode, {'both', 'min'});
            widenMax = ismember(p.Results.mode, {'both', 'max'});
            if widenMax && widenMin
                dropFraction = p.Results.dropFraction / 2;
            else
                dropFraction = p.Results.dropFraction;
            end

            for iA = 1:pset.nAlign
                for iC = 1:pset.nConditions
                    if cMask(iC)
                        % if no time window specified, allow only a certain
                        % fraction of bases to be listed as constrained

                        tMin = p.Results.tMinByAlignCondition(iA,iC);
                        tMax = p.Results.tMaxByAlignCondition(iA,iC);

                        if isnan(tMin)
                            tMin = quantile(tMinValidABC(iA, pset.basisValid, iC), 1-dropFraction, 2);
                        end
                        if isnan(tMax)
                            tMax = quantile(tMaxValidABC(iA, pset.basisValid, iC), dropFraction, 2);
                        end

                        % find bases which can't provide the full tMin:tMax window
                        mask = falsevec(pset.nBases);
                        if widenMin
                            mask = mask | tMinValidABC(iA, :, iC)' > tMin;
                        end
                        if widenMax
                            mask = mask | tMaxValidABC(iA, :, iC)' < tMax;
                        end
                        mask = mask & basisValid;
                        c{iA, iC} = find(mask);
                    else
                        % some bases have no valid window on this align,
                        % condition. find those bases
                        c{iA, iC} = find(makecol(squeeze(pset.hasValidTrialAverageByAlignBasisCondition(iA, :, iC))) & basisValid);
                    end
                end
            end
        end

        function [tMinValidAllBasesByAlignCondition, tMaxValidAllBasesByAlignCondition] = ...
                computeNewTimeWindowValidAfterInvalidatingBases(pset, basesMarkInvalid)
            % changed 20160728
            %             cMask = pset.conditionHasValidTrialAverageAllAlignsBases(:);
            cMask = pset.conditionIncludeMask;

            if ~any(cMask)
                tMinValidAllBasesByAlignCondition = nan(pset.nAlign, pset.nConditions);
                tMaxValidAllBasesByAlignCondition = nan(pset.nAlign, pset.nConditions);
            else
                tMinABC = pset.tMinValidByAlignBasisConditionRaw;
                tMaxABC = pset.tMaxValidByAlignBasisConditionRaw;
                basisValid = pset.basisValid;
                basisValid(basesMarkInvalid) = false;
                % but then reexpand this to have the full complement of
                % conditions, using NaNs for non-contributing conditions
                inflate = @(v) TensorUtils.squeezeDims(TensorUtils.inflateMaskedTensor(v, 3, cMask, NaN), 2);
                tMinValidAllBasesByAlignCondition = inflate(nanmax(tMinABC(:, basisValid, cMask), [], 2));
                tMaxValidAllBasesByAlignCondition = inflate(nanmin(tMaxABC(:, basisValid, cMask), [], 2));
            end
        end

        function list = listBasesConstrainingTimeWindowValid(pset, varargin)
            c = pset.listBasesConstrainingTimeWindowValidByAlign(varargin{:});
            list = unique(cat(1, c{:}));
        end

        function listByAlign = listBasesConstrainingTimeWindowValidByAlign(pset, varargin)
            p = inputParser();

            % either specify these
            p.addParameter('tMinByAlign', nanvec(pset.nAlign), @(x) isvector(x) && numel(x) == pset.nAlign);
            p.addParameter('tMaxByAlign', nanvec(pset.nAlign), @(x) isvector(x) && numel(x) == pset.nAlign);

            % or specify tMinByAlignCondition and tMaxBYALignCondition

            % or use this to allow only a certain percent of the
            p.addParameter('dropFractionPerCondition', 0.05, @isscalar);

            p.parse(varargin{:});

            c = pset.listBasesConstrainingTimeWindowValidByAlignCondition(...
                'tMinByAlignCondition', repmat(p.Results.tMinByAlign(:), 1, pset.nConditions), ...
                'tMaxByAlignCondition', repmat(p.Results.tMaxByAlign(:), 1, pset.nConditions), ...
                'dropFraction', p.Results.dropFractionPerCondition, ...
                p.Unmatched);

            listByAlign = cellvec(pset.nAlign);
            for iA = 1:pset.nAlign
                listByAlign{iA} = unique(cat(1, c{iA, :}));
            end
        end
    end

    % methods which set and apply the conditionDescriptor
    methods
        function pset = setConditionDescriptor(pset, cd, varargin)
            % update the conditionDescriptor for all bases within
            % and clear any generated condition average data
            p = inputParser();
            p.addParameter('suppressWarning', false, @islogical);
            p.parse(varargin{:});

            pset.warnIfNoArgOut(nargout);
            assert(ismember(class(cd), {'ConditionDescriptor', 'ConditionInfo'}), ...
                'Must be ConditionDescriptor instance');

            if isa(cd, 'ConditionInfo')
                cd = cd.fixAllValueLists().getConditionDescriptor();
            end

            assert(cd.allAxisValueListsManual, 'ConditionDescriptor must have manual axis value lists. Use .setAxisValueList or .fixValueListsByApplyingToTrialData');

            pset.conditionDescriptor = cd.getConditionDescriptor();

            if pset.dataInfoManual
                assert(pset.conditionDescriptor.nConditions == cd.nConditions, ...
                    'Pset has manual data source, conditionDescriptor can only be replaced so as to preserve nConditions');
                if ~p.Results.suppressWarning
                    warning('Replacing ConditionDescriptor for PopulationTrajectorySet with manual data source. This will only change the labeling of conditions, not the extracted data');
                end

                % need to replace condition descriptor in alignsummaryData
                % too
                for i = 1:numel(pset.alignSummaryData)
                    pset.alignSummaryData{i} = pset.alignSummaryData{i}.setConditionDescriptor(cd);
                end
            else
                prog = ProgressBar(pset.nDataSources, 'Grouping trials in data sources');
                conditionDescriptor = pset.conditionDescriptor;
                for iSrc = 1:pset.nDataSources
                    pset.dataSources{iSrc} = pset.dataSources{iSrc}.setConditionDescriptor(conditionDescriptor);
                    prog.update(iSrc);
                end
                prog.finish();

                pset = pset.invalidateDerivedProperties('conditionDescriptor');
            end
            prog.finish();
        end

        function pset = flattenConditionAxes(pset)
            pset.warnIfNoArgOut(nargout);
            pset = pset.setConditionDescriptor(pset.conditionDescriptor.flattenAxes(), 'suppressWarning', true);
        end

        function pset = setConditionNames(pset, varargin)
            % setConditionNames(names, [namesShort])
            pset.warnIfNoArgOut(nargout);
            pset.conditionDescriptor = pset.conditionDescriptor.setConditionNames(varargin{:});
        end
    end
    
    methods(Access=protected)
        function pset = internalRealignDataSources(pset)
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
        end
    end

    % alignment: methods which directly update the align descriptor set of this pset
    methods
        function pset = setAlignDescriptorSet(pset, adSet)
            pset.warnIfNoArgOut(nargout);

            assert(~pset.dataInfoManual, 'PopulationTrajectorySets with manual data source cannot be realigned');

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

            pset = pset.internalRealignDataSources();

            % changing alignments invalidates everything
            pset = pset.invalidateDerivedProperties({'alignDescriptorSet', 'dataSources'});
        end

        function pset = selectAlign(pset, mask)
            % keep only align listed in or selected by mask
            pset.warnIfNoArgOut(nargout);

            mask = TensorUtils.vectorIndicesToMask(mask, pset.nAlign);

            % apply slicing operation along N dim on all properties
            % note that this will take care of translationNormalization automatically through customDim handling
            pset = pset.transformInternalProperties_selectAlongDimension('A', mask);
            
            % now alignDescriptor should have been masked
%             if ~pset.dataMeanManual
            pset = pset.invalidateDerivedProperties({'tMinValidByAlignBasisConditionRaw', 'tMaxValidByAlignBasisConditionRaw'});
%             end
        end
    end

    % manual time reslicing: converts to manual data source pset and
    % selects along time axis
    methods
        function psetManual = getAsManualDataByTrial(pset, varargin)
            if pset.dataByTrialManual
                psetManual = pset;
            else
                psetManual = PopulationTrajectorySetBuilder.convertToManualDataByTrial(pset, vararin{:});
            end
        end
        
        function psetManual = getAsManualDataInfo(pset)
            if pset.dataInfoManual
                psetManual = pset;
            else
                psetManual = PopulationTrajectorySetBuilder.convertToManualDataInfo(pset, varargin{:});
            end
        end
        
        function psetManual = getAsManualDataMean(pset, varargin)
            if pset.dataMeanManual
                psetManual = pset;
            else
                psetManual = PopulationTrajectorySetBuilder.convertToManualDataMean(pset, vararin{:});
            end
        end

        function psetSliced = manualSliceOrExpandTimeWindow(pset, tMinByAlign, tMaxByAlign)
            % convert to a manual pset, then manually slice timepoints
            % if tMinByAlign and tMaxByAlign exceed the current time
            % bounds, the time window will be expanded to reflect this

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
                tvecCell = pset.tvecDataByTrial;
                for iBasis = 1:pset.nBases
                    prog.update(iBasis);
                    for iAlign = 1:pset.nAlign

                        % slice data by trial
                        tvecCurrent = tvecCell{iBasis, iAlign};
                        tvecNew = tMinByAlign(iAlign):pset.timeDelta:tMaxByAlign(iAlign);
                        pb.dataByTrial{iBasis, iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataByTrial{iBasis, iAlign}, tvecNew, 2);

                        pb.tMinForDataByTrial(iBasis, iAlign) = tMinByAlign(iAlign);
                        pb.tMaxForDataByTrial(iBasis, iAlign) = tMaxByAlign(iAlign);

                        % also adjust the tMinByTrial and tMaxByTrial of valid
                        % time points each trial
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
                pb.tMinValidByAlignBasisConditionRaw(iAlign, :, :) = max(tMinByAlign(iAlign), ...
                    pb.tMinValidByAlignBasisConditionRaw(iAlign, :, :));
                pb.tMaxValidByAlignBasisConditionRaw(iAlign, :, :) = min(tMaxByAlign(iAlign), ...
                    pb.tMaxValidByAlignBasisConditionRaw(iAlign, :, :));

                pb.tMinForDataMean(iAlign) = tMinByAlign(iAlign);
                pb.tMaxForDataMean(iAlign) = tMaxByAlign(iAlign);
            end

            % slice dataByMean
            %             tmaskByAlign = cellvec(pset.nAlign);
            for iAlign = 1:pset.nAlign
                % build masks used below to slice dataMean
                tvecCurrent = pset.tvecDataMean{iAlign};
                tvecNew = tMinByAlign(iAlign):pset.timeDelta:tMaxByAlign(iAlign);
                pb.dataMean{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataMean{iAlign}, tvecNew, 3);
                pb.dataSem{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataSem{iAlign}, tvecNew, 3);

%                 if pset.hasDataRandomized
%                     pb.dataMeanRandomized{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataMeanRandomized{iAlign}, tvecNew, 3);
%                     pb.dataSemRandomized{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataSemRandomized{iAlign}, tvecNew, 3);
%                 end
            end

            % TODO
%             % slice dataDifferenceOfTrialsScaledNoiseEstimate and dataDifferenceOfTrialsScaledNoiseEstimateRandomized
%             if ~isempty(pb.dataDifferenceOfTrialsScaledNoiseEstimate)
%                 for iAlign = 1:pset.nAlign
%                     % build masks used below to slice dataMean
%                     tvecCurrent = pset.tvecDataMean{iAlign};
%                     tvecNew = tMinByAlign(iAlign):pset.timeDelta:tMaxByAlign(iAlign);
%                     pb.dataDifferenceOfTrialsScaledNoiseEstimate{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataDifferenceOfTrialsScaledNoiseEstimate{iAlign}, tvecNew, 2);
%                     if pset.hasDataRandomized
%                         pb.dataDifferenceOfTrialsScaledNoiseEstimateRandomized{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataDifferenceOfTrialsScaledNoiseEstimateRandomized{iAlign}, tvecNew, 2);
%                     end
%                 end
%             end
%
%             % slice dataCachedSampledTrial and dataCachedMeanExcludingSampledTrials
%             if ~isempty(pb.dataCachedSampledTrialsTensor)
%                 for iAlign = 1:pset.nAlign
%                     tvecCurrent = pset.tvecDataMean{iAlign};
%                     tvecNew = tMinByAlign(iAlign):pset.timeDelta:tMaxByAlign(iAlign);
%                     pb.dataCachedSampledTrials{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataCachedSampledTrials{iAlign}, tvecNew, 2);
%                     pb.dataCachedMeanExcludingSampledTrials{iAlign} = TensorUtils.sliceOrExpandToAlignTimeVector(tvecCurrent, pb.dataCachedMeanExcludingSampledTrials{iAlign}, tvecNew, 2);
%                 end
%             end

            if pset.hasDataByTrial
                psetSliced = pb.buildManualWithSingleTrialData();
            else
                psetSliced = pb.buildManualWithTrialAveragedData();
            end
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

        % function pset = setInterAlignGap(pset, gaps)
        %     % set .interAlignGaps, which represent the time gaps between
        %     % successive alignments, mainly when plotting
        %     pset.warnIfNoArgOut(nargout);
        %
        %     if pset.nAlign < 2
        %         error('Inter alignment gap not valid when only one align present');
        %     end
        %     if isscalar(gaps)
        %         gaps = repmat(gaps, pset.nAlign - 1, 1);
        %     else
        %         assert(numel(gaps) == pset.nAlign - 1, 'Gaps must be scalar or be length nAlign-1');
        %     end
        %
        %     pset.interAlignGaps = gaps;
        %
        %     % this is for convenience / avoiding confusion
        %     prog = ProgressBar(pset.nDataSources, 'Updating condition appearanceFn in data sources');
        %     for iSrc = 1:pset.nDataSources
        %         pset.dataSources{iSrc} = pset.dataSources{iSrc}.setInterAlignGap(gaps);
        %         prog.update(iSrc);
        %     end
        %     prog.finish();
        % end

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

            function data = doTransform(data, prop, propMeta, varargin) %#ok<INUSL>
                attr = propMeta.attr;
                if isfield(attr, 'translate') && attr.translate
                    if isfield(attr, 'normalize') && attr.normalize
                        data = propMeta.applyFnToLastLevel(data, @trNorm.applyTranslationNormalizationToData);
                    else
                        error('Translation without normalization not supported');
                    end
                elseif isfield(attr, 'normalize') && attr.normalize
                    data = propMeta.applyFnToLastLevel(data, @trNorm.applyNormalizationToData);
                end
            end

            pset = pset.transformInternalProperties(@doTransform);
        end

        function pset = clearTranslationNormalization(pset)
            pset.warnIfNoArgOut(nargout);

            if isempty(pset.translationNormalization)
                % No translation / normalization applied
                return;
            end

            trNorm = pset.translationNormalization;

            function data = undoTransform(data, prop, propMeta, varargin) %#ok<INUSL>
                attr = propMeta.attr;
                if isfield(attr, 'translate') && attr.translate
                    if isfield(attr, 'normalize') && attr.normalize
                        data = propMeta.applyFnToLastLevel(data, @trNorm.undoTranslationNormalizationToData);
                    else
                        error('Translation without normalization not supported');
                    end
                elseif isfield(attr, 'normalize') && attr.normalize
                    data = propMeta.applyFnToLastLevel(data, @trNorm.undoNormalizationToData);
                end
            end

            pset = pset.transformInternalProperties(@undoTransform);

            pset.translationNormalization = [];
        end

        function [pset, tn] = translateNormalize(pset, offsetByBasis, normalizationByBasis, varargin)
            pset.warnIfNoArgOut(nargout);

            if any(normalizationByBasis == 0)
                warning('Replacing normalization by 0 with 1');
                normalizationByBasis(normalizationByBasis == 0) = 1;
            end
            tn = StateSpaceTranslationNormalization.buildManual(offsetByBasis, normalizationByBasis, varargin{:});
            pset = pset.applyTranslationNormalization(tn);
        end

        function [pset, tn] = translate(pset, offsetByBasis, varargin)
            pset.warnIfNoArgOut(nargout);
            [pset, tn] = pset.translateNormalize(offsetByBasis, onesvec(pset.nBases), varargin{:});
        end

        function [pset, tn] = normalize(pset, normalizationByBasis, varargin)
            pset.warnIfNoArgOut(nargout);
            [pset, tn] = pset.translateNormalize(zerosvec(pset.nBases), normalizationByBasis, varargin{:});
        end

        function [pset, tn] = meanSubtractBases(pset, varargin)
            pset.warnIfNoArgOut(nargout);
            [pset, tn] = pset.translate(-pset.computeMeanByBasis(varargin{:}), 'translationDescription', 'mean-subtracted');
        end

        function [pset, tn] = normalizeBasesByStd(pset, varargin)
            p = inputParser();
            p.addParameter('denominatorOffset', 0, @isscalar); % x = x / (std(x) + offset)
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pset.warnIfNoArgOut(nargout);
            if p.Results.denominatorOffset ~= 0
                desc = 'soft std-normalized';
            else
                desc = 'std-normalized';
            end
            [pset, tn] = pset.normalize(1 ./ (pset.computeStdByBasis(p.Unmatched) + p.Results.denominatorOffset), 'normalizationDescription', desc);
        end

        function [pset, tn] = normalizeBasesByRange(pset, varargin)
            p = inputParser();
            p.addParameter('denominatorOffset', 0, @isscalar); % x = x / (std(x) + offset)
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pset.warnIfNoArgOut(nargout);
            if p.Results.denominatorOffset ~= 0
                desc = 'soft range-normalized';
            else
                desc = 'range-normalized';
            end

            [pset, tn] = pset.normalize(1 ./ (pset.computeRangeByBasis(p.Unmatched) + p.Results.denominatorOffset), 'normalizationDescription', desc);
        end

        function [pset, tn] = zscoreByBasis(pset, varargin)
            pset.warnIfNoArgOut(nargout);
            [pset, tn] = pset.translateNormalize(-pset.computeMeanByBasis(varargin{:}), pset.computeStdByBasis(varargin{:}), ...
                'translationDescription', 'mean-subtracted', ...
                'normalizationDescription', 'std-normalized');
        end

        function pset = normalizeBasesSoftRange(pset, alpha, varargin)
            pset.warnIfNoArgOut(nargout);

            tr = SoftRangeNormalization.buildFromPopulationTrajectorySet(pset, alpha, varargin{:});
            pset = pset.applyTranslationNormalization(tr);
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
    end

    methods % related to build dataMeanRandomized in various ways

        function pset = storeDataMeanResampleTrialsWithinConditions(pset, varargin)
            % see computeDataMeanResampleTrialsWithinConditions for args,
            % this version stores the results in pset.dataMeanRandomized
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();

            [pset.dataMeanRandomized, pset.dataSemRandomized, ...
                pset.conditionDescriptorRandomized, pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, pset.dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanResampleTrialsWithinConditions(varargin{:});
        end

        function [dataMeanRandomized, dataSemRandomized, cd, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                computeDataMeanResampleTrialsWithinConditions(pset, varargin)
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.addParameter('useValidTimepointsFromOriginal', true, @islogical); % this makes sense to default true here
            p.addParameter('normalizeToMaintainTotalVariance', true, @islogical); % this makes sense to default true here
            p.addParameter('translateToMaintainTotalVariance', true, @islogical);
            p.parse(varargin{:});

            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;

            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            listByConditionCellOriginal = cell(pset.nBases, pset.nConditions);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.resampleTrialsWithinConditions();
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);

                listOriginal = src.conditionInfo.listByCondition(:);
                listByConditionCellOriginal(basisIdx, :) = repmat(listOriginal, [numel(basisIdx), 1, 1]);
            end
            prog.finish();

            [dataMeanRandomized, dataSemRandomized, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell, listByConditionCellOriginal, ...
                'normalizeToMaintainTotalVariance', p.Results.normalizeToMaintainTotalVariance, ...
                'translateToMaintainTotalVariance', p.Results.translateToMaintainTotalVariance, ...
                'useValidTimepointsFromOriginal', p.Results.useValidTimepointsFromOriginal);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function pset = storeDataMeanAxisResampleFromSpecifiedValueListIndices(pset, varargin)
            % see computeDataMeanAxisResampleFromSpecifiedValueListIndices
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace)
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();

            [pset.dataMeanRandomized, pset.dataSemRandomized, ...
                pset.conditionDescriptorRandomized, pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, pset.dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanAxisResampleFromSpecifiedValueListIndices(varargin{:});
        end

        function [dataMeanRandomized, dataSemRandomized, cd, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                computeDataMeanAxisResampleFromSpecifiedValueListIndices(pset, ...
                axisIdxOrAttr, resampleFromList, replace, varargin)
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;

            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            listByConditionCellOriginal = cell(pset.nBases, pset.nConditions);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisResampleFromSpecifiedValueListIndices(axisIdxOrAttr, resampleFromList, replace);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);

                listOriginal = src.conditionInfo.listByCondition(:);
                listByConditionCellOriginal(basisIdx, :) = repmat(listOriginal, [numel(basisIdx), 1, 1]);
            end
            prog.finish();

            [dataMeanRandomized, dataSemRandomized, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell, listByConditionCellOriginal, ...
                p.Unmatched);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function pset = storeDataMeanAxisResampleFromSpecifiedValues(pset, varargin)
            % see computeDataMeanAxisResampleFromSpecifiedValuesfal
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace)
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();

            [pset.dataMeanRandomized, pset.dataSemRandomized, ...
                pset.conditionDescriptorRandomized, pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, pset.dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanAxisResampleFromSpecifiedValues(varargin{:});
        end

        function [dataMeanRandomized, dataSemRandomized, cd, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                computeDataMeanAxisResampleFromSpecifiedValues(pset, ...
                axisIdxOrAttr, resampleFromList, replace, varargin)
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;

            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            listByConditionCellOriginal = cell(pset.nBases, pset.nConditions);
            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisResampleFromSpecifiedValues(axisIdxOrAttr, resampleFromList, replace);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);

                listOriginal = src.conditionInfo.listByCondition(:);
                listByConditionCellOriginal(basisIdx, :) = repmat(listOriginal, [numel(basisIdx), 1, 1]);
            end
            prog.finish();

            [dataMeanRandomized, dataSemRandomized, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell, listByConditionCellOriginal, ...
                p.Unmatched);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function pset = storeDataMeanAxisShuffle(pset, varargin)
            % see computeDataMeanAxisShuffled
            % for arguments
            % this version stores the results in pset.dataMeanRandomized
            % calls into: axisResampleFromSpecifiedValueListIndices(ci, axisIdxOrAttr, resampleFromList, replace)
            pset.warnIfNoArgOut(nargout);
            pset.odc = pset.odc.copy();

            [pset.dataMeanRandomized, pset.dataSemRandomized, ...
                pset.conditionDescriptorRandomized, pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized, pset.dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanAxisShuffle(varargin{:});
        end

        function [dataMeanRandomized, dataSemRandomized, cd, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ....
                computeDataMeanAxisShuffle(pset, ...
                axisIdxOrAttr, varargin)
            p = inputParser();
            p.addParameter('nRandomSamples', 100, @isscalar);
            p.addParameter('randomSeed', pset.conditionDescriptor.randomSeed, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            seed = p.Results.randomSeed;
            nRandomSamples = p.Results.nRandomSamples;

            listByConditionCell = cell(pset.nBases, pset.nConditions, nRandomSamples);
            listByConditionCellOriginal = cell(pset.nBases, pset.nConditions);

            prog = ProgressBar(pset.nDataSources, 'Computing resampled from specified condition info by data source');
            for iSrc = 1:pset.nDataSources
                prog.update(iSrc);
                src = pset.dataSources{iSrc};
                ci = src.conditionInfo.axisShuffle(axisIdxOrAttr);
                list = shiftdim(ci.generateMultipleRandomizedListByCondition(nRandomSamples, ...
                    'initialSeed', seed, 'showProgress', p.Results.nRandomSamples > 100), -1);
                basisIdx = pset.getBasisIdxForDataSource(iSrc);
                listByConditionCell(basisIdx, :, :) = repmat(list, [numel(basisIdx), 1, 1]);

                listOriginal = src.conditionInfo.listByCondition(:);
                listByConditionCellOriginal(basisIdx, :) = repmat(listOriginal, [numel(basisIdx), 1, 1]);
            end
            prog.finish();

            [dataMeanRandomized, dataSemRandomized, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                pset.computeDataMeanUsingMultipleListByCondition(listByConditionCell, listByConditionCellOriginal, ...
                p.Unmatched);
            cd = ConditionDescriptor.fromConditionDescriptor(ci);
        end

        function [dataMeanRandomized, dataSemRandomized, dataDifferenceOfTrialsScaledNoiseEstimateRandomized, dataNumTrialsRawRandomized] = ...
                computeDataMeanUsingMultipleListByCondition(pset, listByConditionCell, listByConditionOriginal, varargin)
            % a utility function for recomputing dataMean using a different
            % conditionInfo.listByCondition than the one already present in
            % the data source. This is typically used when shuffling or resampling
            % condition labels for statistical bootstraps.
            %
            % listByConditionCell is a nBases x nConditions x nSamples

            p = inputParser();
            p.addParameter('useValidTimepointsFromOriginal', false, @islogical);
            p.addParameter('normalizeToMaintainTotalVariance', false, @islogical);
            p.addParameter('translateToMaintainTotalVariance', false, @islogical);
            p.parse(varargin{:});

            %             tMinForDataByTrial = pset.tMinForDataByTrial;
            %             tMaxForDataByTrial = pset.tMaxForDataByTrial;
            %             tMinForDataMean = pset.tMinForDataMean;
            %             tMaxForDataMean = pset.tMaxForDataMean;
            %             timeDelta = pset.timeDelta;
            nAlign = pset.nAlign;
            nBases = pset.nBases;
            nConditions = pset.nConditions;
            nTimeDataMean = pset.nTimeDataMean;
            nRandomSamples = size(listByConditionCell, 3);
            minTrialsForTrialAveraging = pset.minTrialsForTrialAveraging;
            minFractionTrialsForTrialAveraging = pset.minFractionTrialsForTrialAveraging;
            %             T = sum(pset.nTimeDataMean);
            %             dataByTrial = pset.dataByTrial;

            % new matrix will contain all data, concatenated by alignment
            [dataMeanRandomized, dataSemRandomized] = deal(cell(nAlign, 1));
            dataNumTrialsRawRandomized = nan(nBases, nConditions);

            dataByTrial = pset.computeDataByTrialCommonTime();

            for iAlign = 1:nAlign
                prog = ProgressBar(pset.nBases, 'Computing re-conditioned trial-averaged data for align %d', iAlign);

                [dataMeanRandomized{iAlign}, dataSemRandomized{iAlign}] = deal(nan(nBases, nConditions, nTimeDataMean(iAlign), nRandomSamples));

                if nRandomSamples == 0
                    continue;
                end

                for iBasis = 1:nBases

                    % don't process invalid bases, leave these as NaNs
                    if ~pset.basisValid(iBasis)
                        continue;
                    end

                    % nConditions x nSamples
                    listByConditionSamples = squeeze(listByConditionCell(iBasis, :, :));

                    % figure out how many trials for averaging we need by condition, based on
                    % the number of trials each condition (assuming all lists
                    % maintain the same number of trials by condition)
                    nTrialsByCondition = cellfun(@numel, listByConditionSamples(:, 1));
                    minTrialsByCondition = max(minTrialsForTrialAveraging, ...
                        ceil(nTrialsByCondition * minFractionTrialsForTrialAveraging));

                    thisDataByTrial = dataByTrial{iBasis, iAlign};
                    [dataMeanRandomizedThisBasis, dataSemRandomizedThisBasis] = deal(nan(nConditions, nTimeDataMean(iAlign), nRandomSamples));
                    byConditionOriginal = cellfun(@(idx) thisDataByTrial(idx, :), ...
                        listByConditionOriginal(iBasis, :)', 'UniformOutput', false);

                    for iSample = 1:nRandomSamples
                        byCondition = cellfun(@(idx) thisDataByTrial(idx,:), ...
                            listByConditionSamples(:, iSample), 'UniformOutput', false);

                        for iCondition = 1:numel(byCondition)
                            % n by t
                            mat = byCondition{iCondition};

                            matOrig = byConditionOriginal{iCondition};
                            if p.Results.useValidTimepointsFromOriginal
                                % use the nTrials by time of the original data so
                                % placement of nans doesn't change
                                nTrialsByTime = sum(~isnan(matOrig), 1);
                            else
                                nTrialsByTime = sum(~isnan(mat), 1);
                            end

                            % clear invalid timesteps
                            mat(:, nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;
                            matOrig(:, nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;

                            m = nanmean(mat, 1);

                            if p.Results.translateToMaintainTotalVariance && numel(m) > 1
                                % maintain at the same level of variance of the
                                % original mean. when operating across
                                % alignments, this means shifting the means
                                % to match the original mean over time
                                meanOrig = nanmean(nanmean(matOrig, 1));
                                mu = nanmean(m);
                                m = m - mu + meanOrig;
                            end

                            if p.Results.normalizeToMaintainTotalVariance && numel(m) > 1
                                % maintain at the same level of variance of the
                                % original mean.
                                meanOrig = nanmean(matOrig, 1);
                                mu = nanmean(m);
                                vOrig = nanvar(meanOrig, 1);
                                vThis = nanvar(m, 1);

                                if vOrig ~= 0 && vThis ~= 0
                                    m = (m - mu) / sqrt(vThis / vOrig) + mu;
                                end

                                multiplier = sqrt(vOrig / vThis);
                            else
                                multiplier = 1;
                            end

                            dataMeanRandomizedThisBasis(iCondition, :, iSample) = m;

                            se = nansem(mat, 1)' * multiplier;
                            se(nTrialsByTime < minTrialsByCondition(iCondition)) = NaN;
                            dataSemRandomizedThisBasis(iCondition, :, iSample) = se;
                        end
                    end

                    dataMeanRandomized{iAlign}(iBasis, :, :, :) = dataMeanRandomizedThisBasis;
                    dataSemRandomized{iAlign}(iBasis, :, :, :) = dataSemRandomizedThisBasis;
                    dataNumTrialsRawRandomized(iBasis, :) = nTrialsByCondition;

                    prog.update(iBasis);
                end
                prog.finish();
            end

            % apply translation normalization if present
            if ~isempty(pset.translationNormalization)
                cellApplyToDataFn = @(data) cellfun(@pset.translationNormalization.applyTranslationNormalizationToData, ...
                    data, 'UniformOutput', false);

                cellApplyNormOnlyToDataFn = @(data) cellfun(@pset.translationNormalization.applyNormalizationToData, ...
                    data, 'UniformOutput', false);
                dataMeanRandomized = cellApplyToDataFn(dataMeanRandomized);
                dataSemRandomized = cellApplyNormOnlyToDataFn(dataSemRandomized);
            end

            % generate difference by trial for randomized data
            if ~pset.hasDataByTrial
                dataDifferenceOfTrialsScaledNoiseEstimateRandomized = [];
            else
                % grab trial counts
                nTrials_NbyC = dataNumTrialsRawRandomized;
                % scale appropriately each timeseries (running along dim 2) taking into account trial counts
                nTrials_Nby1byC = permute(nTrials_NbyC, [1 3 2]);

                % grab random pair of trials for each basis, each condition,
                % concatenate aligns in time. This function does every
                % random sample at once for speed
                debug('Building difference of trials noise estimate for randomized data\n');
                NbyTAbyCby2byRS = pset.arrangeNbyTAbyCbyR('numTrials', 2, ...
                    'nRepeats', nRandomSamples, ...
                    'commonTime', true, ...
                    'sampleMethod', 'randomWithReplacment', ...
                    'ignoreTrialsWithTooFewSamples', true, ...
                    'trialLists', listByConditionCell);

                % diff be N x TA x C x nRandomSamples
                dif_NbyTAbyCbyRS = diff(NbyTAbyCby2byRS, 1, 4);

                % N x TA x C x nRandomSamples (after squeezing out the diff
                % dimension)
                dataDifferenceOfTrialsScaledNoiseEstimateRandomized = ...
                    TensorUtils.squeezeDims(bsxfun(@rdivide, dif_NbyTAbyCbyRS, sqrt(2*nTrials_Nby1byC)), 4);

                % apply translation normalization if present
                if ~isempty(pset.translationNormalization)
                    dataDifferenceOfTrialsScaledNoiseEstimateRandomized = ...
                        pset.translationNormalization.applyNormalizationToData(dataDifferenceOfTrialsScaledNoiseEstimateRandomized);
                end
            end
        end

        function psetNew = withDataRandomSampleAsData(pset, dataRandomIndex, varargin)
            % take random sample idx from dataMeanRandomized as the new
            % dataMean, so that the randomization can be studied as a
            % normal pset
            pset.warnIfNoArgOut(nargout);
            assert(pset.hasDataRandomized == 1);
            cdr = pset.conditionDescriptorRandomized;
            % the ith sample used seed initialSeed+i-1
            cdr = cdr.setRandomSeed(cdr.randomSeed+dataRandomIndex-1);

            b = PopulationTrajectorySetBuilder.copyTrialAveragedOnlyFromPopulationTrajectorySet(pset);
            b.conditionDescriptor = cdr;
            for iAlign = 1:pset.nAlign
                b.dataMean{iAlign} = pset.dataMeanRandomized{iAlign}(:, :, :, dataRandomIndex);
                b.dataSem{iAlign} = pset.dataSemRandomized{iAlign}(:, :, :, dataRandomIndex);
            end
            b.dataNumTrialsRaw = pset.dataNumTrialsRawRandomized;
            b.dataDifferenceOfTrialsScaledNoiseEstimate = pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized(:, :, :, dataRandomIndex);
            b.dataMeanRandomized = [];
            b.dataSemRandomized = [];
            b.dataNumTrialsRawRandomized = [];
            b.conditionDescriptorRandomized = [];
            b.dataDifferenceOfTrialsScaledNoiseEstimateRandomized = [];
            psetNew = b.buildManualWithTrialAveragedData();
        end

        % TODO
%         function pset = dropDataRandom(pset)
%             % drop the .dataRandom and derived fields from the pset
%             pset.warnIfNoArgOut(nargout);
%             pset.dataMeanRandomized = [];
%             pset.dataSemRandomized = [];
%             pset.dataNumTrialsRawRandomized = [];
%             pset.dataDifferenceOfTrialsScaledNoiseEstimateRandomized = [];
%             pset.conditionDescriptorRandomized = [];
%             pset.odc = pset.odc.copy();
%             pset.odc.flushRandomizedTrialAveragedData();
%         end
%
%         function pset = clearCachedSampledTrialsTensor(pset)
%             pset.warnIfNoArgOut(nargout);
%             pset.dataCachedSampledTrials = [];
%             pset.dataCachedSampledTrialCounts = [];
%             pset.dataCachedMeanExcludingSampledTrials = [];
%         end
    end

    methods % Selecting bases
        function pset = selectBases(pset, mask)
            % keep only bases listed in or selected by mask
            pset.warnIfNoArgOut(nargout);

            mask = TensorUtils.vectorIndicesToMask(mask, pset.nBases);

            % apply slicing operation along N dim on all properties
            % note that this will take care of translationNormalization automatically through customDim handling
            pset = pset.transformInternalProperties_selectAlongDimension('N', mask);

            % basisDataSourceIdx has been filtered by mask already
            % but this leaves some elements of dataSources unused
            % so we want to select the elements of dataSources that are used in basisDataSourceIdx,
            % and then update the lookup indices in basisDataSurceIdx to match.
            [pset.basisDataSourceIdx, pset.dataSources, selectedDataSources] = PopulationTrajectorySet.filterUsedUpdateLookup(pset.basisDataSourceIdx, pset.dataSources);
            pset = pset.transformInternalProperties_selectAlongDimension('nDataSources', selectedDataSources);

            % same squashing required with alignSummaryData and basisAlignSummaryLookup
            [pset.basisAlignSummaryLookup, pset.alignSummaryData, selectedAlignSummaryData] = PopulationTrajectorySet.filterUsedUpdateLookup(pset.basisAlignSummaryLookup, pset.alignSummaryData);
            pset = pset.transformInternalProperties_selectAlongDimension('nAlignSummaryData', selectedAlignSummaryData);

            if ~pset.dataMeanManual
                % when data mean is computed from dataByTrial, we should recompute it and all of the time vectors,
                % to allow for the possibility that the time vectors will expand. If not, we leave it as is.
                pset = pset.invalidateDerivedProperties({'dataByTrialRaw'});
            end
        end

        function [pset, mask] = selectBasesMissingTrialAverageForNonEmptyConditionAligns(pset)
            pset.warnIfNoArgOut(nargout);
            mask = ~pset.basesMissingTrialAverageForNonEmptyConditionAligns;
            pset = pset.selectBases(mask);
        end
    end

    methods % Filtering conditions
        function [pset, maskC] = selectConditions(pset, varargin)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);

            [cd, maskC] = pset.conditionDescriptor.selectConditions(varargin{:});
            pset = pset.setConditionDescriptor(cd);
        end

        function [pset, maskC] = selectConditionsAlongAxis(pset, axisAttr, mask)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);

            [cd, maskC] = pset.conditionDescriptor.selectConditionsAlongAxis(axisAttr, mask);

            if ~pset.dataInfoManual
                % just update the condition descriptor and the data will be
                % re-trialaveraged accordingly
                pset = pset.setConditionDescriptor(cd);
            else
                % manually orchestrate a selection of the conditions
                % this will drop any single trial data since these won't be
                % preserved by
                % PopulationTrajectorySetCrossConditionUtilities, this
                % could be fixed if useful.
                aIdx = pset.conditionDescriptor.axisLookupByAttributes(axisAttr);
                nOld = pset.conditionDescriptor.conditionsSize(aIdx);
                inds = TensorUtils.vectorMaskToIndices(mask);
                nNew = numel(inds);
                wNewByOld = zeros(nNew, nOld);
                for iN = 1:nNew
                    wNewByOld(iN, inds(iN)) = 1;
                end

                pset = PopulationTrajectorySetCrossConditionUtilities.applyLinearCombinationAlongConditionAxis(pset, aIdx, wNewByOld);
                pset = pset.setConditionDescriptor(cd);
            end
        end

        function [pset, maskC] = matchSelectConditionsAlongAxis(pset, varargin)
            % select specific conditions by linear index or mask
            pset.warnIfNoArgOut(nargout);

            [cd, maskC] = pset.conditionDescriptor.matchSelectConditionsAlongAxis(varargin{:});
            pset = pset.setConditionDescriptor(cd);
        end
    end

    methods(Access=protected)
        % utility method to subselect conditions from data for a manual
        % data source
        function pset = manualApplyConditionMask(pset, maskC) %#ok<INUSD>
            pset.warnIfNoArgOut(narout);
            error('Not implemented');
        end
    end

    methods % Marking bases as invalid
        function pset = updateValid(pset)
            pset.warnIfNoArgOut(nargout);

            if pset.dataInfoManual
                mask = ~pset.basisValid;
                clearFn = @(in, dim) TensorUtils.assignValueMaskedSelectionAlongDimension(in, dim, mask, NaN);

                % mask out data from every invalid basis
                maskDim1 = {'basisNames', 'basisUnits', ...
                    'dataByTrial', 'tMinForDataByTrial', ...
                    'tMaxForDataByTrial', 'tMinByTrial', 'tMaxByTrial'};
                maskDim2 = {'dataNumTrialsRaw', 'dataNumTrialsRawRandomized', 'tMinValidByAlignBasisConditionRaw', 'tMaxValidByAlignBasisConditionRaw'};
                maskCellDim1 = {'dataMean', 'dataSem', 'dataMeanRandomized', 'dataSemRandomized', ...
                    'dataIntervalHigh', 'dataIntervalLow', ...
                    'dataCachedSampledTrials', 'dataCachedSampledTrialCounts', 'dataCachedMeanExcludingSampledTrials', ...
                    'dataDifferenceOfTrialsScaledNoiseEstimate', 'dataDifferenceOfTrialsScaledNoiseEstimateRandomized'};

                for i = 1:numel(maskDim1)
                    fld = maskDim1{i};
                    if ~isempty(pset.(fld))
                        % some have up to 4 dimensions, add a couple of extra
                        % colons just in case
                        pset.(fld) = clearFn(pset.(fld), 1);
                    end
                end

                for i = 1:numel(maskDim2)
                    fld = maskDim2{i};
                    if ~isempty(pset.(fld))
                        pset.(fld) = clearFn(pset.(fld), 2);
                    end
                end

                for i = 1:numel(maskCellDim1)
                    fld = maskCellDim1{i};
                    if ~isempty(pset.(fld))
                        for j = 1:numel(pset.(fld))
                            pset.(fld){j} = clearFn(pset.(fld){j}, 1);
                        end
                    end
                end

                % filter alignSummaryData
                if ~isempty(pset.alignSummaryData)
                    dsMaskKeep = TensorUtils.vectorIndicesToMask(pset.basisAlignSummaryLookup(pset.basisValid), pset.nAlignSummaryData);
                    pset.alignSummaryData(~dsMaskKeep, :) = {[]};
                end
                pset.alignSummaryAggregated = [];

                % filter translationNormalization
                if ~isempty(pset.translationNormalization)
                    pset.translationNormalization = pset.translationNormalization.setBasesInvalid(mask);
                end
            end
        end

        function pset = markBasesPermanentlyInvalid(pset, mask, cause)
            pset.warnIfNoArgOut(nargout);

            permValid = pset.basisValidPermanent;
            permCause = pset.basisInvalidCausePermanent;

            mask = makecol(TensorUtils.vectorIndicesToMask(mask, pset.nBases));
            maskOrig = mask; % cache for later

            mask = mask & permValid;
            permValid(mask) = false;

            if nargin < 3
                cause = 'unspecified';
            end
            assert(iscellstr(cause) || ischar(cause), 'Cause must be a cellstr or a string');
            if ischar(cause)
                cause = repmat({cause}, nnz(maskOrig), 1);
            end
            cause = makecol(cause);

            if ~any(mask)
                return;
            end

            if numel(cause) == nnz(maskOrig)
                cause = TensorUtils.inflateMaskedTensor(cause, 1, maskOrig, {''});
            end
            assert(numel(cause) == numel(mask), 'Length of cellstr cause must match nnz(mask) or numel(mask)');

            permCause(mask) = cause(mask);

            pset.basisValidPermanent = permValid;
            pset.basisInvalidCausePermanent = permCause;

            % we'll update the properties
            pset = pset.invalidateDerivedProperties({'basisValidPermanent', 'basisInvalidCausePermanent'}, {}, updateFn);

            function value = updateFn(value, prop, propMeta, container)
                % when marking additional bases invalid temporarily, we will clear the values if they are computed on demand
                % if they are manual, we will mark them as NaN
                if strcmp(container, 'pset') && pset.isODCPropUsingManual(prop)
                    value = propMeta.assignValueMaskedSelectionAlongDimByName(value, 'N', pset.basisValidTemporary, NaN, true);
                else
                    value = [];
                end
            end
        end

        function pset = markBasesTemporarilyInvalid(pset, mask, cause)
            pset.warnIfNoArgOut(nargout);

            tempValid = pset.basisValidTemporary;
            tempCause = pset.basisInvalidCauseTemporary;

            mask = makecol(TensorUtils.vectorIndicesToMask(mask, pset.nBases));
            mask = mask & tempValid;
            tempValid(mask) = false;

            if nargin < 3
                cause = 'unspecified';
            end
            assert(iscellstr(cause) || ischar(cause), 'Cause must be a cellstr or a string');
            if ischar(cause)
                cause = repmat({cause}, nnz(mask), 1);
            end
            cause = makecol(cause);
            if numel(cause) == nnz(mask)
                cause = TensorUtils.inflateMaskedTensor(cause, 1, mask, {''});
            end
            assert(numel(cause) == numel(mask), 'Length of cellstr cause must match nnz(mask) or numel(mask)');

            tempCause(mask) = cause(mask);

            pset.basisValidTemporary = tempValid;
            pset.basisInvalidCauseTemporary = tempCause;

            % we'll update the properties
            pset = pset.invalidateDerivedProperties({'basisValidTemporary', 'basisInvalidCauseTemporary'}, {}, updateFn);

            function value = updateFn(value, prop, propMeta, container)
                % when marking additional bases invalid temporarily, we will clear the values if they are computed on demand
                % if they are manual, we will mark them as NaN
                if strcmp(container, 'pset') && pset.isODCPropUsingManual(prop)
                    value = propMeta.assignValueMaskedSelectionAlongDimByName(value, 'N', pset.basisValidTemporary, NaN, true);
                else
                    value = [];
                end
            end
        end

        function pset = restoreBasesTemporarilyInvalid(pset)
            pset.warnIfNoArgOut(nargout);
            if pset.dataInfoManual
                error('Data source has manually-provided data. Bases can only be marked invalid; basisValid cannot be reset');
            end
            pset.basisValidTemporary = truevec(pset.nBases);
            pset.basisInvalidCauseTemporary = cellstrvec(pset.nBases);

        end

        function [pset, mask] = markBasesTemporarilyInvalidMissingTrialAverages(pset)
            pset.warnIfNoArgOut(nargout);
            mask = pset.basesMissingTrialAverageForNonEmptyConditionAligns;
            pset = pset.markBasesTemporarilyInvalid(mask, 'missing trial average for non empty condition aligns');
        end

        function [pset, mask] = markBasesPermanentlyInvalidMissingTrialAverages(pset)
            pset.warnIfNoArgOut(nargout);
            mask = pset.basesMissingTrialAverageForNonEmptyConditionAligns;
            pset = pset.markBasesPermanentlyInvalid(mask, 'missing trial average for non empty condition aligns');
        end

        function warnIfAnyBasesMissingTrialAverageForNonEmptyConditionAligns(pset)
            N = nnz(pset.basesMissingTrialAverageForNonEmptyConditionAligns);
            if N > 0
                warning('%d bases are missing valid trial averages on condition/aligns where other bases have valid trial averages. Use .findBasesMissingTrialAverageForNonEmptyConditionAligns() to identify these and .ignoreBasesMissingTrialAverageForNonEmptyConditionAligns() to mark these invalid', N);
            end
        end

        function idx = findBasesMissingTrialAverageForNonEmptyConditionAligns(pset)
            idx = find(pset.basesMissingTrialAverageForNonEmptyConditionAligns);
        end

        function pset = setBasisNames(pset, basisNames)
            assert(iscellstr(basisNames) && numel(basisNames) == pset.nBases);
            pset.warnIfNoArgOut(nargout);

            pset.basisNames = basisNames;
        end

        function pset = setBasisNamesUsingFormatString(pset, prefix)
%             assert(ischar(prefix) && ~isempty(strfind(prefix, '%d')), 'Format string must be char containing ''%d'''); %#ok<STREMP>
            pset.warnIfNoArgOut(nargout);

            names = arrayfun(@(n) sprintf(prefix, n), (1:pset.nBases)', 'UniformOutput', false);
            pset = pset.setBasisNames(names);
        end
    end

    methods(Static) % Static methods which take multiple psets as args
        function [varargout] = equalizeBasesInvalid(varargin)
            % works for arguments of PopulationTrajectorySet and
            % StateSpaceProjection
            assert(nargin == nargout, 'Number of input and output arguments must match');
            N = nargin;

            if N == 1
                varargout{1} = varargin{1};
                return;
            end

            nBasesVec = nan(numel(varargin), 1);
            for i = 1:numel(varargin)
                if isa(varargin{i}, 'PopulationTrajectorySet')
                    nBasesVec(i) = varargin{i}.nBases;
                elseif isa(varargin{i}, 'StateSpaceProjection')
                    nBasesVec(i) = varargin{i}.nBasesSource;
                else
                    error('Not defined for class %s', class(varargin{i}));
                end
            end
            assert(numel(unique(nBasesVec)) == 1, 'All inputs must have the same nBases');
            nBases = nBasesVec(1);

            mask = truevec(nBases);
            cause = cellstrvec(nBases);

            for i = 1:N
                obj = varargin{i};
                mask = mask & obj.basisValid;
                cause(~obj.basisValid) = cellfun(@(c) sprintf('equalizeBasisInvalid: %s', c),  ...
                    obj.basisInvalidCause(~obj.basisValid), 'UniformOutput', false);
            end

            varargout = cellvec(N);
            for i = 1:N
                obj = varargin{i};
                if isa(obj, 'PopulationTrajectorySet')
                    varargout{i} = obj.markBasesPermanentlyInvalid(~mask, cause);
                elseif isa(obj, 'StateSpaceProjection')
                    varargout{i} = obj.markBasesInvalid(~mask, cause);
                end
            end
        end

        function tf = checkSameBasesValid(varargin)
            nBasesVec = cellfun(@(pset) pset.nBases, varargin);
            assert(numel(unique(nBasesVec)) == 1, 'All inputs must have the same nBases');

            basisValidMat = cell2mat(cellfun(@(pset) pset.basisValid, makerow(varargin), 'UniformOutput', false));

            tf = all(all(basisValidMat, 2) | all(~basisValidMat, 2));
        end

        function [tf, reason] = checkBasesMatch(varargin)
            tf = true;
            reason = '';

            nBases = cellfun(@(pset) pset.nBases, varargin);
            if numel(unique(nBases)) ~= 1
                tf = false;
                reason = 'nBases do not match';
                return;
            end

            basisAlignSummaryLookupCell = cellfun(@(pset) pset.basisAlignSummaryLookup, varargin, 'UniformOutput', false);
            if ~allEqual(basisAlignSummaryLookupCell{:})
                tf = false;
                reason = 'basisAlignSummaryLookup do not match';
                return;
            end

            function tf = allEqual(varargin)
                tf = true;
                for i = 2:numel(varargin)
                    if ~isequal(varargin{1}, varargin{2})
                        tf = false;
                        return;
                    end
                end
            end
        end

        function assertBasesMatch(varargin)
            [tf, reason] = PopulationTrajectorySet.checkBasesMatch(varargin{:});
            assert(tf, reason);
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

        function n = get.nTrialsByCondition(pset)
            if ~pset.simultaneous
                n = nan(pset.conditionsSize);
            else
                n = pset.dataSources{1}.nTrialsByCondition;
            end
        end

        function n = get.nDataSources(pset)
            n = numel(pset.dataSources);
        end

        function n = get.nBases(pset)
            if pset.dataInfoManual
                n = size(pset.dataMean{1}, 1);
            else
                n = numel(pset.basisDataSourceIdx);
            end
        end

        function n = get.nBasesValid(pset)
            n = nnz(pset.basisValid);
        end

        function n = get.nBasesValidPermanent(pset)
            n = nnz(pset.basisValidPermanent);
        end

        function n = get.nBasesPermanentlyInvalid(pset)
            n = nnz(~pset.basisValidPermanent);
        end

        function n = get.nBasesTemporarilyInvalid(pset)
            n = nnz(~pset.basisValidTemporary);
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

        function clists = get.conditionIdxByTrial(pset)
            clists = cellvec(pset.nBases);
            dataSourcesByBasis = pset.dataSources(pset.basisDataSourceIdx);
            for iB = 1:pset.nBases
                idxByCondition = pset.trialLists(iB, :);
                clists{iB} = nanvec(dataSourcesByBasis{iB}.nTrials);
                for iC = 1:numel(idxByCondition)
                    clists{iB}(idxByCondition{iC}) = iC;
                end
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

        % function c = get.alignConditionsWithTrials(pset)
        %    % dataNTrials is nAlign x nBases x nConditions
        %     % c is nAlign x nConditions
        %     c = TensorUtils.squeezeDims(any(pset.dataNTrials(:, pset.basisValid, :), 2), 2);
        % end

        function c = get.alignConditionsWithTrialAverage(pset)
            % c is nAlign x nConditions
            c = TensorUtils.squeezeDims(any(pset.hasValidTrialAverageByAlignBasisCondition(:, pset.basisValid, :), 2), 2);
        end

        function c = get.basisConditionsWithValidTrialAverage(pset)
            % c is nBases x nConditions
            c = TensorUtils.squeezeDims(all(pset.hasValidTrialAverageByAlignBasisCondition, 1), 1);
        end
        %
        function basesMissingTrialAverages = get.basesMissingTrialAverageForNonEmptyConditionAligns(pset)
            % determine which bases missing trials in at least one condition for
            % which other bases have trials in that alignment. basesMask is
            % nBases x 1 logical vector

            % nAlign x nBases x nConditions
            shouldHaveTrialAverage = repmat(permute(pset.alignConditionsWithTrialAverage, [1 3 2]), [1 pset.nBasesValidPermanent 1]);
            hasTrialAverage = pset.hasValidTrialAverageByAlignBasisCondition(:, pset.basisValidPermanent, :);

            % nBases x 1
            basesValidMissingTrialAverages = makecol(TensorUtils.squeezeDims(any(any(shouldHaveTrialAverage & ~hasTrialAverage, 3), 1), [1 3]));
            basesMissingTrialAverages = TensorUtils.inflateMaskedTensor(basesValidMissingTrialAverages, 1, pset.basisValidPermanent, false);
        end

        function v = get.basisValidWithTrialAverageAllNonEmptyConditionAligns(pset)
            v = pset.basisValid & ~pset.basesMissingTrialAverageForNonEmptyConditionAligns;
        end

        function bMask = get.basesNonEmpty(pset)
            % nBases x 1 logical
            bMask = pset.basisValid & makecol(any(any(pset.hasValidTrialAverageByAlignBasisCondition, 1), 3));
        end

        function cMask = get.conditionsWithValidTrialAverageOnNonEmptyBases(pset)
            % look at valid bases that have a trial average on at least one
            % condition. which conditions have trial averages on ALL bases
            % satisfying this?
            bMask = pset.basesNonEmpty;
            % any alignment, all bases
            cMask = squeeze(all(any(pset.hasValidTrialAverageByAlignBasisCondition(:,bMask,:), 1), 2));
        end
    end

    methods % Simple statistics, should be NaN for invalid bases
        % some of the coding here is designed to accept the arguments
        % 'type', 'meanRandom', which will mean the results will operate on
        % dataMeanRandomized, and thus have size nRandomSamples along the
        % last dimension, just like the arrangeCTAbyN-like methods

        function dataStd = computeDataStd(pset)
            % dataNumTrialsRaw is nAlign x nBases x nConditions
            % dataSem is nAlign cell of nBases x nConditions x nTime
            % dataStd will be same size as dataSem

            dataStd = cell(pset.nAlign, 1);
            for iAlign = 1:pset.nAlign
                dataStd{iAlign} = bsxfun(@times, pset.dataSem{iAlign}, sqrt(squeeze(pset.dataNumTrialsRaw(iAlign, :, :))));
            end
        end

        function dataVar = computeDataVar(pset)
            % dataNumTrialsRaw is nAlign x nBases x nConditions
            % dataSem is nAlign cell of nBases x nConditions x nTime
            % dataVar will be same size as dataSemRandomized

            dataVar = cell(pset.nAlign, 1);
            for iAlign = 1:pset.nAlign
                dataVar{iAlign} = bsxfun(@times, pset.dataSem{iAlign}, sqrt(squeeze(pset.dataNumTrialsRaw(iAlign, :, :)))).^2;
            end
        end

        function dataStd = computeDataStdRandomized(pset, varargin)
            % dataNumTrialsRaw is nAlign x nBases x nConditions
            % dataSem is nAlign cell of nBases x nConditions x nTime
            % dataStd will be same size as dataSem

            p = inputParser();
            p.addParameter('dataRandomIndex', 1:pset.nRandomSamples, @isvector);
            p.parse(varargin{:});

            assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');

            dataStd = cell(pset.nAlign, 1);
            for iAlign = 1:pset.nAlign
                dataStd{iAlign} = bsxfun(@times, pset.dataSemRandomized{iAlign}, sqrt(squeeze(pset.dataNumTrialsRawRandomized(iAlign, :, :, p.Results.dataRandomIndex))));
            end
        end

        function dataVar = computeDataVarRandomized(pset, varargin)
            % dataNumTrialsRaw is nAlign x nBases x nConditions
            % dataSem is nAlign cell of nBases x nConditions x nTime
            % dataVar will be same size as dataSemRandomized
            p = inputParser();
            p.addParameter('dataRandomIndex', 1:pset.nRandomSamples, @isvector);
            p.parse(varargin{:});

            assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');

            dataVar = cell(pset.nAlign, 1);
            for iAlign = 1:pset.nAlign
                dataVar{iAlign} = bsxfun(@times, pset.dataSemRandomized{iAlign}, sqrt(squeeze(pset.dataNumTrialsRawRandomized(iAlign, :, :, p.Results.dataRandomIndex)))).^2;
            end
        end

        function snrByBasis = computeSnrByBasis(pset, varargin)
            semCTAbyN = pset.arrangeCTAbyN(varargin{:}, 'type', 'sem');
            noiseByBasis = TensorUtils.squeezeDims(nanmax(semCTAbyN, [], 1), 1);
            rangeByBasis = pset.computeRangeByBasis(varargin{:}, 'type', 'mean');

            snrByBasis = rangeByBasis ./ noiseByBasis;
            snrByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(snrByBasis, 1, ~pset.basisValid, NaN);
        end

        function minByBasis = computeMinByBasis(pset, varargin)
            CTAbyN = pset.arrangeCTAbyN(varargin{:});
            minByBasis = TensorUtils.squeezeDims(nanmin(CTAbyN, [], 1), 1);
            minByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(minByBasis, 1, ~pset.basisValid, NaN);
        end

        function maxByBasis = computeMaxByBasis(pset, varargin)
            CTAbyN = pset.arrangeCTAbyN(varargin{:});
            maxByBasis = TensorUtils.squeezeDims(nanmax(CTAbyN, [], 1), 1);
            maxByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(maxByBasis, 1, ~pset.basisValid, NaN);
        end

        function meanByBasis = computeMeanByBasis(pset, varargin)
            CTAbyN = pset.arrangeCTAbyN(varargin{:});
            meanByBasis = TensorUtils.squeezeDims(nanmean(CTAbyN, 1), 1);
            meanByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(meanByBasis, 1, ~pset.basisValid, NaN);
        end

        function varByBasis = computeVarUncorrectedByBasis(pset, varargin)
            varByBasis = TensorUtils.squeezeDims(nanvar(pset.arrangeCTAbyN(varargin{:}), 0, 1), 1);
            varByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(varByBasis, 1, ~pset.basisValid, NaN);
        end

        function varPerCTA = computeTotalVarUncorrectedPerCTA(pset, varargin)
            % note, this should be the same as doing
            % pset.meanSubtractBases.computeTotalSSPerCTA

            ctaByN = pset.arrangeCTAbyN(varargin{:});
            ctaByN_meanSub = bsxfun(@minus, ctaByN, nanmean(ctaByN, 1));
            varPerCTAByBasis = nanmean(ctaByN_meanSub.^2, 1);
            varPerCTA = squeeze(nansum(varPerCTAByBasis, 2));
        end

        function ssq = computeTotalSS(pset, varargin)
            tensor = pset.arrangeCTAbyN(varargin{:});
            ssq = squeeze(TensorUtils.nansumMultiDim(tensor.^2, [1 2]));
        end

        function ssq = computeTotalSSPerCTA(pset, varargin)
            tensor = pset.arrangeCTAbyN(varargin{:});
            ssq = squeeze(nansum(nanmean(tensor.^2, 1), 2));
        end

        function ssqPer = computeTotalSSPerNCTA(pset, varargin)
            % take squared values averaged over neurons and time
            tensor = pset.arrangeCTAbyN(varargin{:});
            ssqPer = squeeze(TensorUtils.nanmeanMultiDim(tensor.^2, [1 2]));
        end

        function normByBasis = computeNormByBasis(pset, varargin)
            ctaByN = pset.arrangeCTAbyN(varargin{:});
            normByBasis = TensorUtils.squeezeDims(nansum(ctaByN.^2, 1), 1);
            normByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(normByBasis, 1, ~pset.basisValid, NaN);
        end

        function stdByBasis = computeStdByBasis(pset, varargin)
            stdByBasis = TensorUtils.squeezeDims(nanstd(pset.arrangeCTAbyN(varargin{:}), 0, 1), 1);
            stdByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(stdByBasis, 1, ~pset.basisValid, NaN);
        end

        function rangeByBasis = computeRangeByBasis(pset, varargin)
            CTAbyN = pset.arrangeCTAbyN(varargin{:});
            rangeByBasis = TensorUtils.squeezeDims(nanmax(CTAbyN, [], 1) - nanmin(CTAbyN, [], 1), 1);
            rangeByBasis = TensorUtils.assignValueMaskedSelectionAlongDimension(rangeByBasis, 1, ~pset.basisValid, NaN);
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
        function addConditionLabels(pset, varargin)
            pset.conditionDescriptor.addColoredLabels(varargin{:});
        end

        function [offsets, lims] = getAlignPlottingTimeOffsets(pset, tvecCell, varargin)
            % when plotting multiple alignments of data simultaneously,
            % timeseries are plotted side by side along an axis, separated
            % by gaps given by .interAlignGaps. This returns the list of
            % time offsets where each alignment's zero event would be
            % plotted. The time vectors in tvec Cell are used to determine
            % the "width" of each alignments and to find the appropriate
            % location to start the next alignment.
            %
            % tvecCell is a nAlign x 1 cell of:
            %     time vectors, for common time across trials, or
            %     nTrials x 1 cell of time vectors, for different times
            %     across trials.
            %
            % Parameters:
            %     alignIdx : vector indicating which alignments to include

            p = inputParser();
            p.addParameter('alignIdx', 1:pset.nAlign, @isvector);
            p.parse(varargin{:});

            % select alignment indices
            alignIdx = 1:pset.nAlign;
            alignIdx = alignIdx(p.Results.alignIdx);
            nAlign = numel(alignIdx);

            offsets = nan(nAlign, 1);
            offsets(1) = 0;
            currentOffset = 0;

            % compute start/stop of each alignment
            [mins, maxs] = deal(nanvec(nAlign));
            for iAlign = 1:nAlign
                if iscell(tvecCell{iAlign}) && ~isempty(tvecCell{iAlign})
                    mins(iAlign) = nanmin(cellfun(@nanminNanEmpty, tvecCell{iAlign}));
                    maxs(iAlign) = nanmax(cellfun(@nanmaxNanEmpty, tvecCell{iAlign}));
                elseif isempty(tvecCell{iAlign})
                    error('Time cell for alignment %d is empty', iAlign);
                else
                    mins(iAlign) = nanmin(tvecCell{iAlign});
                    maxs(iAlign) = nanmax(tvecCell{iAlign});
                end
            end

            if nAlign == 1
                % no gaps needed
                offsets = 0;
                lims = [mins(1), maxs(1)];
                return;
            end

            % determine the inter alignment gaps
            % if isempty(pset.interAlignGaps)
                % no inter align gap specified, determine automatically as
                % 2% of total span
                T = sum(maxs - mins);
                gaps = repmat(0.02 * T, nAlign-1, 1);
            % else
            %     gaps = pset.interAlignGaps;
            % end

            for iAlign = 2:nAlign
                currentOffset = currentOffset + maxs(iAlign-1) + ...
                    gaps(iAlign-1) - mins(iAlign);
                offsets(iAlign) = currentOffset;
            end

            lims = [mins(1), maxs(end) + offsets(end)];
        end

        function plotBases(pset, varargin)
            % plot bases one above the next
            p = inputParser;
            p.addParameter('individualTrials', false, @islogical);

            p.addParameter('dataRandomIndex', [], @(x) isempty(x) || isscalar(x));

            p.addParameter('alignIdx', 1:pset.nAlign, @isnumeric);
            p.addParameter('basisIdx', [], @(x) isempty(x) || isvector(x));
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
            p.addParameter('style', 'line', @ischar); % line, stairs

            p.addParameter('alpha', 1, @isscalar); % alpha for main traces
            p.addParameter('showSem', false, @islogical); % show standard error traces?
            p.addParameter('errorAlpha', 0.5, @isscalar); % alpha for surrounding error fills

            % show less by default than in TDCA
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('markShowOnAxis', true, @islogical);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 2, @isscalar);
            p.addParameter('markShowInLegend', true, @islogical);

            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalShowOnAxis', true, @islogical);
            p.addParameter('intervalAlpha', 0.5, @isscalar);

            p.addParameter('showRangesOnData', false, @islogical); % show ranges for marks on traces
            p.addParameter('showRangesOnAxis', true, @islogical); % show ranges for marks below axis

            p.addParameter('timeAxisStyle', 'tickBridge', @ischar); % 'tickBridge' or 'marker'

            % make room for labels using AutoAxis
            p.addParameter('axisMarginLeft', 2.5, @isscalar);

            p.parse(varargin{:});

            indTrial = p.Results.individualTrials;
            if indTrial
                pset.assertSimultaneous();
            end

            axh = newplot;
            wasHolding = ishold(axh);
            hold(axh, 'on');

            if isempty(p.Results.basisIdx)
                basisIdx = 1:min(pset.nBases, 20);
            else
                basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            end
            nBasesPlot = numel(basisIdx);

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

            % if isempty(pset.interAlignGaps)
                % compute absolute x-gap between alignments
                alignGaps = repmat(alignGapFraction*sum(timeWidthByAlign) / (1 - alignGapFraction*nAlignUsed), nAlignUsed-1, 1);
            % else
            %     alignGaps = pset.interAlignGaps;
            % end

            % keep track of start and stop of each align in time
            tAlignZero = nanvec(nAlignUsed);

            % concatenate the data across alignments in order to compute
            % normalization constants. data is N x CTA

            if ~indTrial
                if ~isempty(p.Results.dataRandomIndex)
                    % use a sample from dataMeanRandomized instead
                    [data, indexInfo] = pset.arrangeCTAbyN('type', 'meanRandom', 'validBasesOnly', true, ...
                        'dataRandomIndex', p.Results.dataRandomIndex, ...
                        'basisIdx', basisIdx, 'conditionIdx', p.Results.conditionIdx, 'alignIdx', alignIdx);
                else
                    [data, indexInfo] = pset.arrangeCTAbyN('validBasesOnly', true, 'basisIdx', basisIdx, ...
                        'conditionIdx', p.Results.conditionIdx, 'alignIdx', alignIdx);
                end
                data = data';

                if p.Results.showSem
                    if ~isempty(p.Results.dataRandomIndex)
                        % use a sample from dataSemRandomized instead
                        dataSem = pset.arrangeCTAbyN('type', 'semRandom', 'validBasesOnly', true, ...
                            'dataRandomIndex', p.Results.dataRandomIndex, ...
                            'basisIdx', basisIdx, 'conditionIdx', p.Results.conditionIdx, 'alignIdx', alignIdx);
                    else
                        dataSem = pset.arrangeCTAbyN('type', 'sem', 'validBasesOnly', true, 'basisIdx', basisIdx, ...
                            'conditionIdx', p.Results.conditionIdx, 'alignIdx', alignIdx);
                    end
                    dataSem = dataSem';
                end
            else
                [data, indexInfo] = pset.simultaneous_arrangeCTAbyN('basisIdx', basisIdx, 'validBasesOnly', p.Results.validBasesOnly, 'conditionIdx', p.Results.conditionIdx, 'alignIdx', alignIdx, 'trialIdx', trialIdx, 'validTrialsOnly', true);
                dataSem = zeros(size(data)); % to ease normalization calculations below only
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

            appearances = pset.conditionDescriptor.appearances;

            nConditionsPlot = numel(indexInfo.condition);
            conditionIdx = indexInfo.condition;
            conditionNames = pset.conditionNames(conditionIdx);
            alignIdx = indexInfo.align;
            basisIdx = indexInfo.basis; nBasesPlot = numel(basisIdx);

            hData = cell(nConditionsPlot, nAlignUsed);
            for iAlign = 1:nAlignUsed
                idxAlign = indexInfo.align(iAlign);
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
                if ~indTrial
                    if isempty(p.Results.dataRandomIndex)
                        data = pset.dataMean{idxAlign}(indexInfo.basis, indexInfo.condition, :);
                    else
                        data = pset.dataMeanRandomized{idxAlign}(indexInfo.basis, indexInfo.condition, :, p.Results.dataRandomIndex);
                    end
                else
                    error('not yet implemented');
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

                stairsXOffset = pset.spikeFilter.binAlignmentMode.getBinStartOffsetForBinWidth(pset.timeDelta); % for stairs plotting only

                % draw error bars
                if p.Results.showSem
                    for iCondition = 1:nConditionsPlot
                        idxCondition = conditionIdx(iCondition);
                        app = appearances(idxCondition);
                        dataC = TensorUtils.squeezeDims(data(:, iCondition, :), 2);
                        dataSemC = TensorUtils.squeezeDims(dataSem(:, iCondition, :), 2);
                        for iBasis = 1:nBasesPlot
                            if strcmp(p.Results.style, 'stairs')
                                % offset the plot so as to resemble the binning
                                % mode used
                                [~, hShade] = TrialDataUtilities.Plotting.stairsError(...
                                     tvecPlot + stairsXOffset, dataC(iBasis, :), dataSemC(iBasis, :), ...
                                     'axh', axh, 'errorAlpha', p.Results.errorAlpha, 'color', app(idxCondition).Color, 'alpha', p.Results.alpha, ...
                                     'errorStyle', 'fill', 'errorColor', app.Color, 'lineWidth', p.Results.lineWidth * app.LineWidth, 'showLine', false);
                            else
                                hShade = TrialDataUtilities.Plotting.errorshade(tvecPlot, dataC(iBasis, :), ...
                                    dataSemC(iBasis, :), app.Color, 'axh', axh, ...
                                    'alpha', p.Results.errorAlpha, 'z', 0, 'showLine', false);
                            end
                            TrialDataUtilities.Plotting.hideInLegend(hShade);
                        end
                    end
                end

                % draw means
                for iCond = 1:nConditionsPlot
                    idxCondition = conditionIdx(iCond);
                    dataC = TensorUtils.squeezeDims(data(:, iCond, :), 2);
                    app = appearances(idxCondition);
                    plotArgsC = app.getPlotArgsCombinedWithDefaults('LineWidth', p.Results.lineWidth, 'Color', [0 0 0], 'Alpha', p.Results.alpha);
                    if strcmp(p.Results.style, 'stairs')
                        hData{iCond, iAlign} = TrialDataUtilities.Plotting.stairs(...
                                 tvecPlot + stairsXOffset, dataC', ...
                                 'axh', axh, plotArgsC{:});
                    else
                        hData{iCond, iAlign} = plot(tvecPlot, dataC', '-', ...
                            'Parent', axh, plotArgsC{:});
                    end

                    if iAlign==1
                        TrialDataUtilities.Plotting.showFirstInLegend(hData{iCond, iAlign}, conditionNames{iCond});
                    else
                        TrialDataUtilities.Plotting.hideInLegend(hData{iCond, iAlign});
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
                            'showInLegend', p.Results.markShowInLegend, ...
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
            au.yUnits = pset.dataUnitsCommon;

            % add scale bars to right side of axis
            au.yUnits = pset.dataUnitsCommon;
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
                    label = sprintf('%g %s', actualValue, pset.dataUnitsCommon);
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
                            label = sprintf('%g %s', actualValue, pset.dataUnitsCommon);
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

            pset.setupTimeAxisMultiAlign('axh', axh, 'doUpdate', true, ...
                'showMarks', p.Results.markShowOnAxis, ...
                'showIntervals', p.Results.intervalShowOnAxis, ...
                'showRanges', p.Results.showRangesOnAxis, ...
                'timeAxisStyle', p.Results.timeAxisStyle, ...
                'tOffsetByAlign', tAlignZero, ...
                'alignIdx', alignIdx, 'basisIdx', basisIdx);

            hold(axh, 'off');

            set(axh, 'SortMethod', 'childorder');
        end

        function setupTimeAxisMultiAlign(pset, varargin)
            p = inputParser();
            p.addParameter('tOffsetByAlign', pset.offsetsTimeDataMeanForPlotting, @isvector);
            p.addParameter('axh', gca, @ishandle);
            p.addParameter('doUpdate', true, @islogical);
            p.addParameter('showMarks', true, @islogical);
            p.addParameter('showIntervals', true, @islogical);
            p.addParameter('showRanges', true, @islogical); % show ranges for marks below axis
            p.addParameter('timeAxisStyle', 'tickBridge', @ischar); % 'tickBridge' or 'marker'
            p.addParameter('alignIdx', 1:pset.nAlign, @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.parse(varargin{:});

            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);

            % if we're only using one basis, use that basis' align summary
            % otherwise just use the aggregate over all bases to save time
            if numel(basisIdx) == 1 && ~isempty(pset.basisDataSourceIdx)
                asSet = pset.alignSummaryData(pset.basisDataSourceIdx(basisIdx), alignIdx);
            else
                asSet = pset.alignSummaryAggregated(alignIdx);
            end

            AlignSummary.setupTimeAutoAxisForMultipleAligns(asSet, pset.tvecDataMean, ...
                'tOffsetByAlign', p.Results.tOffsetByAlign, ...
                'axh', p.Results.axh, ...
                'doUpdate', p.Results.doUpdate, ...
                'showMarks', p.Results.showMarks, ...
                'showIntervals', p.Results.showIntervals, ...
                'showRanges', p.Results.showRanges, ...
                'timeAxisStyle', p.Results.timeAxisStyle, ...
                'tUnits', pset.timeUnitName);
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
            p.addParameter('timeDelta', pset.timeDelta, @isscalar); % allow for data to be resampled in time

            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('lineWidth', 1, @isscalar);

            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markOutline', true, @islogical);
            p.addParameter('markOutlineAlpha', 1, @isscalar);
            p.addParameter('markSize', 10, @isscalar);
            p.addParameter('useThreeVector', false, @islogical);
            p.addParameter('threeVectorLength', 1, @isscalar);
            p.addParameter('useTranslucentMark3d', true, @islogical);
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('zeroShowOnData', false, @islogical);
            p.addParameter('startShowOnData', false, @islogical);
            p.addParameter('stopShowOnData', false, @islogical);
            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);
            p.addParameter('clipping', 'off', @ischar);

            p.addParameter('showRangesOnData', true, @islogical); % show ranges for marks on traces

            p.addParameter('spliceAlignments', false, @islogical); % use sppline interpolation to splice
            p.addParameter('spliceOptions', struct(), @isstruct);

            p.addParameter('plotTubes', false, @islogical);
            p.addParameter('tubeRadius', 1, @isscalar);
            p.addParameter('tubeOptions', struct(), @isstruct);

            p.CaseSensitive = false;
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

            if p.Results.spliceAlignments
                data = pset.arrangeNbyCbyTA('spliceAlignments', true, ...
                    'alignIdx', alignIdx, 'basisIdx', basisIdx, ...
                    p.Results.spliceOptions, 'timeDelta', p.Results.timeDelta);

                for iCondition = 1:nConditions
                    c = conditionIdx(iCondition);
                    app = pset.conditionDescriptor.appearances(c);
%                     plotArgsC = appear.getPlotArgs();
                    plotArgsC = {};

                    dataMat = squeeze(data(basisIdx, c, :));

                    if p.Results.plotTubes
                        h = TrialDataUtilities.Plotting.tubeplot(dataMat(1:3, :), p.Results.tubeRadius, 'FaceAlpha', p.Results.alpha, ...
                            'EdgeColor', 'none', 'FaceColor', appear.Color, 'Clipping', p.Results.clipping);
                    else
                        if use3d
                            h = plot3(axh, dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                '-', 'LineWidth', p.Results.lineWidth * app.LineWidth, 'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        else
                            h = plot(axh, dataMat(1, :), dataMat(2, :), '-', ...
                                'LineWidth', p.Results.lineWidth * app.LineWidth, ...
                                'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        end
                        if p.Results.alpha < 1
                            TrialDataUtilities.Plotting.setLineOpacity(h, p.Results.alpha);
                        end
                    end
                    TrialDataUtilities.Plotting.showFirstInLegend(h, pset.conditionDescriptor.namesShort{c});
                end
            else
                [dataMean, indexInfo, tvecCell] = pset.arrangeNbyCbyTA('timeDelta', p.Results.timeDelta, 'splitByAlign', true, ...
                    'basisIdx', basisIdx, 'conditionIdx', conditionIdx, 'alignIdx', alignIdx);
                %                 dataMean = squeeze(TensorUtils.splitAlongDimension(dataMean, 3, pset.nTimeDataMean));

                for iAlign = 1:nAlignUsed
                    idxAlign = alignIdx(iAlign);
                    %                     tvec = pset.tvecDataMean{idxAlign};
                    data = dataMean{idxAlign};
                    for iCondition = 1:nConditions
                        idxCondition = conditionIdx(iCondition);
                        app = pset.conditionDescriptor.appearances(idxCondition);
                        plotArgsC = app.getPlotArgsCombinedWithDefaults('Color', [0 0 0], 'LineWidth', p.Results.lineWidth, 'Alpha', p.Results.alpha);
                        dataMat = squeeze(data(:, iCondition, :));

                        if use3d
                            h = plot3(axh, dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                '-', plotArgsC{:}, 'Clipping', p.Results.clipping, ...
                                plotArgs{:});
                        else
                            %                             dataMat = [dataVec1 dataVec2];
                            h = plot(axh, dataMat(1, :), dataMat(2, :), '-', ...
                                plotArgsC{:}, ...
                                'Clipping', p.Results.clipping, ...
                                plotArgs{:});
                        end
%                         if p.Results.alpha < 1
%                             TrialDataUtilities.Plotting.setLineOpacity(h, p.Results.alpha);
%                         end

                        if iAlign == 1
                            TrialDataUtilities.Plotting.showFirstInLegend(h, pset.conditionDescriptor.namesShort{idxCondition});
                        else
                            TrialDataUtilities.Plotting.hideInLegend(h);
                        end
                    end
                end
            end

            % plot marks and intervals on each condition
            for iAlign = 1:nAlignUsed
                idxAlign = alignIdx(iAlign);
                tvec = tvecCell{iAlign};
                data = dataMean{iAlign};
                % draw marks and intervals on the data traces
                as = pset.alignSummaryAggregated{idxAlign};
                % data is nBases x C x T; drawOnDataByConditions needs T x nBasesPlot x C
                dataForDraw = permute(data(:, conditionIdx, :), [3 1 2]);
                as.drawOnDataByCondition(tvec, dataForDraw, ...
                    'conditionIdx', conditionIdx, 'markAlpha', p.Results.markAlpha, ...
                    'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
                    'showZero', p.Results.zeroShowOnData, ...
                    'showStart', p.Results.startShowOnData, ...
                    'showStop', p.Results.stopShowOnData, ...
                    'markAlpha', p.Results.markAlpha, ...
                    'markOutline', p.Results.markOutline, ...
                    'markOutlineAlpha', p.Results.markOutlineAlpha, ...
                    'intervalAlpha', p.Results.intervalAlpha, ...
                    'showRanges', p.Results.showRangesOnData, ...
                    'markSize', p.Results.markSize, 'clipping', p.Results.clipping);
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
            else
                axis(axh, 'off');
            end

            if p.Results.useThreeVector
                tv = ThreeVector(axh);
                tv.vectorLength = p.Results.threeVectorLength;
            end

            hold(axh, 'off');
            set(get(axh, 'Parent'), 'CurrentAxes', axh);
            axis(axh, 'vis3d');
        end

        function plotStateSpaceIndividualTrials(pset, varargin)
            % plot a 2d or 3d basis1 x basis2 x basis3 trajectory plot
            p = inputParser;
            p.addParameter('basisIdx', 1:min(pset.nBases, 3), @(x) isvector(x) && ...
                all(inRange(x, [1 pset.nBases])));
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', true(pset.nAlign), @isnumeric);
            p.addParameter('trialIdx', truevec(pset.nTrials), @isvector);
            p.addParameter('timeDelta', pset.timeDelta, @isscalar); % allow for data to be resampled in time

            p.addParameter('maxTrials', Inf, @isscalar);

            p.addParameter('plotOptions', {}, @(x) iscell(x));
            p.addParameter('lineWidth', 2, @isscalar);
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('markAlpha', 1, @isscalar);
            p.addParameter('markSize', 10, @isscalar);
            p.addParameter('useThreeVector', true, @islogical);
            p.addParameter('threeVectorLength', 1, @isscalar);
            p.addParameter('useTranslucentMark3d', true, @islogical);
            p.addParameter('markShowOnData', false, @islogical);
            p.addParameter('intervalShowOnData', false, @islogical);
            p.addParameter('intervalAlpha', 1, @isscalar);
            p.addParameter('clipping', 'off', @ischar);

            p.addParameter('spliceAlignments', false, @islogical); % use sppline interpolation to splice
            p.addParameter('spliceOptions', struct(), @isstruct);

            p.addParameter('plotTubes', false, @islogical);
            p.addParameter('tubeRadius', 1, @isscalar);
            p.addParameter('tubeOptions', struct(), @isstruct);

            p.CaseSensitive = false;
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

            if p.Results.spliceAlignments
                [data, indexInfo, tvecCell]  = pset.simultaneous_arrangeNbyRbyTA('spliceAlignments', true, ...
                    'alignIdx', alignIdx, 'basisIdx', basisIdx, 'validBasesOnly', true, 'conditionIdx', conditionIdx, ...
                    'trialIdx', p.Results.trialIdx, 'validTrialsOnly', true, 'timeDelta', p.Results.timeDelta, ...
                    'maxTrials', p.Results.maxTrials);

                for iTrial = 1:numel(indexInfo.trial)
                    r = indexInfo.trial(iTrial);
                    c = pset.dataSources{1}.conditionIdx(r);
                    appear = pset.conditionDescriptor.appearances(c);
                    plotArgsC = appear.getPlotArgs();

                    dataMat = squeeze(data(:, iTrial, :));

                    if p.Results.plotTubes
                        h = TrialDataUtilities.Plotting.tubeplot(dataMat(1:3, :), p.Results.tubeRadius, 'FaceAlpha', p.Results.alpha, ...
                            'EdgeColor', 'none', 'FaceColor', appear.Color, 'Clipping', p.Results.clipping);
                    else
                        if use3d
                            h = plot3(axh, dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                '-', 'LineWidth', p.Results.lineWidth * app(idxCondition).LineWidth, 'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        else
                            h = plot(axh, dataMat(1, :), dataMat(2, :), '-', ...
                                'LineWidth', p.Results.lineWidth * app(idxCondition).LineWidth, ...
                                'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        end
                        if p.Results.alpha < 1
                            TrialDataUtilities.Plotting.setLineOpacity(h, p.Results.alpha);
                        end
                    end
                    TrialDataUtilities.Plotting.showFirstInLegend(h, pset.conditionDescriptor.namesShort{c});
                end
            else
                [dataByAlign, indexInfo, tvecCell]  = pset.arrangeNbyRbyTA('splitByAlign', true, ...
                    'alignIdx', alignIdx, 'basisIdx', basisIdx, 'validBasesOnly', true, 'conditionIdx', conditionIdx, ...
                    'trialIdx', p.Results.trialIdx, 'validTrialsOnly', true, 'timeDelta', p.Results.timeDelta, ...
                    'maxTrials', p.Results.maxTrials);


                for iAlign = 1:nAlignUsed

                    data = dataByAlign{iAlign};
                    for iTrial = 1:numel(indexInfo.trial)
                        r = indexInfo.trial(iTrial);
                        c = pset.dataSources{1}.conditionIdx(r);
                        appear = pset.conditionDescriptor.appearances(c);
                        plotArgsC = appear.getPlotArgs();
                        dataMat = squeeze(data(:, iTrial, :));

                        if use3d
                            h = plot3(axh, dataMat(1, :), dataMat(2, :), dataMat(3, :), ...
                                '-', 'LineWidth', p.Results.lineWidth * app(idxCondition).LineWidth, 'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        else
                            %                             dataMat = [dataVec1 dataVec2];
                            h = plot(axh, dataMat(1, :), dataMat(2, :), '-', ...
                                'LineWidth', p.Results.lineWidth * app(idxCondition).LineWidth, ...
                                'Clipping', p.Results.clipping, ...
                                plotArgsC{:}, plotArgs{:});
                        end
                        if p.Results.alpha < 1
                            TrialDataUtilities.Plotting.setLineOpacity(h, p.Results.alpha);
                        end

                        if iAlign == 1
                            TrialDataUtilities.Plotting.showFirstInLegend(h, pset.conditionDescriptor.namesShort{c});
                        else
                            TrialDataUtilities.Plotting.hideInLegend(h);
                        end
                    end
                end

%                 % plot marks and intervals on each condition
%                 for iAlign = 1:nAlignUsed
%                     idxAlign = alignIdx(iAlign);
%                     tvec = tvecCell{iAlign};
%                     data = dataMean{iAlign};
%
%                     for iCondition = 1:nConditions
%                         c = conditionIdx(iCondition);
%
%                         % draw marks and intervals on the data traces
%                         as = pset.alignSummaryAggregated{idxAlign};
%                         % data is nBases x C x T; drawOnDataByConditions needs T x nBasesPlot x C
%                         dataForDraw = permute(data(basisIdx, c, :), [3 1 2]);
%                         as.drawOnDataByCondition(tvec, dataForDraw, ...
%                             'conditionIdx', c, 'markAlpha', p.Results.markAlpha, ...
%                             'showMarks', p.Results.markShowOnData, 'showIntervals', p.Results.intervalShowOnData, ...
%                             'useTranslucentMark3d', p.Results.useTranslucentMark3d, ...
%                             'alpha', p.Results.alpha, ...
%                             'intervalAlpha', p.Results.intervalAlpha, ...
%                             'showRanges', p.Results.showRangesOnData, ...
%                             'markSize', p.Results.markSize, 'clipping', p.Results.clipping);
%                     end
%                 end
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
            else
                axis(axh, 'off');
            end

            if p.Results.useThreeVector
                tv = ThreeVector(axh);
                tv.vectorLength = p.Results.threeVectorLength;
            end

            hold(axh, 'off');
            set(get(axh, 'Parent'), 'CurrentAxes', axh);
            axis(axh, 'vis3d');
        end
    end

    methods % Arrange data matrices trial averaged
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

        function [NbyCbyTA, indexInfo, tvec, avec] = arrangeNbyCbyTA(pset, varargin)
            % [NbyCbyTA, tvec, avec] = arrangeNbyCbyTA(pset, ...)
            p = inputParser();
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('validTimepointsAllConditionsBasesOnly', false, @islogical); % keep only timepoints where all bases have data
            p.addParameter('type', 'mean', @ischar); % mean, sem, meanRandom
%             p.addParameter('dataRandomIndex', 1:pset.nRandomSamples, @isvector);
            p.addParameter('spliceAlignments', false, @islogical); % use sppline interpolation to splice
            p.addParameter('spliceOptions', struct(), @isstruct);
            p.addParameter('timeDelta', pset.timeDelta, @isscalar); % allow for data to be resampled in time
            p.addParameter('splitByAlign', false, @islogical);
            p.addParameter('offsetTimeForPlotting', false, @islogical);

            p.parse(varargin{:});
            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            nAlign = numel(alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            conditionIdx = TensorUtils.vectorMaskToIndices(p.Results.conditionIdx);
%             dataRandomIndex = TensorUtils.vectorMaskToIndices(p.Results.dataRandomIndex);
%             checkDataRandomIndex = @() assert(all(dataRandomIndex >= 1 & dataRandomIndex <= pset.nRandomSamples), ...
%                 'dataRandomIndex must be in range 1:%d (nRandomSamples)', pset.nRandomSamples);
            if p.Results.validBasesOnly
                mask = pset.basisValid(basisIdx);
                basisIdx = basisIdx(mask);
            end

            if strcmp(p.Results.type, 'std')
                dataStd = pset.computeDataStd();
            elseif strcmp(p.Results.type, 'stdRandom')
                dataStd = pset.computeDataStdRandomized('dataRandomIndex', dataRandomIndex);
            end

            if strcmp(p.Results.type, 'var')
                dataVar = pset.computeDataVar();
            elseif strcmp(p.Results.type, 'varRandom')
                dataVar = pset.computeDataVarRandomized('dataRandomIndex', dataRandomIndex);
            end

            data = cellvec(nAlign);
            for iAlign = 1:numel(alignIdx)
                idxAlign = alignIdx(iAlign);
                switch p.Results.type
                    case 'mean'
                        data{iAlign} = pset.dataMean{idxAlign}(basisIdx, conditionIdx, :);
                    case 'sem'
                        data{iAlign} = pset.dataSem{idxAlign}(basisIdx, conditionIdx, :);
                    case 'std'
                        data{iAlign} = dataStd{idxAlign}(basisIdx, conditionIdx, :);
                    case 'var'
                        data{iAlign} = dataVar{idxAlign}(basisIdx, conditionIdx, :);

                    case 'meanRandom'
                        assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');
                        checkDataRandomIndex();
                        data{iAlign} = pset.dataMeanRandomized{idxAlign}(basisIdx, conditionIdx, :, dataRandomIndex);
                    case 'semRandom'
                        assert(pset.hasDataRandomized, 'Must generate randomized data using storeDataMean* method first');
                        checkDataRandomIndex();
                        data{iAlign} = pset.dataSemRandomized{idxAlign}(basisIdx, conditionIdx, :, p.Results.dataRandomIndex);
                    case 'stdRandom'
                        data{iAlign} = dataStd{idxAlign}(basisIdx, conditionIdx, :, :); % dataRandomIndex handled above
                    case 'varRandom'
                        data{iAlign} = dataVar{idxAlign}(basisIdx, conditionIdx, :, :); % dataRandomIndex handled above

                    otherwise
                        error('Unknown data type %s', p.Results.type);
                end
            end

            basisValidMask = pset.basisValid(basisIdx);
            tvecCell = pset.tvecDataMean(alignIdx);

            if p.Results.spliceAlignments
                % data is N x C x T --> N x T x C
                data = cellfun(@(d) permute(d, [1 3 2]), data, 'UniformOutput', false);

                data = TrialDataUtilities.Data.spliceTrajectories(data, 'basisMask', basisValidMask, p.Results.spliceOptions);
                data = ipermute(data, [1 3 2]);

                % reslice into cell
                data = squeeze(TensorUtils.splitAlongDimension(data, 3, pset.nTimeDataMean));
            end

            if p.Results.timeDelta ~= pset.timeDelta
                % interp to new time delta
                for iA = 1:nAlign
                    [data{iA}, tvecCell{iA}] = TrialDataUtilities.Data.resampleTensorInTime(data{iA}, 3 ,...
                        tvecCell{iA}, 'timeDelta', p.Results.timeDelta);
                end
            end

            if p.Results.offsetTimeForPlotting
                offsets = pset.getAlignPlottingTimeOffsets(tvecCell, 'alignIdx', indexInfo.align);
                for iAlign = 1:pset.nAlign
                    tvecCell{iAlign} = tvecCell{iAlign} + offsets(iAlign);
                end
            end

            if p.Results.splitByAlign
                for iA = 1:nAlign
                    data{iA}(~basisValidMask, :, :) = NaN;
                    if p.Results.validTimepointsAllConditionsBasesOnly
                        Tmask = all(all(~isnan(data{iA}(basisValidMask, :, :)), 1), 2);
                        data{iA} = data{iA}(:, :, Tmask);
                        tvecCell{iA} = tvecCell{iA}(Tmask);
                    end
                end
                NbyCbyTA = data;
                tvec = tvecCell;
                avec = [];
            else
                [NbyCbyTA, avecRaw] = TensorUtils.catWhich(3, data{:});
                avec = makecol(alignIdx(avecRaw));
                tvec = cat(1, tvecCell{:});
                % nan out invalid bases
                NbyCbyTA(~basisValidMask, :, :) = NaN;

                if p.Results.validTimepointsAllConditionsBasesOnly
                    TAmask = all(all(~isnan(NbyCbyTA(basisValidMask, :, :)), 1), 2);
                    NbyCbyTA = NbyCbyTA(:, :, TAmask);
                    tvec = tvec(TAmask);
                end
                % nan out invalid tiempoints
            end

            indexInfo.align = alignIdx;
            indexInfo.condition = conditionIdx;
            indexInfo.basis = basisIdx;
        end

        function [CTAbyN, indexInfo, cvec, tvec, avec, nvec] = arrangeCTAbyN(pset, varargin)
            % out is C*T*A x N concatenated matrix
            % see arrangeNbyCbyTA for additional options
            p = inputParser();
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('validTimepointsAllBasesOnly', false, @islogical); % keep only timepoints where all bases have data
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [NbyCbyTA, indexInfo, tvec, avec] = pset.arrangeNbyCbyTA('validBasesOnly', p.Results.validBasesOnly, p.Unmatched);
            labels = {indexInfo.basis, indexInfo.condition, [tvec, avec]};
            % we include dim 4 as this would hold randomized samples
            [CTAbyN, labelsOut] = TensorUtils.reshapeByConcatenatingDims(NbyCbyTA, {[3 2], 1}, labels);

            if p.Results.validTimepointsAllBasesOnly
                if p.Results.validBasesOnly
                    CTAmask = all(~isnan(CTAbyN), 2);
                else
                    % ignore the invalid bases - these will always be nan
                    % anyway
                    basisValid = pset.basisValid(indexInfo.basis);
                    CTAmask = all(~isnan(CTAbyN(:, basisValid)), 2);
                end
                CTAbyN = CTAbyN(CTAmask, :);
            end

            cvec = labelsOut{1}(:, 1);
            tvec = labelsOut{1}(:, 2);
            avec = labelsOut{1}(:, 3);
            nvec = labelsOut{2};
        end

        function [NbyTAbyAttr, indexInfo, tvec, avec] = arrangeNbyTAbyConditionAttributes(pset, varargin)
            % build a tensor of N by T by nValsAttr1 by nValsAttr2 x ...
            % this tensor is the format used by dpca_covs
            [NbyCbyTA, indexInfo, tvec, avec] = pset.arrangeNbyCbyTA(varargin{:});
            NbyTAbyC = permute(NbyCbyTA, [1 3 2]);
            N = size(NbyTAbyC, 1);
            TA = size(NbyTAbyC, 2);
            condSize = pset.conditionDescriptor.conditionsSize;
            NbyTAbyAttr = reshape(NbyTAbyC, [N TA makerow(condSize)]);
        end
    end

    methods % Individual trial data
        function assertSimultaneous(pset)
            assert(pset.simultaneous, 'This method may only be called on simultaneously collected data sets, i.e. where there is only one data source');
        end

        function [dataByTrialSubset, indexInfo] = computeDataByTrialSubset(pset, varargin)
            % this function can be used directly, but it is mostly used as
            % the workhorse for single trial data access methods below. It samples trials from
            % each basis' data source based on a specified sampling method.
            %
            % dataByTrialSubset is nBases x nAlign cell of numTrials (or shorter) x nTime.
            %
            % sampling mode parameter:
            %  - takeFirstValid [default] take the first numTrials valid trials.
            %      If a given basis has too few trials, the remainder will not be included
            %  - takeFirst take the first numTrials trials
            %      including valid and invalid trials. If a given basis has too few trials, the remainder will be NaN.
            %  - fixedOrder: take the trials specified by 'trialIdx'. If a
            %      given basis has too few trials, the remainder will not be included
            %  - 'randomWithoutReplacement': pick a random subset of
            %      trials without replacement. If a given basis has too few trials,
            %      the remainder will be NaN.
            %  - 'randomWithReplacement': pick a random subset of trials
            %      with replacement. All bases with at least 1 trial will have
            %      a full set of numTrials
            %
            %  'sameTrialsEachBasis' - defaults to pset.simultaneous. If true, the same
            %    trialIdx will be used across all bases. If false, different
            %    idx will be used for each basis. For simultaneous datasets,
            %    this shuffles the trial-to-trial correlation.
            %
            % 'fillNanToNumTrials' - pads each cell's contents to be size
            %    numTrials along dim 1 by padding with NaN. Wen combined with
            %   commonTime == true, every cell will be guaranteed to have the
            %   same size as needed for concatenation

            p = inputParser();
            p.addParameter('sameTrialsEachBasis', pset.simultaneous, @islogical);
            p.addParameter('samplingMode', 'takeFirstValid', @(x) any(validatestring(x, {'takeFirstValid', ...
                'takeFirst', 'fixedOrder', 'randomWithoutReplacement', 'randomWithReplacement'}, ...
                'computeDataByTrialSubset', 'samplingMode'))); % see above.
            p.addParameter('numTrials', NaN, @isscalar); % the exact number of trials to sample, must be specified manually if not simultaneous
            p.addParameter('truncateNumTrialsToMaxOverBases', false, @islogical); % if no basis has enough trials, lower numTrials to the max across bases
            p.addParameter('trialIdx', [], @(x) isempty(x) || isvector(x)); % when fixedOrder is the sampling mode
            p.addParameter('fillNanToNumTrials', false, @islogical); % ensure each cell has numTrials by filling with NaN if not enough trials are present

            p.addParameter('nRepeats', 1, @isscalar);
            p.addParameter('randomSeed', 0, @isscalar); % randomSeed will

            % give preference for trials that aren't missing too many samples
            % when using samplingMode takeFirstValid, randomWithoutReplacement, or randomWithReplacement
            p.addParameter('ignoreTrialsWithTooFewSamples', true, @islogical);
            p.addParameter('minFractionTrialSamples', 0.8, @isscalar);

            % subsetting the data considered
            p.addParameter('commonTime', true, @islogical); % use the dataMean time vector so that all bases have the same vectors
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector); % used to mask which conditions the trials will be drawn from
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical); % select only valid bases
            p.parse(varargin{:});

            sameTrialsEachBasis = p.Results.sameTrialsEachBasis;
            if sameTrialsEachBasis && ~pset.simultaneous
                 warning('sameTrialsEachBasis is true but Pset is not simultaneously recorded. Trials may not match.');
            end

            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            nAlign = numel(alignIdx);

            conditionIdx = TensorUtils.vectorMaskToIndices(p.Results.conditionIdx);

            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            if p.Results.validBasesOnly
                basisIdx = basisIdx(pset.basisValid(basisIdx));
            end
            nBases = numel(basisIdx);

            randomSeed = p.Results.randomSeed;
            trialIdxManual = TensorUtils.vectorMaskToIndices(p.Results.trialIdx);
            samplingMode = p.Results.samplingMode;
            nRepeats = p.Results.nRepeats;

            % check sampling mode
            if strcmp(samplingMode, 'fixedOrder')
                assert(~isempty(p.Results.trialIdx), 'Must specify trialIdx when using samplingMode fixedOrder');
            end
            if ~ismember(samplingMode, {'randomWithoutReplacement', 'randomWithReplacement'})
                assert(nRepeats == 1, 'nRepeats > 1 only makes sense with samplingMode randomWithoutReplacement or randomWithoutReplacement');
            end

            % decide on numTrials
            numTrials = p.Results.numTrials;
            if isnan(numTrials)
                switch samplingMode
                    case 'fixedOrder'
                        numTrials = numel(trialIdxManual);
                    case 'takeFirst'
                        assert(pset.simultaneous, 'For non-simultaneous datasets, numTrials must be specified with this samplingMode');
                        numTrials = pset.nTrials;
                    case 'takeFirstValid'
                        assert(pset.simultaneous, 'For non-simultaneous datasets, numTrials must be specified with this samplingMode');
                        numTrials = pset.nTrialsValid;
                    case {'randomWithoutReplacement', 'randomWithReplacement'}
                        error('numTrials must be specified with this samplingMode');
                end
            end

            % we sample from the full set of dataByTrial
            if p.Results.commonTime
                dataByTrial = pset.computeDataByTrialCommonTime();
            else
                dataByTrial = pset.dataByTrial;
            end

            % compute minimum samples required
            ignoreTrialsWithTooFewSamples = p.Results.ignoreTrialsWithTooFewSamples;
            minFractionTrialSamples = p.Results.minFractionTrialSamples;
            dataMean = pset.dataMean;

            function validTrialMask = getValidTrialMaskForBasis(iBasis)
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                validTrialMask = src.valid & ismember(src.conditionIdx, conditionIdx);

                % check for a sufficient set of samples that overlap with the corresponding dataMean
                if ignoreTrialsWithTooFewSamples
                    for iA = 1:nAlign
                        idxA = alignIdx(iA);

                        % nTrials x nTime matrix of where the trial's mean has a valid non-nan entry
                        % do this calculation within the set of validTrialMask to make it faster
                        validSamplesInMeanMatrix = shiftdim(~isnan(dataMean{idxA}(iBasis, src.conditionIdx(validTrialMask), :)), 1); % nnz(validTrialMask) x time
                        fracValidSamplesEachTrial = sum(~isnan(dataByTrial{iBasis, idxA}(validTrialMask, :)) .* validSamplesInMeanMatrix, 2) ./ sum(validSamplesInMeanMatrix, 2); % nnz(validTrialMask) x 1
                        validTrialMask(fracValidSamplesEachTrial < minFractionTrialSamples) = false; % fracValidSamplesEachTrial indexes into validTrialMask
                    end
                end
            end

            function [trialIdx, validTrialMask] = pickTrialsForBasis(iBasis, iRepeat, validTrialMask)
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};

                switch samplingMode
                    case 'fixedOrder'
                        trialIdx = trialIdxManual;

                    case 'takeFirst'
                        trialIdx = 1:min(src.nTrials, numTrials);

                    case 'takeFirstValid'
                        trialIdx = find(validTrialMask, numTrials); % will truncate automatically

                    case 'randomWithReplacement'
                        validTrialIdx = find(validTrialMask);
                        s = RandStream('mt19937ar', 'Seed', randomSeed + iRepeat + iBasis - 1);
                        trialIdx = TrialDataUtilities.Data.randsamplePopulation(s, validTrialIdx, numTrials, true); %#ok<FNDSB> % we have to do this indirect randsample

                    case 'randomWithoutReplacement'
                        validTrialIdx = find(validTrialMask);
                        s = RandStream('mt19937ar', 'Seed', randomSeed + iRepeat + iBasis - 1);
                        trialIdx = TrialDataUtilities.Data.randsamplePopulation(s, validTrialIdx, numTrials, false); %#ok<FNDSB> % we have to do this indirect randsample

                end
                trialIdx = makecol(trialIdx);
            end

            dataByTrialSubset = cell(nBases, nAlign, nRepeats);
            trialIdxInfo = cell(nBases, nRepeats);
            basisValid = pset.basisValid;

            % select common trialIdx up front across repeats
            if sameTrialsEachBasis
                trialIdxCommon = cellvec(nRepeats);
                firstValidBasis = find(basisValid, 1);
                validTrialMask = getValidTrialMaskForBasis(firstValidBasis);
                for iRepeat = 1:nRepeats
                    trialIdxCommon{iRepeat} = pickTrialsForBasis(firstValidBasis, iRepeat, validTrialMask);
                end
            end

            for iBasis = 1:nBases
                idxBasis = basisIdx(iBasis);

                validTrialMask = getValidTrialMaskForBasis(idxBasis);

                for iRepeat = 1:nRepeats
                    % trialIdxCommon is used when sameTrialsEachBasis == true
                    % trialIdx is the original picked trials
                    % trialIdxSelect are the valid trials we pick from dataByTrial

                    if ~basisValid(idxBasis)
                        trialIdxSelect = [];
                        trialIdx = [];

                    elseif ~sameTrialsEachBasis
                        % pick a new set of trials for this basis
                        trialIdx = pickTrialsForBasis(idxBasis, iRepeat, validTrialMask);
                        trialIdxSelect = trialIdx;

                    else
                        % use the trialIdxCommon
                        trialIdx = trialIdxCommon{iRepeat};
                        if ~pset.simultaneous
                            % we may only have some of the trials for this basis
                            % the rest will be NaN to maintain equivalency
                            validTrialMask = getValidTrialMaskForBasis(idxBasis);
                            trialIdxMask = trialIdxCommon{iRepeat} <= size(dataByTrial{idxBasis, alignIdx(1)}, 1);
                            trialIdxSelect = trialIdxCommon{iRepeat}(trialIdxMask);
                        else
                            trialIdxSelect = trialIdxCommon{iRepeat};
                        end
                    end

                    for iAlign = 1:nAlign
                        idxAlign = alignIdx(iAlign);
                        dataByTrialSubset{iBasis, iAlign, iRepeat} = dataByTrial{idxBasis, idxAlign}(trialIdxSelect, :);
                        markInvalid = ~validTrialMask(trialIdxSelect);
                        dataByTrialSubset{iBasis, iAlign, iRepeat}(markInvalid, :) = NaN;
                        if sameTrialsEachBasis && ~pset.simultaneous
                            % need to reinflate to keep the same size
                            dataByTrialSubset{iBasis, iAlign, iRepeat} = TensorUtils.inflateMaskedTensor(dataByTrialSubset{iBasis, iAlign}, 1, trialIdxMask);
                        end
                    end
                    trialIdxInfo{iBasis, iRepeat} = trialIdx;
                end
            end

            % compute the number of trials actually realized
            numTrialsActual = nan(size(trialIdxInfo, 1), 1);
            for i = 1:size(trialIdxInfo, 1)
                numTrialsActual(i) = numel(trialIdxInfo{i, 1});
            end

            if p.Results.truncateNumTrialsToMaxOverBases
                numTrials = min(numTrials, max(numTrialsActual(:)));
            end

            if p.Results.fillNanToNumTrials
                % expand each cell to numTrials with NaN. The original
                % number of trials can be found using
                % indexInfo.trial
                for i = 1:numel(dataByTrialSubset)
                    sz = size(dataByTrialSubset{i});
                    if sz(1) < numTrials
                        dataByTrialSubset{i} = cat(1, dataByTrialSubset{i}, nan(numTrials - sz(1), sz(2)));
                    end
                end
            end

            indexInfo.numTrials = numTrials;
            indexInfo.numTrialsActual = numTrialsActual;
            indexInfo.basis = basisIdx;
            indexInfo.basisValidMask = pset.basisValid(basisIdx);
            indexInfo.align = alignIdx;
            indexInfo.condition = conditionIdx;

            indexInfo.trial = trialIdxInfo;
            if sameTrialsEachBasis
                indexInfo.trialShared = trialIdxInfo{1, :}';
            else
                indexInfo.trialShared = [];
            end
        end

        function [dataByTrialSubset, indexInfo] = computeDataByTrialSubsetByCondition(pset, varargin)
            % this function can be used directly, but it is mostly used as
            % the workhorse for single trial data access methods below.
            % It is similar to computeDataByTrialSubset except numTrials trials are sampled per condition.
            %
            % dataByTrialSubset is nBases x nAlign x nConditions cell of numTrials (or shorter) x nTime.
            %
            % sampling mode parameter:
            %  - takeFirstValid [default] take the first numTrials valid trials.
            %      If a given basis has too few trials, the remainder will not be included
            %  - 'randomWithoutReplacement': pick a random subset of
            %      trials without replacement. If a given basis has too few trials,
            %      the remainder will be NaN.
            %  - 'randomWithReplacement': pick a random subset of trials
            %      with replacement. All bases with at least 1 trial will have
            %      a full set of numTrials
            %
            %  'sameTrialsEachBasis' - defaults to pset.simultaneous. If true, the same
            %    trialIdx will be used across all bases. If false, different
            %    idx will be used for each basis. For simultaneous datasets,
            %    this shuffles the trial-to-trial correlation.
            %
            % 'fillNanToNumTrials' - pads each cell's contents to be size
            %    numTrials along dim 1 by padding with NaN. Wen combined with
            %   commonTime == true, every cell will be guaranteed to have the
            %   same size as needed for concatenation
            %
            % if trialLists is specified, this will overwrite the value of pset.trialLists when assigning trials into

            p = inputParser();
            p.addParameter('sameTrialsEachBasis', pset.simultaneous, @islogical);
            p.addParameter('samplingMode', 'takeFirstValid', @(x) any(validatestring(x, ...
                {'takeFirstValid', 'randomWithoutReplacement', 'randomWithReplacement'}, ...
                'computeDataByTrialSubsetByCondition', 'samplingMode'))); % see above.
            p.addParameter('numTrials', NaN, @isscalar); % the exact number of trials to sample, must be specified manually if not simultaneous
            p.addParameter('truncateNumTrialsToMaxOverBases', false, @islogical); % if no basis has enough trials, lower numTrials to the max across bases
            p.addParameter('fillNanToNumTrials', false, @islogical); % ensure each cell has numTrials by filling with NaN if not enough trials are present

            p.addParameter('nRepeats', 1, @isscalar);
            p.addParameter('randomSeed', 0, @isscalar);

            % if specified, this will be used instead of pset.trialLists. Its dimensions
            % are nBases x nConditions x (1 or nRepeats). It is most useful when
            % optionally including a third dimension with length nRepeats, so
            % that slice of the trial lists will be used for the ith sample
            p.addParameter('trialLists', [], @iscell);

            % give preference for trials that aren't missing too many samples
            % when using samplingMode takeFirstValid, randomWithoutReplacement, or randomWithReplacement
            p.addParameter('ignoreTrialsWithTooFewSamples', true, @islogical);
            p.addParameter('minFractionTrialSamples', 0.8, @isscalar);

            % subsetting the data considered
            p.addParameter('commonTime', true, @islogical); % use the dataMean time vector so that all bases have the same vectors
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector); % used to mask which conditions the trials will be drawn from
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical); % select only valid bases
            p.parse(varargin{:});

            sameTrialsEachBasis = p.Results.sameTrialsEachBasis;
            if sameTrialsEachBasis && ~pset.simultaneous
                 error('sameTrialsEachBasis is not supported for byCondition when pset is not simultaneously recorded');
            end

            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            nAlign = numel(alignIdx);

            conditionIdx = TensorUtils.vectorMaskToIndices(p.Results.conditionIdx);
            nConditions = numel(conditionIdx);

            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            if p.Results.validBasesOnly
                basisIdx = basisIdx(pset.basisValid(basisIdx));
            end
            nBases = numel(basisIdx);

            randomSeed = p.Results.randomSeed;
            nRepeats = p.Results.nRepeats;
            samplingMode = p.Results.samplingMode;

            if isempty(p.Results.trialLists)
                trialLists = pset.trialLists;
            else
                trialLists = p.Results.trialLists;
                assert(size(trialLists, 1) == pset.nBases, 'triallists must have size nBases along dim 1');
                assert(size(trialLists, 2) == pset.nConditions, 'triallists must have size nConditions along dim 2');
                assert(ismember(size(trialLists, 3), [1 nRepeats]), 'triallists must have size 1 or nRepeats along dim 3');
            end
            % keep trialLists at nBases x nConditions x (1 or nRepeats)

            % check settings
            if ~ismember(samplingMode, {'randomWithoutReplacement', 'randomWithReplacement'})
                assert(nRepeats == 1, 'nRepeats > 1 only makes sense with samplingMode randomWithoutReplacement or randomWithoutReplacement');
            end

            % decide on numTrials
            numTrials = p.Results.numTrials;
            if isnan(numTrials)
                % default num trials is max over nTrials by
                % condition
                switch samplingMode
                    case 'takeFirstValid'
                        assert(pset.simultaneous, 'For non-simultaneous datasets, numTrials must be specified with this samplingMode');
                        numTrials = max(pset.nTrialsByCondition(:));
                    case {'randomWithoutReplacement', 'randomWithReplacement'}
                        error('numTrials must be specified with this samplingMode');
                    otherwise
                        error('samplingMode not supported for byCondition');
                end
            end

            % we sample from the full set of dataByTrial
            if p.Results.commonTime
                dataByTrial = pset.computeDataByTrialCommonTime();
            else
                dataByTrial = pset.dataByTrial;
            end

            % compute minimum samples required
            ignoreTrialsWithTooFewSamples = p.Results.ignoreTrialsWithTooFewSamples;
            minFractionTrialSamples = p.Results.minFractionTrialSamples;
            dataMean = pset.dataMean;

            function validTrialMask = getValidTrialMaskForBasis(iBasis)
                src = pset.dataSources{pset.basisDataSourceIdx(iBasis)};
                validTrialMask = src.valid & ismember(src.conditionIdx, conditionIdx);

                % check for a sufficient set of samples that overlap with the corresponding dataMean
                if ignoreTrialsWithTooFewSamples
                    for iA = 1:nAlign
                        idxA = alignIdx(iA);

                        % nTrials x nTime matrix of where the trial's mean has a valid non-nan entry
                        % do this calculation within the set of validTrialMask to make it faster
                        validSamplesInMeanMatrix = shiftdim(~isnan(dataMean{idxA}(iBasis, src.conditionIdx(validTrialMask), :)), 1); % nnz(validTrialMask) x time
                        fracValidSamplesEachTrial = sum(~isnan(dataByTrial{iBasis, idxA}(validTrialMask, :)) .* validSamplesInMeanMatrix, 2) ./ sum(validSamplesInMeanMatrix, 2); % nnz(validTrialMask) x 1
                        validTrialMask(fracValidSamplesEachTrial < minFractionTrialSamples) = false; % fracValidSamplesEachTrial indexes into validTrialMask
                    end
                end
            end

            function trialIdx = pickTrialsForBasisByCondition(iBasis, iRepeat)
                if size(trialLists, 3) == 1
                    listByCondition = trialLists(iBasis, conditionIdx)';
                else
                    listByCondition = trialLists(iBasis, conditionIdx, iRepeat)';
                end

                trialIdx = cell(nConditions, 1);
                switch samplingMode
                    case 'takeFirstValid'
                        for iC = 1:nConditions
                            if numel(listByCondition{conditionIdx(iC)}) > numTrials
                                trialIdx{iC} = listByCondition{conditionIdx(iC)}(1:numTrials);
                            else
                                trialIdx{iC} = listByCondition{conditionIdx(iC)};
                            end
                        end

                    case 'randomWithReplacement'
                        s = RandStream('mt19937ar', 'Seed', randomSeed + iRepeat + iBasis - 1);
                        for iC = 1:nConditions
                            trialIdx{iC} = TrialDataUtilities.Data.randsamplePopulation(s, listByCondition{conditionIdx(iC)}, numTrials, true); % we have to do this indirect randsample
                        end

                    case 'randomWithoutReplacement'
                        s = RandStream('mt19937ar', 'Seed', randomSeed + iRepeat + iBasis - 1);
                        for iC = 1:nConditions
                            trialIdx{iC} = TrialDataUtilities.Data.randsamplePopulation(s, listByCondition{conditionIdx(iC)}, numTrials, false); % we have to do this indirect randsample
                        end

                    otherwise
                        error('samplingMode not supported for byCondition');
                end
            end

            % byCondition
            dataByTrialSubset = cell(nBases, nAlign, nConditions, nRepeats);
            trialIdxInfo = cell(nBases, nConditions, nRepeats);
            basisValid = pset.basisValid;

            if sameTrialsEachBasis
                % select all idx at the beginning for all repeats to save time
                assert(pset.simultaneous); % should have already been checked above
                firstValidBasis = find(basisValid, 1);

                % subset valid trials from trialLists for first basis
                validTrialMask = getValidTrialMaskForBasis(firstValidBasis);
                trialLists(firstValidBasis, :, :) = cellfun(@(list) list(validTrialMask(list)), trialLists(firstValidBasis, :, :), 'UniformOutput', false);

                trialIdxCommon = cell(nRepeats, 1);
                for iRepeat = 1:nRepeats
                    trialIdxCommon{iRepeat} = pickTrialsForBasisByCondition(firstValidBasis);
                end

            else
                % subset valid trials from trialLists all up front
                for iBasis = 1:nBases
                    idxBasis = basisIdx(iBasis);
                    if ~basisValid(idxBasis), continue, end
                    validTrialMask = getValidTrialMaskForBasis(idxBasis);
                    trialLists(iBasis, :, :) = cellfun(@(list) list(validTrialMask(list)), trialLists(idxBasis, :, :), 'UniformOutput', false);
                end
            end

            for iBasis = 1:nBases
                idxBasis = basisIdx(iBasis);
                for iRepeat = 1:nRepeats
                    if ~basisValid(idxBasis)
                        trialIdx = cell(nConditions, 1);
                    elseif ~sameTrialsEachBasis
                        % pick a new set of trials each basis
                        trialIdx = pickTrialsForBasisByCondition(idxBasis, iRepeat);
                    else
                        trialIdx = trialIdxCommon{iRepeat};
                    end

                    for iAlign = 1:nAlign
                        idxAlign = alignIdx(iAlign);
                        for iCondition = 1:nConditions
                            dataByTrialSubset{iBasis, iAlign, iCondition, iRepeat} = dataByTrial{idxBasis, idxAlign}(trialIdx{iCondition}, :);
                        end
                    end
                    trialIdxInfo(iBasis, :, iRepeat) = trialIdx;
                end
            end

            numTrialsActual = nan(size(trialIdxInfo(:, :, 1)));
            for i = 1:size(trialIdxInfo, 1)
                for j = 1:size(trialIdxInfo, 2)
                    numTrialsActual(i, j) = numel(trialIdxInfo{i, j, 1});
                end
            end

            if p.Results.truncateNumTrialsToMaxOverBases
                numTrials = min(numTrials, max(numTrialsActual(:)));
            end

            if p.Results.fillNanToNumTrials
                % expand each cell to numTrials with NaN. The original
                % number of trials can be found using
                % indexInfo.trial
                for i = 1:numel(dataByTrialSubset)
                    sz = size(dataByTrialSubset{i});
                    if sz(1) < numTrials
                        dataByTrialSubset{i} = cat(1, dataByTrialSubset{i}, nan(numTrials - sz(1), sz(2)));
                    end
                end
            end

            indexInfo.numTrials = numTrials;
            indexInfo.numTrialsActual = numTrialsActual;
            indexInfo.basis = basisIdx;
            indexInfo.basisValidMask = pset.basisValid(basisIdx);
            indexInfo.align = alignIdx;
            indexInfo.condition = conditionIdx;
            indexInfo.trial = trialIdxInfo;

            if sameTrialsEachBasis
                indexInfo.trialShared = TensorUtils.squeezeDim(trialIdxInfo{1, :, :}, 1);
            else
                indexInfo.trialShared = [];
            end
        end

        % utilities for building single trial tensor not grouped by condition
        function [eachAlign_NbyRbyT, indexInfo, tvecCell] = arrangeEachAlign_NbyRbyT(pset, varargin)
            % eachAlign_NbyRbyT will be nAlign x 1 cell of nBases x numTrials x nTimeDataMean(iAlign) tensors
            % tvecCell will be tvecDataMean by default, unless timeDelta is specified
            p = inputParser();
            p.addParameter('validTimepointsAllBasesOnly', false, @islogical); % keep only timepoints where all valid bases have data
            p.addParameter('timeDelta', pset.timeDelta, @isscalar); % allow for data to be resampled in time
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [dataByTrial, indexInfo] = pset.computeDataByTrialSubset(p.Unmatched, ...
                'byCondition', false, ...
                'commonTime', true, ...
                'fillNanToNumTrials', true);

            tvecCell = pset.tvecDataMean;

            % dataByTrial is nBases x nAlign cell of R x T_a
            % data will be nAlign x 1 cell of N x R x T_a
            eachAlign_NbyRbyT = cellvec(indexInfo.nAlign);
            for iAlign = 1:numel(alignIdx)
                idxAlign = alignIdx(iAlign);

                % R x T x N --> N x R x T
                eachAlign_NbyRbyT{iAlign} = permute(cat(3, dataByTrial{:, idxAlign}), [3 1 2]);

                if p.Results.validTimepointsAllBasesOnly
                    % determine timepoints where all valid bases have data for all actual trials in the cell
                    rMaskSomeNonNaN = any(eachAlign_NbyRbyT{iAlign}, 2);
                    Tmask = all(all(~isnan(eachAlign_NbyRbyT{iAlign}(indexInfo.basisValidMask, rMaskSomeNonNaN, :)), 1), 2);
                    eachAlign_NbyRbyT{iAlign} = eachAlign_NbyRbyT{iAlign}(:, :, Tmask);
                end
            end

            if p.Results.timeDelta ~= pset.timeDelta
                % interp to new time delta
                for iA = 1:indexInfo.nAlign
                    [eachAlign_NbyRbyT{iA}, tvecCell{iA}] = TrialDataUtilities.Data.resampleTensorInTime(eachAlign_NbyRbyT{iA}, 3 ,...
                        tvecCell{iA}, 'timeDelta', p.Results.timeDelta);
                end
            end
        end

        function [NbyRbyTA, indexInfo, tvec, avec] = arrangeNbyRbyTA(pset, varargin)
            p = inputParser();
            % splicing and resampling
            p.addParameter('spliceAlignments', false, @islogical);
            p.addParameter('spliceOptions', struct(), @isstruct);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            [data, indexInfo, tvecCell] = pset.arrangeEachAlign_NbyRbyT(p.Unmatched);

            if p.Results.spliceAlignments
                % data is N x R x T --> N x T x R
                data = cellfun(@(d) permute(d, [1 3 2]), data, 'UniformOutput', false);
                data = TrialDataUtilities.Data.spliceTrajectories(data, 'basisMask', pset.basisValid(indexInfo.basis), p.Results.spliceOptions);
                NbyRbyTA = ipermute(data, [1 3 2]);
            else
                NbyRbyTA = TensorUtils.catWhich(3, data{:});
            end

            [tvec, avec] = TensorUtils.catWhich(1, tvecCell{:});
            avec = indexInfo.align(avec);
        end

        function [NbyTAbyR, indexInfo, tvec, avec] = arrangeNbyTAbyR(pset, varargin)
            [NbyRbyTA, indexInfo, tvec, avec] = pset.arrangeNbyRbyTA(pset, varargin{:});
            NbyTAbyR = permute(NbyRbyTA, [1 3 2]);
        end

        function [RTAbyN, indexInfo, rvec, tvec, avec, nvec] = arrangeRTAbyN(pset, varargin)
            % out is R*T*A x N concatenated matrix
            [NbyRbyTA, indexInfo, tvec, avec] = pset.arrangeNbyRbyTA(p.Unmatched);
            labels = {indexInfo.basis, indexInfo.trialShared, [tvec, avec]};
            % we include dim 4 as this would hold randomized samples
            [RTAbyN, labelsOut] = TensorUtils.reshapeByConcatenatingDims(NbyRbyTA, {[3 2], 1}, labels);

            rvec = labelsOut{1}(:, 1);
            tvec = labelsOut{1}(:, 2);
            avec = labelsOut{1}(:, 3);
            nvec = labelsOut{2};
        end

        % utilities for building single trial tensor grouped by condition
        function [eachAlign_RbyTbyNbyC, indexInfo] = arrangeEachAlign_RbyTbyNbyC(pset, varargin)
            [dataByTrialSubset, indexInfo] = pset.computeDataByTrialSubsetByCondition(varargin{:}, ...
                'commonTime', true, 'fillNanToNumTrials', true);
            nAlign = numel(indexInfo.align);
            eachAlign_RbyTbyNbyC = cell(nAlign, 1);
            for iAlign = 1:nAlign
                % start by making cell 1 x 1 x N x C x Rep
                tempCell = permute(dataByTrialSubset(:, iAlign, :, :), [5 6 1 3 4 2]);
                % collapse into R x T x N x C x Rep matrix
                eachAlign_RbyTbyNbyC{iAlign} = cell2mat(tempCell);
            end
        end

        % used in particular for dpca-type noise floor
        function [eachAlign_NbyTbyCbyR, indexInfo] = arrangeEachAlign_NbyTbyCbyR(pset, varargin)
            [eachAlign, indexInfo] = pset.arrangeEachAlign_RbyTbyNbyC(varargin{:});
            for iAlign = 1:numel(eachAlign)
                eachAlign{iAlign} = permute(eachAlign{iAlign}, [3 2 4 1 5]);
            end
            eachAlign_NbyTbyCbyR = eachAlign;
        end

        function [meansExcluding_NbyTAbyCbyR, trials_NbyTAbyCbyR, nTrials_NbyC, tvec, avec] = computeDataMeansExcludingSampledTrials(pset, varargin)
            p = inputParser();
            p.addParameter('maxTrials', Inf, @isscalar);
            p.addParameter('conditionIdx', truevec(pset.nConditions), @isvector);
            p.addParameter('alignIdx', truevec(pset.nAlign), @isvector);
            p.addParameter('basisIdx', truevec(pset.nBases), @isvector);
            p.addParameter('validBasesOnly', false, @islogical);
            p.addParameter('nRepeats', 1, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            alignIdx = TensorUtils.vectorMaskToIndices(p.Results.alignIdx);
            basisIdx = TensorUtils.vectorMaskToIndices(p.Results.basisIdx);
            conditionIdx = makecol(p.Results.conditionIdx);

            if p.Results.validBasesOnly
                mask = pset.basisValid(basisIdx);
                basisIdx = basisIdx(mask);
            end

            [trials_NbyTAbyCbyR, nTrials_NbyC, tvec, avec] = pset.arrangeNbyTAbyCbyTrials(varargin{:});

            if ~pset.hasDataByTrial
                % retrieve means excluding trial manually
                % then mask bases, align, and condition
                meansExcluding_NbyTAbyCbyR = pset.dataCachedMeanExcludingSampledTrialsTensor(basisIdx, :, conditionIdx, :);
                dsplit = TensorUtils.splitAlongDimension(meansExcluding_NbyTAbyCbyR, 2, pset.nTimeDataMean);
                dsplit = dsplit(alignIdx);
                meansExcluding_NbyTAbyCbyR = cat(2, dsplit{:});

                if p.Results.maxTrials < size(meansExcluding_NbyTAbyCbyR, 4)
                    meansExcluding_NbyTAbyCbyR = meansExcluding_NbyTAbyCbyR(:, :, :, 1:p.Results.maxTrials, :);
                end
            else
                % manually compute means without the trials chosen by
                % finding the sum, subtracting that trial, and
                % renormalizing

                % nAlign cellvec with N x C x T --> N x C x TA --> N x TA x C
                dataMeanTensor = permute(cat(3, pset.dataMean{:}), [1 3 2]);

                % A x N x C --> TA x N x C --> N x TA x C
                dataNumTrialsRaw = permute(repmat(max(pset.dataNumTrialsRaw, [], 1), [sum(pset.nTimeDataMean) 1 1]), [2 1 3]);

                dataSumTensor = dataMeanTensor .* dataNumTrialsRaw;

                % compute new mean by taking Sum - Sampled Trial / (Ntrials - 1)
                meansExcluding_NbyTAbyCbyR = bsxfun(@rdivide, bsxfun(@minus, dataSumTensor, trials_NbyTAbyCbyR), dataNumTrialsRaw-1);
            end
        end

        % added primarily for dpca-type noise floor
        function [NbyTAbyAttrbyR, nTrials_NbyAttr, tvec, avec, whichTrials_NbyAttr] = arrangeNbyTAbyConditionAttributesbyTrials(pset, varargin)
            [NbyTAbyCbyR, nTrialsTensor, tvec, avec, whichTrials] = pset.arrangeNbyTAbyCbyTrials(varargin{:});
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

        function [NbyTAbyRbyAttr, nTrials_NbyAttr, tvec, avec, whichTrials_NbyAttr] = arrangeNbyTAbyTrialsbyConditionAttributes(pset, varargin)
            [NbyTAbyCbyR, nTrialsTensor, tvec, avec, whichTrials] = pset.arrangeNbyTAbyCbyTrials(varargin{:});
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

        function nTrials_NbyAttr = computeTrialCountsNbyConditionAttr(pset)
            nTrials_NbyC = pset.dataNumTrials;
            condSize = pset.conditionDescriptor.conditionsSize;
            nTrials_NbyAttr = reshape(nTrials_NbyC, [size(nTrials_NbyC, 1) makerow(condSize)]);
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
