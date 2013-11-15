classdef PopulationTrajectorySet < handle & matlab.mixin.Copyable 
% class for storing a set of population trajectories, one for each of multiple conditions
% The bases may be units (e.g. trial-averaged firing rates by condition by unit), or 
% components of a low-dimensional projection (e.g. PCA)
%
% The set of conditions is specified by a ConditionDescriptor
% Common time alignment is specified by an AlignDescriptor

    properties(SetAccess=protected)
        % if true, the raw sources used to build the data will be maintained within
        % this instance, which allows things like resampling and shuffling
        % to operate on the raw data sources
        hasDataSources 

        % a cell array of unknown length which stores the union of all trial data 
        % instances used by all bases
        trialDataSources

        % a pointer list of size nBases x 1 which indicates which element of 
        % trialDataSources a particular basis derives its data 
        trialDataIndByBasis

    end

    properties %(SetAccess=?StateSpaceProjection)
        % unit names 
        basisNames = {};

        % unit meta data
        basisMeta = {};

        % nBases x nConditions x nAlign cell array of aligned time vectors 
        timeData = {};
        
        % nBases x nConditions x nAlign cell array to firing rate traces
        data = {};
        
        % nBases x nConditions x nAlign logical array indicating whether to include
        % the conditions across 
        dataValid
        
        % nBases x nConditions x nAlign scalar array indicating how many
        % trials contributed to data
        nTrialsData

        
        % active time window, respected by getDataTimeWindowed, but doesn't
        % actually truncate data
        % nBases x nConditions x nAlign
        tMinDataManual
        tMaxDataManual
    end
    
    properties % these properties have set methods which check the new value
        % alignment information (AlignDescriptor cell array)
        alignDescriptorSet = {};

        % condition information (ConditionDescriptor)
        conditionDescriptor
    end

    properties(Dependent)
        % number of separate alignDescriptors
        nAlign

        nBases

        nConditions
        
        dataSize % [nBases nConditions nAlign]
        
        nConditionsValid % number of conditions which are valid in all bases, all aligns

        conditionsSize % pass-thru to .conditionDescriptor
        
        alignNames

        conditionNames 
        
        % nConditions x nAlign scalar array indicating whether this element
        % is valid across all bases
        conditionAlignsValidAllBases
        
        % nConditions array indicating whether this condition has valid
        % data across all bases and all alignments
        conditionsValidAllBasesAlign
    end
    
    properties(Constant) % how .data, .timeData are organized
        DIM_BASES = 1;
        DIM_CONDITIONS = 2;
        DIM_ALIGN = 3;
    end

    methods % Constructor 
        function pset = PopulationTrajectorySet(varargin)
            % will need to fix this to copy ad's 
            p = inputParser;
            p.addParamValue('alignDescriptor', [], @(x) isa(x, 'AlignDescriptor'));
            p.addParamValue('alignDescriptorSet', [], @iscell);
            p.addParamValue('conditionDescriptor', [], @(x) isa(x, 'ConditionDescriptor'));
            p.addParamValue('storeDataSources', false, @islogical);
            p.parse(varargin{:});

            if ~isempty(p.Results.alignDescriptor)
                % use only one alighDescriptor
                pset.alignDescriptorSet = { p.Results.alignDescriptor };
            elseif ~isempty(p.Results.alignDescriptorSet)
                pset.alignDescriptorSet = p.Results.alignDescriptorSet;
            end

            pset.conditionDescriptor = p.Results.conditionDescriptor;

            pset.storeDataSources = p.Results.storeDataSources;
        end
    end

    methods(Access=protected) % copyElement for deep copy
        % deep copy all handle class properties
        function pset2 = copyElement(pset)
            pset2 = copyElement@matlab.mixin.Copyable(pset);
        end
    end

    methods % Filtering bases, conditions
        function filterAlign(pset, idx)
            p = inputParser;
            p.addRequired('alignIdx', @isvector);
            p.parse(idx);
            alignIdx = makecol(p.Results.alignIdx);

            pset.data = pset.data(:, :, alignIdx);
            pset.timeData = pset.timeData(:, :, alignIdx);
            pset.alignTimeInfoData = pset.alignTimeInfoData(:, :, alignIdx);
            pset.alignDescriptorSet = pset.alignDescriptorSet(alignIdx);
            pset.dataValid = pset.dataValid(:, :, alignIdx);
            pset.nTrialsData = pset.nTrialsData(:, :, alignIdx);
            pset.tMinDataManual = pset.tMinDataManual(:, :, alignIdx);
            pset.tMinDataManual = pset.tMaxDataManual(:, :, alignIdx);
           
            if pset.storeDataSources
                pset.dataSources = pset.dataSources(:, alignIdx);
                pset.dataSourcesOrig = pset.dataSourcesOrig(:, alignIdx);
            end
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
            pset.nTrialsData = pset.nTrialsData(basisIdx, :, :);
            
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
%             pset.nTrialsData = pset.nTrialsData(:, mask, :);
% 
%             pset.tMinDataManual = pset.tMinDataManual(:, mask, :);
%             pset.tMaxDataManual = pset.tMaxDataManual(:, mask, :);
%         end
        
%         function updateConditionDescriptor(pset)
%             pset.conditionDescriptor = pset.conditionDescriptor.updateCache();
%         end
    end

    methods % Dependent properties
        function n = get.nAlign(pset)
            n = length(pset.alignDescriptorSet);
        end

        function n = get.nBases(pset)
            n = size(pset.data, 1);
        end
        
        function sz = get.conditionsSize(pset)
            if isempty(pset.conditionDescriptor)
                sz = [0 0];
            else
                sz = pset.conditionDescriptor.conditionsSize;
            end
        end

        function n = get.nConditions(pset)
            if isempty(pset.conditionDescriptor)
                n = 0;
            else
                n = pset.conditionDescriptor.nConditions;
            end
        end
        
        function sz = get.dataSize(pset)
            sz = [pset.nBases pset.nConditions pset.nAlign];
        end
        
        function n = get.nConditionsValid(pset)
            if isempty(pset.conditionDescriptor)
                n = 0;
            else
                n = nnz(pset.conditionsValidAllBasesAlign);
            end
        end
        
        function valid = get.conditionAlignsValidAllBases(pset)
            valid = shiftdim(any(pset.dataValid, 1), 1);
        end
        
        function valid = get.conditionsValidAllBasesAlign(pset)
            valid = any(pset.conditionAlignsValidAllBases, 2);
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
        
        function set.alignDescriptorSet(pset, adCell)
            assert(iscell(adCell), 'alignDescriptorSet must be a cell array');
            
            nAlign = length(adCell);
            if pset.nBases > 0 && pset.nAlign > 0
                assert(nAlign == pset.nAlign, 'Number of alignments must match');
            end
            
            pset.alignDescriptorSet = adCell;
        end
        
        function set.conditionDescriptor(pset, cd)
%             if ~isempty(pset.conditionDescriptor) && pset.nBases > 0
%                 assert(isequal(cd.conditionsSize, pset.conditionsSize), 'Condition tensor size must match existing ConditionDescriptor');
%             end
            
            pset.conditionDescriptor = cd;
        end
        
%         function t = get.tMinDataManual(pset)
%             if isempty(pset.tMinDataManual)
%                 pset.tMinDataManual = nan(pset.dataSize);
%             end
%             t = pset.tMinDataManual;
%         end
%         
%         function t = get.tMaxDataManual(pset)
%             if isempty(pset.tMaxDataManual)
%                 pset.tMaxDataManual = nan(pset.dataSize);
%             end
%             t = pset.tMaxDataManual;
%         end
    end

    methods % Adding bases / units
        function initialize(pset)
            % create appropriately sized empty arrays for storing data,
            % timeData, alignTimeInfo, etc.
            if isempty(pset.conditionDescriptor)
                error('Set .conditionDescriptor before calling initialize()');
            end
            if isempty(pset.alignDescriptorSet)
                error('Set .alignDescriptorSet before calling initialize()');
            end
            
            nAlign = length(pset.alignDescriptorSet);
            [pset.timeData, pset.data, pset.alignTimeInfoData] = ...
                deal(cell(0, pset.nConditions, nAlign));
            pset.dataValid = false(0, pset.nConditions, nAlign);
            pset.nTrialsData = nan(0, pset.nConditions, nAlign);
            pset.tMinDataManual = nan(0, pset.nConditions, pset.nAlign);
            pset.tMaxDataManual = nan(0, pset.nConditions, pset.nAlign);
                        
            % initialize data source store if asked
            if pset.storeDataSources
                pset.dataSources = cell(0, nAlignAdd);
                pset.dataSourcesOrig = cell(0, nAlignAdd);
            end
            
            pset.basisNames = cell(0, 1);
            pset.basisMeta = cell(0, 1);
        end
        
%         function addEmptyBases(pset, nBasesAdd)
%             if nargin < 2
%                 nBasesAdd = 1;
%             end
%             
%             expand = @(in) TensorUtils.expandAlongDims(in, pset.DIM_BASES, nBasesAdd);
%             
%             pset.data = expand(pset.data);
%             pset.timeData = expand(pset.timeData);
%             pset.alignTimeInfoData = expand(pset.alignTimeInfoData);
%             pset.dataValid = expand(pset.dataValid);
%             pset.nTrialsData = cat(1, pset.nTrialsData, nTrialsData);
%                 
%             pset.basisNames = cat(1, pset.basisNames, names);
%             pset.basisMeta = cat(1, pset.basisMeta, meta);   
%         end 
%         
        function initializeFromSpikeRaster(pset, srCell)
            % set alignDescriptorSet and conditionDescriptor based off a
            % SpikeRaster or cell vector of SpikeRasters
            % then call intiialize
            if ~iscell(srCell)
                srCell = {srCell};
            end
            
            srCell = makecol(srCell);
            pset.conditionDescriptor = srCell{1}.conditionInfo.getConditionDescriptor();
            pset.alignDescriptorSet = cellfun(@(sr) sr.alignDescriptor, srCell, 'UniformOutput', false);
            
            pset.initialize();
        end
        
        function addUnitFromSpikeRaster(pset, srCell, varargin)
            % srCell should be a cell array of nUnitsToAdd x nAlign of spike rasters
            % if there is only one alignment, ensure srCell is a column vector of units
            % if only one alignment and only one unit to add, no need to wrap in {}
            p = inputParser;
            p.addRequired('srCell', @(x) iscell(x) || isa(x, 'SpikeRaster'));
            p.addParamValue('name', [], @ischar); 
            p.addParamValue('meta', [], @(x) true);
            p.parse(srCell, varargin{:});

            if ~iscell(srCell)
                srCell = {srCell};
            end

            nUnitsAdd = size(srCell, 1);
            nAlignAdd = size(srCell, 2);

            names = p.Results.name;
            if isempty(names)
                names = cell(nUnitsAdd, 1);
                for i = 1:nUnitsAdd
                    names{i} = srCell{i}.name;
                    if isempty(names{i});
                        names{i} = sprintf('Unit %d', pset.nBases + i);
                    end
                end
            end

            if isempty(pset.alignDescriptorSet)
                % haven't been initialized yet, learn the nAlign from the srCell
                pset.alignDescriptorSet = cell(nAlignAdd, 1);
            else
                assert(nAlignAdd == pset.nAlign, ...
                    'Please provide a cell array of SpikeRasters corresponding to the %d AlignDescriptors used in this PopulationTrajectorySet', ...
                    pset.nAlign);
            end
            
            nBases = pset.nBases;
                
            % do error checking and validation first, then store data
            prog = ProgressBar(nUnitsAdd, 'Validating compatible AlignDescriptors, ConditionDescriptors...');
            for iUnit = 1:nUnitsAdd
                prog.update(iUnit);
                for iAlign = 1:nAlignAdd
                    sr = srCell{iUnit, iAlign};

                    % fetch or validate conditionDescriptor
                    cdFromSR = sr.conditionInfo.getConditionDescriptor();
                    if isempty(pset.conditionDescriptor)
                        % use the condition Descriptor from this raster
                        pset.conditionDescriptor = cdFromSR;
                    else
                        % check that the condition descriptor matches
                        assert(isequal(cdFromSR, pset.conditionDescriptor), ...
                            'ConditionDescriptor does not match this SpikeRaster''s');
                    end

                    % fetch or validate alignDescriptor 
                    adFromSR = sr.alignInfo;
         
                    if isempty(adFromSR.name)
                        adFromSR.name = sprintf('Align %d', iAlign);
                    end
                    if isempty(pset.alignDescriptorSet{iAlign})
                        % use the align descriptor from this raster
                        pset.alignDescriptorSet{iAlign} = adFromSR;

                    else
                        % check that the align descriptor matches
                        assert(adFromSR.isCompatibleWith(pset.alignDescriptorSet{iAlign}), ...
                            'AlignDescriptor does not match this SpikeRaster''s');
                    end

                end
            end
            prog.finish();

            prog = ProgressBar(nUnitsAdd, 'Adding bases from SpikeRasters');
            for iUnit = 1:nUnitsAdd
                prog.update(iUnit);
                iBasis = nBases + iUnit;
                for iAlign = 1:nAlignAdd
                    sr = srCell{iUnit, iAlign};
                    [psthByCondition, ~, timeBins] = sr.getPSTHByCondition();
                    
                    if iBasis == 1 && iAlign == 1
                        % TODO eventually want to instantiate this in its own function
                        [pset.timeData pset.data pset.alignTimeInfoData] = deal(cell(1, pset.nConditions, nAlignAdd));
                        pset.dataValid = false(1, pset.nConditions, nAlignAdd);
                        pset.nTrialsData = nan(1, pset.nConditions, nAlignAdd);
                        
                        % initialize data source store if asked
                        if pset.storeDataSources
                            pset.dataSources = cell(1, nAlignAdd);
                            pset.dataSourcesOrig = cell(1, nAlignAdd);
                        end
                    end

                    % add the align time info for each trial to facilitate proper time axis drawing
                    % could make this more efficient by changing AlignDescriptor to accept the medians directly
                    alignTimeInfoByCondition = sr.getAlignTimeInfoByCondition();
                    nTrialsByCondition = sr.nTrialsByCondition;
                    for iCondition = 1:pset.nConditions
                        if nTrialsByCondition(iCondition) > 1
                            % crucial that everything stays a column vector for cell2mat to work appropriately
                            pset.timeData{iBasis, iCondition, iAlign} = makecol(timeBins);
                            pset.data{iBasis, iCondition, iAlign} = makecol(psthByCondition(iCondition, :));
                        end
                        
                        pset.alignTimeInfoData{iBasis, iCondition, iAlign} = makecol(alignTimeInfoByCondition{iCondition});    
                        
                        % store whether this spike raster contributed any trials
                        pset.dataValid(iBasis, :, iAlign) = nTrialsByCondition(:) > 0;
                        pset.nTrialsData(iBasis, :, iAlign) = nTrialsByCondition(:);
                    end

                    name = p.Results.name;
                    if isempty(name)
                        name = srCell{1}.name;
                    end
                    pset.basisNames{iBasis} = name;
                    pset.basisMeta{iBasis} = p.Results.meta;

                    % hold on to the SpikeRaster if requested
                    if pset.storeDataSources
                        pset.dataSources{iBasis, iAlign} = sr;
                        pset.dataSourcesOrig{iBasis, iAlign} = sr;
                    end
                end
                
                pset.tMinDataManual = cat(1, pset.tMinDataManual, nan(1, pset.nConditions, pset.nAlign));
                pset.tMaxDataManual = cat(1, pset.tMaxDataManual, nan(1, pset.nConditions, pset.nAlign));
            end
            prog.finish();
            
        end
        
        function addBasisFromRawData(pset, dataCell, timeCell, alignTimeInfoCell, varargin)
            % Add raw data directly to PopulationTrajectorySet
            %
            % dataCell, timeCell, alignTimeInfoCell :
            %    nBases x nConditions x nAlign cell tensors containing the
            %    appropriate data vector for that basis x condition x
            %    alignment.
            p = inputParser;
            p.addRequired('dataCell', @(x) iscell(x));
            p.addRequired('timeCell', @(x) iscell(x));
            % nUnitsToAdd x nAlign x nCondition cell
            p.addRequired('alignTimeInfoCell', @(x) iscell(x));
            p.addParamValue('names', [], @(x) ischar(x) || iscell(x)); 
            p.addParamValue('meta', [], @iscell);
            p.parse(dataCell, timeCell, alignTimeInfoCell, varargin{:});

            nUnitsAdd = size(dataCell, 1);
            nConditionsAdd = size(dataCell, 2);
            nAlignAdd = size(dataCell, 3);
            
            names = p.Results.names;
            if isempty(names) || ~iscell(names)
                names = cell(nUnitsAdd, 1);
                for i = 1:nUnitsAdd         
                    names{i} = sprintf('Basis %d', pset.nBases + i);
                end
            end
            names = makecol(names);
            assert(iscellstr(names), isvector(names) && length(names) == nUnitsAdd, ...
                'Names must be a cellstr vector whose length matches size(dataCell, 1)');
            
            meta = p.Results.meta;
            meta = makecol(meta);
            if isempty(meta)
                meta = cell(nUnitsAdd, 1);
            end
            assert(iscell(meta) && isvector(meta) && length(meta) == nUnitsAdd, ...
                'Meta must be a cell vector whose length matches size(dataCell, 1)');
            
            if ~isempty(pset.alignDescriptorSet)
                assert(nAlignAdd == pset.nAlign, ...
                    'size(dataCell, 3) must match .nAlign==%d', ...
                    pset.nAlign);
            else
                error('Set .alignDescriptorSet first');
            end
            
            if isempty(pset.conditionDescriptor)
                error('Set .conditionDescriptor first');
            else
                assert(nConditionsAdd == pset.nConditions, ...
                    'size(dataCell, 2) must match .nConditions==%d', pset.nConditions);
            end  
            
            nBases = pset.nBases;
            
            assert(isequal(size(dataCell), size(timeCell)), ...
                'Size of dataCell and timeCell must match');
            assert(isequal(size(dataCell), size(alignTimeInfoCell)), ...
                'Size of dataCell and alignTimeInfoCell must match');
            
            % make sure everything is a column vector for cell2mat to work
            % correctly
            dataCell = cellfun(@makecol, dataCell, 'UniformOutput', false);
            timeCell = cellfun(@makecol, timeCell, 'UniformOutput', false);
            alignTimeInfoCell = cellfun(@makecol, alignTimeInfoCell, 'UniformOutput', false);
            dataValid = cellfun(@(x) any(~isnan(x)), dataCell);
            nTrialsData = nan(nUnitsAdd, nConditionsAdd, nAlignAdd);
            
            pset.data = cat(1, pset.data, dataCell);
            pset.timeData = cat(1, pset.timeData, timeCell);
            pset.alignTimeInfoData = cat(1, pset.alignTimeInfoData, alignTimeInfoCell);
            pset.dataValid = cat(1, pset.dataValid, dataValid);
            pset.nTrialsData = cat(1, pset.nTrialsData, nTrialsData);
            
            pset.tMinDataManual = cat(1, pset.tMinDataManual, nan(nUnitsAdd, pset.nConditions, pset.nAlign));
            pset.tMaxDataManual = cat(1, pset.tMaxDataManual, nan(nUnitsAdd, pset.nConditions, pset.nAlign));
                
            pset.basisNames = cat(1, pset.basisNames, names);
            pset.basisMeta = cat(1, pset.basisMeta, meta);   
        end
    end

    methods(Static)
        % take all units from each trial data as bases, not collected simultaneously
        function pset = buildFromTrialDataSetAllUnits(tdCell, varargin)
            p = inputParser();
            p.parse(varargin{:});

            nTD = numel(tdCell);
            unitsByTD = cellfun(@(td) td.listAllUnits(), tdCell);
            nBasesByTD = cellfun(@numel, unitsByTD);

            iBasis = 0;
            for iTD = 1:nTD
                for iUnitThisTD = 1:nBasesByTD(iTD)
                    iBasis = iBasis + 1;
                    td.trialDataSet{iBasis} = tdCell{iTD};
                end 
            end
        end

        function pset = buildfromRawData(dataCell)    
            timeCell = cellfun(@(d) (1:numel(d))', dataCell, 'UniformOutput', false);
            
            alignTimeInfoCell = cellfun(@(d) struct('valid', true, 'start', 0, ...
                'stop', numel(d)-1, 'startPad', 0, 'stopPad', numel(d)-1, 'zero', 0, ...
                'mark', {{}}, 'interval', {{}}), dataCell, 'UniformOutput', false);
            
            nConditions = size(dataCell, 2);
            cd = ConditionDescriptor();
            cd  = cd.addAttribute('condition', 'valueList', 1:nConditions, 'groupBy', true);
            
            nAlign = size(dataCell, 3);
            for i = 1:nAlign
                alignDescriptorSet{i} = AlignDescriptor('start:stop');
            end
            
            pset = PopulationTrajectorySet('conditionDescriptor', cd, 'alignDescriptorSet', alignDescriptorSet);
            pset.addBasisFromRawData(dataCell, timeCell, alignTimeInfoCell);
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

        function updateFromDataSources(pset)
            % TODO this should work with sources besides SpikeRaster, and should probably
            % be utilized by the addUnitsFrom spike rasters
            
            for iBasis = 1:pset.nBases
                for iAlign = 1:pset.nAlign
                    % get psth across conditions
                    sr = pset.dataSources{iBasis, iAlign};
                    [psthByCondition, ~, time] = sr.getPSTHByCondition();

                    pset.data(iBasis, :, iAlign) = mat2cell(psthByCondition, ones(pset.nConditions, 1), length(time));
                    [pset.timeData{iBasis, :, iAlign}] = deal(makecol(time));
                    pset.alignTimeInfoData(iBasis, :, iAlign) = TensorUtils.flatten(sr.getAlignTimeInfoByCondition());
                end
            end
        end

        function restoreFromOriginalDataSources(pset)
            % TODO this should work with sources besides SpikeRaster, and should probably
            % be utilized by the addUnitsFrom spike rasters
            pset.dataSources = pset.dataSourcesOrig;
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
            
            tMinData = TensorUtils.repmatSliceAlongDims(tMinByAlign, pset.dataSize, pset.DIM_ALIGN);
            tMaxData = TensorUtils.repmatSliceAlongDims(tMaxByAlign, pset.dataSize, pset.DIM_ALIGN);
            
            pset.constrainTimeWindowManual(tMinData, tMaxData)
        end
        
        function constrainTimeWindowManual(pset, tMinData, tMaxData)
            % further constrain the existing .tMinDataManual and .tMaxDataManual
            pset.tMinDataManual = nanmin(pset.tMinDataManual, tMinData);
            pset.tMaxDataManual = nanmin(pset.tMaxDataManual, tMaxData);
        end
        
        function removeTimeWindowManual(pset)
            [pset.tMinDataManual, pset.tMaxDataManual] = deal(nan(pset.dataSize));
        end
        
        function applyTimeWindowData(pset, varargin)
            % permanently truncate data outside of tMinData(...) and
            % tMaxData(...)
            [pset.data, pset.timeData]= pset.getDataTimeWindowed(varargin{:});
        end
        
        function [tMinData tMaxData] = getTimeValidData(pset, varargin)
            % return the largest time window which excludes nan values in psthData
            % by basis x condition x align
            tMinData = pset.alignBasisConditionDataFun(@getTimeValidStart, 'asMatrix', true);
            tMaxData = pset.alignBasisConditionDataFun(@getTimeValidEnd, 'asMatrix', true);
        
            % constrain to lie within manual time windows
            tMinData = nanmax(tMinData, pset.tMinDataManual);
            tMaxData = nanmin(tMaxData, pset.tMaxDataManual);
            
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

        function [tMinData tMaxData] = getTimeValidAcrossBasesConditions(pset, varargin)
            % finds the global time valid window for each alignment, 
            % ignoring bases that are entirely nan for the interval
            
            [tMinData tMaxData] = pset.getTimeValidData();

            for iAlign = 1:pset.nAlign 
                tMin = tMinData(:, :, iAlign);
                tMinData(:, :, iAlign) = nanmax(tMin(:));  % nanmax ignores all nan bases
                tMax = tMaxData(:, :, iAlign);
                tMaxData(:, :, iAlign) = nanmin(tMax(:));
            end
        end

        function [tMinData tMaxData] = getTimeValidAcrossConditions(pset, varagin)
            % tMin, tMax are N x C x A
            [tMinData tMaxData] = pset.getTimeValidData();

            % determine the largest valid window that works across all conditions (2)
            tMinAcrossBases = nanmax(tMinData, [], 2);
            tMaxAcrossBases = nanmin(tMaxData, [], 2);

            % re-expand the time windows
            tMinData = repmat(tMinAcrossBases, [1 pset.nConditions 1]);
            tMaxData = repmat(tMaxAcrossBases, [1 pset.nConditions 1]);
            
        end

        function [tMinData tMaxData] = getTimeValidAcrossBases(pset, varargin)
            % tMin, tMax are N x C x A
            [tMinData tMaxData] = pset.getTimeValidData();

            % determine the largest valid window that works across all bases (dim 1)
            tMinAcrossBases = nanmax(tMinData, [], 1);
            tMaxAcrossBases = nanmin(tMaxData, [], 1);

            % re-expand the time windows
            tMinData = repmat(tMinAcrossBases, [pset.nBases 1 1]);
            tMaxData = repmat(tMaxAcrossBases, [pset.nBases 1 1]);
        end

        function [data timeData] = getDataTimeWindowed(pset, varargin)
            % getDataTimeWindowed(tMinData, tMaxData)
            % window the data within tMinData(iBasis, iCondition, iAlign) : tMaxData(b,c,a)
            % by default tMinData and tMaxData are determined by getTimeValidData()
            % if scalar, tM??Data will be expanded to nBases x nConditions x nAlign
            %
            % tMinData and tMaxData are arrays of nBasis x nCondition x nAlign
            % data and timeData are cell arrays of nBasis x nCondition x nAlign
            p = inputParser;
            p.addOptional('tMinData', [], @(x) isnumeric(x));
            p.addOptional('tMaxData', [], @(x) isnumeric(x));
            
            % if tMin / tMax are outside the bounds or data is invalid, 
            % fill with appropriately sized NaNs to populate?
            p.addOptional('fillWithNaN', false, @islogical); 
            
            p.parse(varargin{:});
            tMinData = p.Results.tMinData;
            tMaxData = p.Results.tMaxData;
            fillWithNaN = p.Results.fillWithNaN;

            % by default use the valid window on each basis x condition x align
            if isempty(tMinData) || isempty(tMaxData)
                [tMinTemp tMaxTemp] = getTimeValidData(pset);
                if isempty(tMinData)
                    tMinData = tMinTemp;
                end
                if isempty(tMaxData)
                    tMaxData = tMaxTemp;
                end
            end

            % expand to size in case scalar
            if size(tMinData) == 1 
                tMinData = repmat(tMinData, [pset.nBases, pset.nConditions, pset.nAlign]);
            end
            if size(tMaxData) == 1 
                tMaxData = repmat(tMaxData, [pset.nBases, pset.nConditions, pset.nAlign]);
            end

            [data timeData]= pset.alignBasisConditionDataFun(@getPSTHTimeWindowed, 'checkInvalid', false);
            
            return;

            function [d t] = getPSTHTimeWindowed(time, data, iBasis, iCondition, iAlign)
                tMin = tMinData(iBasis, iCondition, iAlign);
                tMax = tMaxData(iBasis, iCondition, iAlign);
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
            [tMinData tMaxData] = pset.getTimeValidAcrossConditions();

            [data timeData] = pset.getDataTimeWindowed(tMinData, tMaxData);

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
            
            [tMinData tMaxData] = pset.getTimeValidAcrossBasesConditions();

            [data timeData] = pset.getDataTimeWindowed(tMinData, tMaxData);

            % all tvecs in timeData for iAlign are the same
            tvecByAlign = squeezedim(timeData(1,1,:), [1 2]);
        end

        function [data timeData tvecByConditionAlign] = getDataTimeWindowedValidAcrossBases(pset, varargin)
            % window the data such that all bases share a common time window within each condition x align
            % data and timeData are cell arrays of nBasis x nCondition x nAlign
            % tvecByConditionAlign is a nCondition x nAlign cell vector with the time vector for each condition x align

            [tMinData tMaxData] = getTimeValidAcrossBases(pset, varargin);
            [data timeData] = pset.getDataTimeWindowed(tMinData, tMaxData); 

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
