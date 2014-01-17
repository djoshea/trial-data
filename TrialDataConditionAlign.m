classdef(ConstructOnLoad) TrialDataConditionAlign < TrialData

    % ConditionInfo and AlignInfo stores
    properties(SetAccess=protected)
        % ConditionInfo instance
        conditionInfo
        
        % AlignInfo instance
        alignInfo
    end
    
    % On Demand Cache handle
    properties(Transient, Access=protected)
        odc % TrialDataConditionAlignOnDemandCache instance
    end
    
    % properties which are stored in odc, see get/set below
    properties(Transient, Dependent, SetAccess=protected)
        alignSummary
    end

    % Properties which read through to ConditionInfo
    properties(Dependent, SetAccess=protected)
        nConditions
        listByCondition
        conditionIdx
        conditionSubs
        conditions
        conditionsAsStrings
        conditionNames
        conditionAppearances
        
        axisValueLists
        axisValueListsAsStrings
        
        conditionAppearanceFn
    end

    % Properties which read through to AlignInfo
    properties(Dependent, SetAccess=protected)
        minTimeDelta
    end
    
    % Initializing and building
    methods
        function td = TrialDataConditionAlign(varargin)
            td = td@TrialData(varargin{:});
            td.odc = TrialDataConditionAlignOnDemandCache();
            td = td.initializeConditionInfo();
            td = td.initializeAlignInfo();
            td = td.updateValid();
        end
    end
    
    % Get / set accessors that write through to ODC
    methods
        function v = get.alignSummary(td)
            v = td.odc.alignSummary;            
            if isempty(v)
                td.buildAlignSummary();
                v = td.odc.alignSummary;
            end
        end
        
        function td = set.alignSummary(td, v)
            td.odc = td.odc.copy();
            td.odc.alignSummary = v;
        end
    end
    
    % build* methods for properties stored in odc
    methods
       function buildAlignSummary(td)
            alignSummary = AlignSummary.buildFromConditionAlignInfo(...
                td.conditionInfo, td.alignInfo);
            
            c = td.odc;
            c.alignSummary = alignSummary;
       end
    end
    
    % General utilites
    methods
        % print a short description
        function disp(td)
            td.printDescriptionShort();
            
            td.alignInfo.printOneLineDescription();
            td.conditionInfo.printOneLineDescription();
            
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end

        % synchronize valid between AlignInfo and ConditionINfo
        function td = updateValid(td)
            td.warnIfNoArgOut(nargout);
            cvalid = td.conditionInfo.computedValid;
            avalid = td.alignInfo.computedValid;
            if isempty(td.manualValid)
                td.manualValid = truevec(td.nTrials);
            end

            valid = td.manualValid & cvalid & avalid;

            td.conditionInfo = td.conditionInfo.setInvalid(~valid);
            td.alignInfo = td.alignInfo.setInvalid(~valid);
        end
        
        function valid = buildValid(td)
            valid = td.conditionInfo.valid & td.alignInfo.valid;
            if ~isempty(td.manualValid)
                valid = valid & td.manualValid;
            end
        end
        
        function td = dropChannels(td, names)
            names = wrapCell(names);
            
            % check whether any of the alignInfo events and error if so
            alignEvents = td.alignInfo.getEventList();
            mask = ismember(names, alignEvents);
            if any(mask)
                error('TrialData alignment depends on event %s', ...
                    strjoin(alignEvents(mask)));
            end
            
            % remove from condition info
            conditionParams = td.conditionInfo.attributeNames;
            mask = ismember(names, conditionParams);
            if any(mask)
                warning('TrialData condition depends on params %s, removing from grouping axes', ...
                    strjoin(conditionParams(mask)));
                td.conditionInfo = td.conditionInfo.removeAttribute(names(mask));
            end
            
            td = dropChannels@ConditionDescriptor(td, names);
            
            % in case we lost some event channels, update the alignInfo
            td = td.applyAlignInfo();
        end
    end

    % ConditionInfo control
    methods
        % parameters that are either scalar or strings
        function names = listConditionInfoCompatibleParamChannels(td)
            channelDescriptors = td.getChannelDescriptorArray();
            mask = arrayfun(@(cd) isa(cd, 'ParamChannelDescriptor') && ...
                (cd.isScalarByField{1} || cd.isStringByField{1}), ...
                channelDescriptors);
            names = {channelDescriptors(mask).name}';
        end
        
        function paramStruct = getConditionInfoCompatibleParamStruct(td)
            names = td.listConditionInfoCompatibleParamChannels();
            paramStruct = keepfields(td.data, names);
        end
        
        function td = initializeConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            if isempty(td.conditionInfo)
                paramStruct = emptyStructArray(td.nTrials);
                td.conditionInfo = ConditionInfo.fromStruct(paramStruct);
            end
        end
        
        function td = addEvent(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = addEvent@TrialData(td, varargin{:});
            td = td.applyAlignInfo();
        end
        
        function td = addChannel(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = addChannel@TrialData(td, varargin{:});
            
            % Don't do this anymore, we no longer auto-add attributes to
            % conditionInfo until they are needed
            
%             % detect whether any new condition info compatible params
%             % have been added
%             namesOld = td.listConditionInfoCompatibleParamChannels();
%             td = addChannel@TrialData(td, varargin{:});
%             names = td.listConditionInfoCompatibleParamChannels();
%             
%             % if so, add them to the condition info with valueLists
%             % specified
%             newAttr = setdiff(names, namesOld);
%             for iA = 1:numel(newAttr)
%                 td.conditionInfo = td.conditionInfo.addAttribute(newAttr{iA}, ...
%                     'values', td.getParam(newAttr{iA}));
%             end
        end
        
        function td = setConditionDescriptor(td, cd)
            td.warnIfNoArgOut(nargout);
            
            % manually accept a condition descriptor instance
            assert(isequal(class(cd), 'ConditionDescriptor'), 'Must be a ConditionDescriptor instance');
            
            % grab the param data to feed to the condition descriptor
            paramData = td.getRawChannelDataAsStruct(cd.attributeRequestAs);
            
            % build condition info from condition descriptor
            td.conditionInfo = ConditionInfo.fromConditionDescriptor(cd, paramData);
            td = td.postUpdateConditionInfo();
        end

        function td = selectTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            td = selectTrials@TrialData(td, mask);
            td.conditionInfo = td.conditionInfo.selectTrials(mask);
            td.alignInfo = td.alignInfo.selectTrials(mask);
            td.alignSummary = [];
        end
        
        function td = addAttribute(td, names)
            % add attributes in names that aren't already in ConditionInfo
            new = setdiff(names, td.conditionInfo.attributeNames);
            
            for i = 1:numel(new)
                td.conditionInfo = td.conditionInfo.addAttribute(new{i}, ...
                    'values', td.getParamRaw(new{i}));
            end
        end
        
        function td = groupBy(td, varargin)
            td.warnIfNoArgOut(nargout);
            
            % add any needed attributes to condition info
            for i = 1:numel(varargin)
                td = td.addAttribute(varargin{i});
            end
            
            td.conditionInfo = td.conditionInfo.groupBy(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = ungroup(td)
            td.warnIfNoArgOut(nargout);
            td = td.groupBy();
        end
        
        % will undo any filtering by attribute value lists and removing
        % binning
        function td = setAllAttributeValueListsAuto(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAllAttributeValueListsAuto();
            td = td.postUpdateConditionInfo();
        end
        
        function td = setAttributeValueListAuto(td, attr)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.setAttributeValueListAuto(attr);
            td = td.postUpdateConditionInfo();
        end
        
        function td = binAttribute(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttribute(varargin{:});
            td = td.postUpdateConditionInfo();
        end
            
        function td = binAttributeUniform(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttributeUniform(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = binAttributeQuantiles(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(varargin{1});
            td.conditionInfo = td.conditionInfo.binAttributeQuantiles(varargin{:});
            td = td.postUpdateConditionInfo();
        end
        
        function td = setAttributeValueList(td, attrName, valueList)
            td.warnIfNoArgOut(nargout);
            td = td.addAttribute(attrName);
            td.conditionInfo = td.conditionInfo.setAttributeValueList(attrName, valueList);
        end
        
        function valueList = getAttributeValueList(td, attrName)
            td = td.addAttrbute(attrName);
            valueList = td.conditionInfo.getAttributeValueList(attrName);
        end
        
        function td = postUpdateConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
            td.alignSummary = [];
        end
        
        % filter trials that are valid based on ConditionInfo
        function td = filterValidTrialsConditionInfo(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.selectTrials(td.conditionInfo.valid);
        end

        % filter trials where a specific attribute matches a specific value
        % pass-thru to ConditionInfo
        function td = filteredByAttribute(td, varargin)
            % call on ConditionInfo, and filter based on its valid trials
            td.warnIfNoArgOut(nargout);
            % no need to copy, filteredBy already copies the conditionInfo
            td.conditionInfo = td.conditionInfo.filteredByAttribute(varargin{:});
            td = td.filterValidTrialsConditionInfo();
        end

        % filter trials based on matching attribute values, 
        % pass thru to ConditionInfo
        function td = filteredByAttributeStruct(td, varargin)
            % call on ConditionInfo, and filter based on its valid trials
            td.warnIfNoArgOut(nargout);
            % no need to copy, filteredBy already copies the conditionInfo
            td.conditionInfo = td.conditionInfo.filteredByAttributeStruct(varargin{:});
            td = td.filterValidTrialsConditionInfo();
        end

        function tdCell = getTrialDataByCondition(td)
            tdCell = arrayfun(@(condition) td.filteredByAttributeStruct(condition), ...
                td.conditions, 'UniformOutput', false);
        end

        function v = get.conditions(td)
            v = td.conditionInfo.conditions;
        end
        
        function v = get.conditionsAsStrings(td)
            v = td.conditionInfo.conditionsAsStrings;
        end

        function n = get.nConditions(td)
            n = td.conditionInfo.nConditions;
        end
        
        function v = get.axisValueLists(td)
            v = td.conditionInfo.axisValueLists;
        end
        
        function v = get.axisValueListsAsStrings(td)
            v = td.conditionInfo.axisValueListsAsStrings;
        end

        function v = get.conditionNames(td)
            v = td.conditionInfo.names;
        end

        function v = get.conditionAppearances(td)
            v = td.conditionInfo.appearances;
        end
        
        function td = set.conditionAppearanceFn(td, v)
            td.conditionInfo.appearanceFn = v;
        end
        
        function v = get.conditionAppearanceFn(td)
            v = td.conditionInfo.appearanceFn;
        end
        
        function v = get.listByCondition(td)
            v = td.conditionInfo.listByCondition;
        end

        function v = get.conditionIdx(td)
            v = td.conditionInfo.conditionIdx;
        end

        function v = get.conditionSubs(td)
            v = td.conditionInfo.conditionSubs;
        end
        
        function minTimeDelta = get.minTimeDelta(td)
            minTimeDelta = td.alignDescriptor.minTimeDelta;
        end
    end

    % Data access by group via ConditionInfo
    methods
        % given a cellvec or nmeric vector, group its elements
        function varargout = groupElements(td, varargin)
            varargout = cell(nargout, 1);
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) data(idx,:), td.listByCondition, ...
                    'UniformOutput', false);
            end
        end
        
        % given data with dimension 1 with size nTrials, group by condition
        % and map out{i} = fn(group{i})
        function out = mapByGroup(td, fn, varargin)
            dataByGroup = td.groupElements(varargin{:});
            out = cellfun(fn, dataByGroup, 'UniformOutput', false);
        end

        function [dCell, tCell] = getAnalogGrouped(td, name)
            [dataCell, timeCell] = td.getAnalog(name);
            [dCell, tCell] = td.groupElements(dataCell, timeCell);
        end
        
        function [dataCell, timeCell] = getAnalogSampleGrouped(td, name, varargin)
            [dataVec, timeVec] = td.getAnalogSample(name);
            [dataCell, timeCell] = td.groupElements(dataVec, timeVec);
        end

        function dCell = getEventGrouped(td, name)
            dCell = td.groupElements(td.getAnalog(name));
        end

        function dCell = getParamGrouped(td, name)
            dCell = td.groupElements(td.getParam(name));
        end
    end

    % AlignInfo control
    methods
        function td = initializeAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            if isempty(td.alignInfo)
                td.alignInfo = AlignInfo();
                td = td.applyAlignInfo();
            end
        end
        
        function td = applyAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.applyToTrialData(td);
            td = td.postUpdateAlignInfo();
        end
        
        function td = postUpdateAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
            td.alignSummary = [];
        end

        function td = align(td, ad)
            td.warnIfNoArgOut(nargout);
            if ischar(ad)
                ad = AlignInfo(ad);
            else
                assert(isa(ad, 'AlignDescriptor'));
                if ~isa(ad, 'AlignInfo')
                    ad = AlignInfo.fromAlignDescriptor(ad);
                end
            end

            td.alignInfo = ad;
            td = td.applyAlignInfo();
        end
        
        function td = unalign(td)
            td.warnIfNoArgOut(nargout);
            td = td.align('TrialStart:TrialEnd');
        end
        
        % the following methods pass-thru to alignInfo:
        
        function td = pad(td, varargin)
            % add a padding window to the AlignInfo
            % may change which trials are valid
            % usage: pad([pre post]) or pad(pre, post)
            % pre > 0 means add padding before the start (typical case)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.pad(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = padForSpikeFilter(td, sf)
            % Pad trial data alignment for spike filter
            td.warnIfNoArgOut(nargout);
            td = td.pad([sf.preWindow sf.postWindow]);
        end
        
        function td = round(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.round(varargin{:});
        end
        
        function td = noRound(td)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.noRound();
        end
        
        function td = start(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.start(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = stop(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.stop(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = zero(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.zero(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = mark(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.mark(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = interval(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.interval(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = truncateBefore(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.truncateBefore(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = truncateAfter(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.truncateAfter(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = invalidateOverlap(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.invalidateOverlap(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = setOutsideOfTrialTruncate(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.setOutsideOfTrialTruncate(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = setOutsideOfTrialInvalidate(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.setOutsideOfTrialInvalidate(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        function td = setOutsideOfTrialIgnore(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.setOutsideOfTrialIgnore(varargin{:});
            td = td.postUpdateAlignInfo();
        end
        
        % filter trials that are valid based on AlignInfo
        function td = filterValidTrialsAlignInfo(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.selectTrials(td.alignInfo.computedValid);
        end
        
        % get the time window for each trial
        function durations = getValidDurations(td)
            durations = td.alignInfo.getValidDurationByTrial();
        end
        
        function [tMinByTrial, tMaxByTrial] = getTimeStartStopEachTrial(td)
            [tMinByTrial, tMaxByTrial] = td.alignInfo.getStartStopRelativeToZeroByTrial();
        end
    end
    
    % Aligned data access via AlignInfo
    methods
        function offsets = getTimeOffsetsFromZeroEachTrial(td)
            offsets = td.alignInfo.getZeroByTrial();
        end
        
        % return aligned analog channel
        function [data, time] = getAnalog(td, name)
            [data, time] = getAnalog@TrialData(td, name);
            [data, time] = td.alignInfo.getAlignedTimeseries(data, time, false);
        end
        
        function [mat, tvec] = getAnalogAsMatrix(td, name, varargin)
            % return aligned analog channel, resampled and interpolated to
            % a uniformly spaced time vector around t=0 such that the
            % result can be embedded in a nTrials x nTime matrix. time will
            % be chosen to encapsulate the min / max timestamps across all
            % trials. Missing samples will be returned as NaN
            
            p = inputParser;
            p.addParameter('timeDelta', [], @isscalar);
            p.parse(varargin{:});
            
            timeDelta = p.Results.timeDelta;
            if isempty(timeDelta)
                timeDelta = td.alignInfo.minTimeDelta;
                if isempty(timeDelta)
                    timeDelta = td.getAnalogTimeDelta(name);
                    warning('timeDelta auto-computed from analog timestamps. Specify manually or call .round for consistent results');
                end
            end
            
            [dataCell, timeCell] = td.getAnalog(name);
            
            [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeCell, ...
                'timeDelta', timeDelta, 'timeReference', 0, ...
                'interpolate', true);
        end 
        
        % return aligned event times
        function timesCell = getEvent(td, name)
            timesCell = getEvent@TrialData(td, name);
            timesCell = td.alignInfo.getAlignedTimes(timesCell, false);
        end

        % return aligned unit spike times
        function [timesCell] = getSpikeTimes(td, unitName)
            timesCell = getSpikeTimes@TrialData(td, unitName);
            timesCell = td.alignInfo.getAlignedTimes(timesCell, true);
        end
    end

    % Spike data
    methods
        function sr = buildSpikeRaster(td, unitName)
            sr = SpikeRaster(td, unitName, 'conditionInfo', td.conditionInfo, 'alignInfo', td.alignInfo);
            sr.useWidestCommonValidTimeWindow = false;
        end
        
        function [rateCell, timeCell] = getSpikeRateFiltered(td, unitName, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            spikeCell = td.getSpikeTimes(unitName);
            timeInfo = td.alignInfo.timeInfo;
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSec);
        end
           
        function [rates, tvec] = getSpikeRateFilteredAsMatrix(td, unitName, varargin)
            p = inputParser;
            p.addParameter('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.addParameter('timeDelta', 1, @isscalar);
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            timeDelta = p.Results.timeDelta;
            
            % Pad trial data alignment for spike filter
            td = td.pad([sf.preWindow sf.postWindow]);
            
            spikeCell = td.getSpikeTimes(unitName);
            timeInfo = td.alignInfo.timeInfo;
            
            % convert to .zero relative times since that's what spikeCell
            % will be in (when called in this class)
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rates, tvec] = sf.filterSpikeTrainsWindowByTrialAsMatrix(spikeCell, ...
                tMinByTrial, tMaxByTrial, td.timeUnitsPerSecond, ...
                'timeDelta', timeDelta);
        end

        function [rateCell, timeCell] = getSpikeRateFilteredGrouped(td, unitName, varargin)
            [rateCell, timeCell] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElements(rateCell);
            timeCell = td.groupElements(timeCell);
        end
        
        function [rateCell, tvec] = getSpikeRateFilteredAsMatrixGrouped(td, unitName, varargin)
            [rates, tvec] = td.getSpikeRateFiltered(unitName, varargin{:});
            rateCell = td.groupElements(rates);
        end
        
        function [psthMatrix, tvec, semMatrix] = getSpikeRateFilteredMeanByGroup(td, unitName, varargin)
            [rateCell, tvec] = getSpikeRateFilteredGrouped(td, unitName, varargin{:});
            psthMatrix = cell2mat(cellfun(@(r) nanmean(r, 1), rateCell, 'UniformOutput', false));
            semMatrix =  cell2mat(cellfun(@(r) nansem(r, 1),  rateCell, 'UniformOutput', false));
        end
        
        function timesCellofCells = getSpikeTimesGrouped(td, unitName)
            timesCell = td.getSpikeTimes(unitName);
            timesCellofCells = td.groupElements(timesCell);
        end
        
        function countsCell = getSpikeCountsGrouped(td, unitName)
            counts = td.getSpikeCounts(unitName);
            countsCell = td.groupElements(counts);
        end
        
        function rateCell = getSpikeMeanRateGrouped(td, unitName)
            rates = td.getSpikeMeanRate(unitName);
            rateCell = td.groupElements(rates);
        end
    end

    % Plotting
    methods
        function plotAnalogGroupedEachTrial(td, name, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.addParamValue('patchline', false, @islogical); % use patchline to enable transparency
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            [dataByGroup, timeByGroup] = td.getAnalogGrouped(name);     
            app = td.conditionAppearances;

            for iCond = 1:td.nConditions
                dataCell = dataByGroup{iCond};
                timeCell = timeByGroup{iCond};
                for iTrial = 1:numel(dataCell)
                    if ~isempty(timeCell{iTrial}) && ~isempty(dataCell{iTrial})
                        if p.Results.patchline
                           patchline(double(timeCell{iTrial}), dataCell{iTrial}, ...
                               'FaceColor', app(iCond).color, 'FaceAlpha', 0.3, ...
                               'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                        else
                            plot(axh, double(timeCell{iTrial}), dataCell{iTrial}, '-', 'Color', app(iCond).color, ...
                                'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                        end
                    end
                    if iTrial == 1, hold(axh, 'on'); end
                end
            end
            box(axh, 'off');
            axis(axh, 'tight');
            
            xlabel(td.getTimeAxisLabel());
            ylabel(td.getAxisLabelForChannel(name));
        end
        
        function plotAnalogGroupedEachTrial2D(td, name1, name2, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.addParamValue('patchline', false, @islogical); % use patchline to enable transparency
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            dataByGroup1 = td.getAnalogGrouped(name1);  
            dataByGroup2 = td.getAnalogGrouped(name2);  
            app = td.conditionAppearances;

            for iCond = 1:td.nConditions
                dataCell1 = dataByGroup1{iCond};
                dataCell2 = dataByGroup2{iCond};
                for iTrial = 1:numel(dataCell1)
                    if p.Results.patchline
                        patchline(dataCell1{iTrial}, dataCell2{iTrial}, ...
                               'EdgeColor', app(iCond).color, ...
                               'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                    else
                        plot(axh, dataCell1{iTrial}, dataCell2{iTrial}, '-', 'Color', app(iCond).color, ...
                            'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                    end
                    if iTrial == 1, hold(axh, 'on'); end
                end
            end
            box(axh, 'off');
            axis(axh, 'tight');
            
            xlabel(td.getAxisLabelForChannel(name1));
            ylabel(td.getAxisLabelForChannel(name2));
        end
        
        function plotAnalogGroupedEachTrial3D(td, name1, name2, name3, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            dataByGroup1 = td.getAnalogGrouped(name1);  
            dataByGroup2 = td.getAnalogGrouped(name2);
            dataByGroup3 = td.getAnalogGrouped(name3);
            app = td.conditionAppearances;

            for iCond = 1:td.nConditions
                dataCell1 = dataByGroup1{iCond};
                dataCell2 = dataByGroup2{iCond};
                dataCell3 = dataByGroup3{iCond};
                for iTrial = 1:numel(dataCell1)
                    plot3(axh, dataCell1{iTrial}, dataCell2{iTrial}, dataCell3{iTrial}, ...
                        '-', 'Color', app(iCond).color, ...
                        'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
                    if iTrial == 1, hold(axh, 'on'); end
                end
            end
            box(axh, 'off');
            axis(axh, 'tight');
            axis(axh, 'vis3d');
            
            xlabel(td.getAxisLabelForChannel(name1));
            ylabel(td.getAxisLabelForChannel(name2));
            zlabel(td.getAxisLabelForChannel(name3));
        end
    end

end
