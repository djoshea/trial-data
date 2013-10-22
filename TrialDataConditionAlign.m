classdef TrialDataConditionAlign < TrialData

    properties(SetAccess=protected)
        conditionInfo
        alignInfo
    end

    % Properties which read through to ConditionInfo
    properties(Dependent)
        nConditions
        listByCondition
        conditionIdx
        conditionSubs
        conditions
        conditionNames
        conditionAppearances
        
        conditionAppearanceFn
    end

    % Initializing and building
    methods
        function td = TrialDataConditionAlign(varargin)
            td = td@TrialData(varargin{:});
            td = td.initializeConditionInfo();
            td = td.initializeAlignInfo();
            td = td.updateValid();
        end
    end
    
    methods
        function disp(td)
            td.printDescriptionShort();
            fprintf('\n');
            
            td.alignInfo.printOneLineDescription();
            fprintf('\n');
            td.conditionInfo.printDescription();
            
            fprintf('\n');
            td.printChannelInfo();
            fprintf('\n');
        end

        function td = updateValid(td)
            td.warnIfNoArgOut(nargout);
            cvalid = td.conditionInfo.computedValid;
            avalid = td.alignInfo.computedValid;

            % combine and update the validity masks 
            td.valid = cvalid & avalid;
            td.conditionInfo.setInvalid(~td.valid);
            td.alignInfo = td.alignInfo.setInvalid(~td.valid);
        end
    end

    % ConditionInfo control
    methods
        function td = initializeConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = ConditionInfo();
            td.conditionInfo = td.conditionInfo.applyToTrialData(td);
        end

        function td = selectTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            td = selectTrials@TrialData(td, mask);
            td.conditionInfo = td.conditionInfo.selectTrials(mask);
            td.alignInfo = td.alignInfo.selectTrials(mask);
        end
        
        function td = groupBy(td, varargin)
            td.warnIfNoArgOut(nargout);
            
            td.conditionInfo = td.conditionInfo.groupBy(varargin{:});

            td = td.postUpdateConditionInfo();
        end
        
        function td = setAttributeValueList(td, attrName, valueList)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo.setValueList(attrName, valueList);
        end
        
        function valueList = getAttributeValueList(td, attrName)
            valueList = td.conditionInfo.getAttributeValueList(attrName);
        end
        
        function td = postUpdateConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
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

        function n = get.nConditions(td)
            n = td.conditionInfo.nConditions;
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
    end

    % Data access by group via ConditionInfo
    methods
        % given a cellvec or nmeric vector, group its elements
        function varargout = groupElements(td, varargin)
            for i = 1:numel(varargin)
                data = varargin{i};
                assert(size(data,1) == td.nTrials, ...
                    'Data must have size nTrials along 1st dimension');
                varargout{i} = cellfun(@(idx) data(idx,:), td.listByCondition, 'UniformOutput', false);
            end
        end

        function [dCell tCell] = getAnalogGrouped(td, name)
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
            td.alignInfo = AlignInfo();
            td.alignInfo = td.alignInfo.applyToTrialData(td);
        end
        
        function td = postUpdateAlignInfo(td)
            td.warnIfNoArgOut(nargout);
            td = td.updateValid();
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

            td.alignInfo = ad.applyToTrialData(td);
            td = td.postUpdateAlignInfo();
        end
        
        % add a padding window to the AlignInfo
        % may change which trials are valid
        % usage: pad([pre post]) or pad(pre, post)
        % pre > 0 means add padding before the start (typical case)
        function td = pad(td, varargin)
            td.warnIfNoArgOut(nargout);
            td.alignInfo = td.alignInfo.pad(varargin{:});
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
    end
    
    % AlignInfo data access
    methods
        % return aligned analog channel
        function [data time] = getAnalog(td, name)
            [data time] = getAnalog@TrialData(td, name);
            [data time] = td.alignInfo.getAlignedTimeseries(data, time);
        end
        
        % return aligned event times
        function timesCell = getEvent(td, name)
            timesCell = getEvent@TrialData(td, name);
            timesCell = td.alignInfo.getAlignedTimes(timesCell);
        end

        % return aligned unit spike times
        function [timesCell] = getSpikeTimesForUnit(td, unitName)
            timesCell = getSpikeTimesForUnit@TrialData(td, unitName);
            timesCell = td.alignInfo.getAlignedTimes(timesCell);
        end
        
     
        
    end

    % Spike data
    methods
        function sr = buildSpikeRasterForUnit(td, unitName)
            sr = SpikeRaster(td, unitName, 'conditionInfo', td.conditionInfo, 'alignInfo', td.alignInfo);
            sr.useWidestCommonValidTimeWindow = false;
        end
        
        function [rates, tvec] = getFilteredSpikeRateEachTrial(td, unitName, varargin)
            p = inputParser;
            p.addParamValue('spikeFilter', SpikeFilter.getDefaultFilter(), @(x) isa(x, 'SpikeFilter'));
            p.parse(varargin{:});
            
            sf = p.Results.spikeFilter;
            spikeCell = td.getSpikeTimesForUnit(unitName);
            timeInfo = td.alignInfo.timeInfo;
            
            % convert to .zero relative times since that's what spikeCell
            % will be in
            tMinByTrial = [timeInfo.start] - [timeInfo.zero];
            tMaxByTrial = [timeInfo.stop] - [timeInfo.zero];
            [rates, tvec] = sf.filterSpikeTrainsWindowByTrial(spikeCell, tMinByTrial, tMaxByTrial);
        end
        
        function [rateCell, tvec] = getFilteredSpikeRateGroupedEachTrial(td, unitName, varargin)
            [rateMat, tvec] = td.getFilteredSpikeRateEachTrial(unitName, varargin{:});
            rateCell = td.groupElements(rateMat);
        end
        
        function timesCellofCells = getSpikeTimesForUnitGrouped(td, unitName)
            timesCell = td.getSpikeTimesForUnit(unitName);
            timesCellofCells = td.groupElements(timesCell);
        end
        
        function countsCell = getSpikeCountsForUnitGrouped(td, unitName)
            counts = td.getSpikeCountsForUnit(unitName);
            countsCell = td.groupElements(counts);
        end
        
        function rateCell = getSpikeRatePerSecForUnitGrouped(td, unitName)
            rates = td.getSpikeRatePerSecForUnit(unitName);
            rateCell = td.groupElements(rates);
        end
    end

    % Plotting
    methods
        function plotAnalogGroupedEachTrial(td, name, varargin) 
            p = inputParser();
            p.addParamValue('plotOptions', {}, @(x) iscell(x));
            p.KeepUnmatched;
            p.parse(varargin{:});

            axh = td.getRequestedPlotAxis(p.Unmatched);

            [dataByGroup timeByGroup] = td.getAnalogGrouped(name);     
            app = td.conditionAppearances;

            for iCond = 1:td.nConditions
                dataCell = dataByGroup{iCond};
                timeCell = timeByGroup{iCond};
                for iTrial = 1:numel(dataCell)
                    if ~isempty(timeCell{iTrial}) && ~isempty(dataCell{iTrial})
                        plot(axh, double(timeCell{iTrial}), dataCell{iTrial}, '-', 'Color', app(iCond).color, ...
                            'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
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
                    plot(axh, dataCell1{iTrial}, dataCell2{iTrial}, '-', 'Color', app(iCond).color, ...
                        'LineWidth', app(iCond).lineWidth, p.Results.plotOptions{:});
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
