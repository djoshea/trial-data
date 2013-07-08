classdef TrialDataConditionAlign < TrialData

    properties 
        conditionInfo
        alignInfo
    end

    % Properties which read through to ConditionInfo
    properties
        listByCondition
        conditionIdx
        conditionSubs
        conditions
        conditionNames
        conditionAppearances
    end

    % Initializing and building
    methods
        function td = TrialDataConditionAlign(varargin)
            td = td@TrialData(varargin{:});
            td = td.initializeConditionInfo();
        end
    end
    
    methods
        function disp(td)
            td.conditionInfo.printDescription();
            builtin('disp', td);
        end
    end

    % ConditionInfo
    methods
        function td = initializeConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = ConditionInfo(td.data);
            td.conditionInfo.getAttributeValueFn = @TrialDataConditionAlign.getAttributeFn;
        end

        function td = createCopyConditionInfo(td)
            td.warnIfNoArgOut(nargout);
            td.conditionInfo = td.conditionInfo.copy();
        end

        function td = selectTrials(td, mask)
            td.warnIfNoArgOut(nargout);
            td = selectTrials@TrialData(td, mask);
            td = td.createCopyConditionInfo();
            td.conditionInfo.resampleTrials(mask);
        end
        
        function td = addAttributeToConditionInfo(td, paramList, createCopy)
            td.warnIfNoArgOut(nargout);
            if nargin < 3
                createCopy = true;
            end
            if createCopy
                % TrialData are not handle classes, copy conditionInfo
                td = td.createCopyConditionInfo();
            end

            if ischar(paramList)
                paramList = {paramList};
            end

            addMask = ~td.conditionInfo.hasAttribute(paramList);
            paramList = paramList(addMask);

            for i = 1:numel(paramList)
                param = paramList{i};
                matrix = td.channelDescriptorsByName.(param).scalar;
                td.conditionInfo.addAttribute(paramList(addMask), 'matrix', matrix);
            end
        end

        function td = groupBy(td, paramList, varargin)
            td.warnIfNoArgOut(nargout);
            p = inputParser;
            p.addRequired('paramList', @(x) ischar(x) || iscellstr(x));
            p.parse(paramList, varargin{:});

            paramList = p.Results.paramList;
            if ischar(paramList)
                paramList = {paramList};
            end

            % TrialData are not handle classes, copy conditionInfo
            td = td.createCopyConditionInfo();

            td = addAttributeToConditionInfo(td, paramList, false);

            td.conditionInfo.groupBy(paramList);
        end

        % filter trials that are valid based on ConditionInfo
        function td = filterValidTrialsConditionInfo(td, varargin)
            td.warnIfNoArgOut(nargout);
            td = td.createCopyConditionInfo();
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

        function v = get.conditionNames(td)
            v = td.conditionInfo.names;
        end

        function v = get.conditionAppearances(td)
            v = td.conditionInfo.appearances;
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

    % ConditionInfo callback methods, not bound to this class
    methods(Static)
        % return a scalar struct with one field for each attribute containing the attribute values
        % as a cell array or numeric vector
        function valueStruct = getAttributeFn(data, attributeNames, varargin)
            %debug('fetching attributes %s\n', strjoin(attributeNames));
            assert(isstruct(data), 'getAttributeFn expects data as struct()');
            for iAttr = 1:length(attributeNames)
                attr = attributeNames{iAttr};
                [valueStruct(1:numel(data)).(attr)] = data.(attr);
            end
        end
    end

end
