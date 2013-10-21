classdef(HandleCompatible, ConstructOnLoad) ConditionDescriptor 
% ConditionDescriptor is a static representation of a A-dimensional combinatorial
% list of attribute values

    properties
        description = '';
        
        % updates cache on set
        nameFn % function which maps .values(i) struct --> name of condition i
        
        % updates cache on set
        appearanceFn; % function which takes struct('attrName1', attrVal1, 'attrName2', attrVal2)
                      % and returns struct('color', 'k', 'lineWidth', 2, ...);      
    end
     
    properties(SetAccess=protected)
        % A x 1 : by attribute                       
        attributeNames = {}; % A x 1 cell array : list of attributes for each dimension
        attributeRequestAs = {}; % A x 1 cell array : list of names by which each attribute should be requested corresponding to attributeNames
        
        axisAttributes % G x 1 cell : each is cellstr of attributes utilized along that grouping axis
    end
    
    properties(Access=protected)
        attributeNumeric = []; % A x 1 logical array : is this attribute a numeric value? 
        
        attributeValueListsManual = {}; % A x 1 cell array of permitted values (or cells of values) for this attribute
        attributeValueBinsManual = {}; % A x 1 cell array of value Nbins x 2 value bins to use for numeric lists
        attributeValueBinsAutoCount % A x 1 numeric array of Nbins to use when auto computing the bins, NaN if not in use
        attributeValueBinsAutoMode % A x 1 numeric array of either AttributeValueBinsAutoUniform or AttributeValueBinsAutoQuantiles
        
        axisValueListsManual % G x 1 cell of cells: each contains a struct specifying an attribute specification for each element along the axis
        axisValueListsOccupiedOnly % G x 1 logical indicating whether to constrain the combinatorial valueList to only occupied elements (with > 0 trials)
        axisRandomizeMode % G x 1 numeric of constants beginning with Axis* (see below)
    end
    
    properties(SetAccess=protected)
        % When re-arranging the axes, default condition appearances can get shuffled around
        % which can make comparison across figures difficult. These cache
        % the condition appearances to make things easier. Call
        % freezeAppearances to activate()
        appearanceFrozen = false;
        frozenAppearanceConditions
        frozenAppearanceData
    end
    
    % END OF STORED TO DISK PROPERTIES
    
    properties(Transient, Access=protected)
        odc % handle to a ConditionDescriptorOnDemandCache
    end
    
    % THE FOLLOWING PROPERTIES WRAP EQUIVALENT PROPERTIES IN ODC
    % on get: retrieve from odc, if empty {call build<Property>, store in odc, return it}
    % on set: make copy of odc to alleviate dependency, store in odc
    % 
    % Note: we use the build<Property> methods because property getters
    % cannot be inherited, so subclasses can override the build method
    % instead.
    properties(Transient, Dependent, SetAccess=protected)        
        % These are generated on the fly by property get, but cached for speed, see invalidateCache to reset them 
        
        % these are A-dimensional objects where A is nAttributesGroupBy or length(groupByList)
        conditions % A-dimensional struct where values(...idx...).attribute is the value of that attribute on that condition
        conditionsRelevantAttributesOnly
        appearances % A-dimensional struct of appearance values
        names % A-dimensional cellstr array with names of each condition 
        attributeValueLists % A x 1 cell array of values allowed for this attribute
                           % here just computed from attributeValueListManual, but in ConditionInfo
                           % can be automatically computed from the data
        attributeValueListsAsStrings % same as above, but everything is a string
        attributeDescriptions

        axisValueLists % G dimensional cell array of structs which select attribute values for that position along an axis
        axisValueListMode % G dimensional array of 
    end
    
    % how are attribute values determined for a given attribute?
    properties(Constant, Hidden)
        % for attributeValueListMode
        AttributeValueListManual = 1;
        AttributeValueListAuto = 2;
        AttributeValueBinsManual = 3;
        AttributeValueBinsAutoUniform = 4;
        AttributeValueBinsAutoQuantiles = 5;
        
        % for axisRandomizeMode
        AxisOriginal = 1;
        AxisShuffled = 2;
        AxisResampled = 3;
        AxisResampledFromFirst = 4;
        
        % for axisValueListMode
        AxisValueListAutoAll = 1;
        AxisValueListAutoOccupied = 2;
        AxisValueListManual = 3;
    end
        
    properties(Dependent, Transient)
        nAttributes % how many attributes: ndims(values)

        nValuesByAttribute % how many values per attribute: size(values)

        attributeAlongWhichAxis % A x 1 array indicating which axis an attribute contributes to (or NaN)
        attributeValueModes
        
        nConditions % how many total conditions
        
        conditionsSize 

        conditionsAsLinearInds % linear index corresponding to each condition if flattened 
        

        nAxes % how many dimensions of grouping axe
        nValuesAlongAxes % X x 1 array of number of elements along the axis
        axisDescriptions % strcell describing each axis
    end
    
    % Constructor, load, save methods
    methods
        function ci = ConditionDescriptor()
            ci.odc = ci.buildOdc();
        end
        
        function odc = buildOdc(varargin)
            odc = ConditionDescriptorOnDemandCache();
        end
    end

    % get, set data stored inside odc
    methods 
        function v = get.conditions(ci)
            v = ci.odc.conditions;
            if isempty(v)
                ci.odc.conditions = ci.buildConditions();
                v = ci.odc.conditions;
            end
        end
        
        function ci = set.conditions(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditions = v;
        end
        
        function v = get.conditionsRelevantAttributesOnly(ci)
            v = ci.odc.conditionsRelevantAttributesOnly;
            if isempty(v)
                ci.odc.conditionsRelevantAttributesOnly = ci.buildConditionsRelevantAttributesOnly();
                v = ci.odc.conditionsRelevantAttributesOnly;
            end
        end
        
        function ci = set.conditionsRelevantAttributesOnly(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionsRelevantAttributesOnly = v;
        end
        
        function v = get.appearances(ci)
            v = ci.odc.appearances;
            if isempty(v)
                ci.odc.appearances = ci.buildAppearances();
                v = ci.odc.appearances;
            end
        end
        
        function ci = set.appearances(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.appearances = v;
        end
        
        function v = get.names(ci)
            v = ci.odc.names;
            if isempty(v)
                ci.odc.names = ci.buildNames();
                v = ci.odc.names;
            end
        end 
        
        function ci = set.names(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.names = v;
        end
        
        function v = get.attributeValueLists(ci)
            v = ci.odc.attributeValueList;
            if isempty(v)
                ci.odc.attributeValueList = ci.buildAttributeValueLists();
                v = ci.odc.attributeValueList;
            end
        end
        
        function ci = set.attributeValueLists(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.attributeValueList = v;
        end
        
        function v = get.attributeValueListsAsStrings(ci)
            v = ci.odc.attributeValueListAsStrings;
            if isempty(v)
                ci.odc.attributeValueListAsStrings = ci.buildAttributeValueListAsStrings();
                v = ci.odc.attributeValueListAsStrings;
            end
        end
        
        function ci = set.attributeValueListsAsStrings(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.attributeValueListAsStrings = v;
        end

        function v = get.axisValueLists(ci)
            v = ci.odc.axisValueLists;
            if isempty(v)
                ci.odc.axisValueLists = ci.buildAxisValueLists();
                v = ci.odc.axisValueLists;
            end
        end

        function ci = set.axisValueLists(ci, v)
            ci.odc.axisValueLists = v;
        end
    end

    % build data stored inside odc
    methods 
        function values = buildConditionsRelevantAttributesOnly(ci)
            if ci.nAxes == 0
                values = struct();
            else
                valueLists = ci.axisValueLists; 
                values = TensorUtils.mapFromAxisLists(@structMergeMultiple,...
                    valueLists, 'asCell', false);
            end
        end
        
        function values = buildConditions(ci)
            values = ci.conditionsRelevantAttributesOnly;
            
            % and add "wildcard" match for all other attributes
            whichAxis = ci.attributeAlongWhichAxis;
            for iA = 1:ci.nAttributes
                if isnan(whichAxis(iA))
                    values = assignIntoStructArray(values, ci.attributeNames{iA}, ...
                        ci.attributeValueLists(iA));
                end
            end
        end

        function valueListByAxes = buildAxisValueLists(ci)
            valueListByAxes = cellvec(ci.nAxes);
            for iX = 1:ci.nAxes  
                % build a list of values as a struct array for this axis
                % valueList{iAxis}(iEl).attribute = values describes the attribute values
                % allowed for attribute `attribute` at position iEl along axis iAxis
               
                % G x 1 cell of cells: each contains a struct specifying an attribute specification for each element along the axis
                if isempty(ci.axisValueListsManual{iX})
                    % build auto list of attributes
                    valueListByAxes{iX} = makecol(ci.buildAutoValueListForAttributeSet(ci.axisAttributes{iX}));
                else
                    valueListByAxes{iX} = makecol(ci.axisValueListsManual{iX});
                end
            end
        end

        % build a struct array for a set of attributes that walks all possible combinations of the attribute values 
        function values = buildAutoValueListForAttributeSet(ci, attributes)
            if ischar(attributes)
                attributes = {attributes};
            end
            attrIdx = ci.getAttributeIdx(attributes);
            valueLists = ci.attributeValueLists(attrIdx);

            values = TensorUtils.mapFromAxisLists(@buildStruct, valueLists, ...
                'asCell', false);

            function s = buildStruct(varargin)
                for i = 1:numel(varargin)
                    s.(attributes{i}) = varargin{i};
                end
            end

        end

        function names = buildNames(ci)
            % pass along values(i) and the subscripts of that condition in case useful 
            if ci.nConditions > 0
                fn = ci.nameFn;
                if isempty(fn)
                    fn = @ConditionDescriptor.defaultNameFn;
                end

                wrapFn = @(varargin) fn(ci, varargin{:});
                names = TensorUtils.mapIncludeSubs(wrapFn, ci.conditions);
            else
                names = {};
            end
        end

        function appearances = buildAppearances(ci)
            if ci.nConditions > 0
                appearFn = ci.appearanceFn;
                defaultFn = eval(sprintf('@%s.defaultAppearanceFn', class(ci)));
                defaults = defaultFn(ci);

                if isempty(appearFn)
                    % use the default function built into ConditionDescriptor
                    % or whatever subclass version of
                    % defaultAppearanceFn there is (namely ConditionInfo)
                    appearances = defaults;
                else
                    appearances = appearFn(ci, defaults);
                    % ensure that no fields have been lost from the
                    % defaults
                    appearances = structMerge(defaults, appearances, 'warnOnOverwrite', false);
                end
            else
                appearances = struct([]);
            end
        end

        function valueList = buildAttributeValueLists(ci)
            % just pull the manual lists (ConditionInfo will deal
            modes = ci.attributeValueModes;
            valueList = cellvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                switch modes(i) 
                    case ci.AttributeValueListManual
                        valueList{i} = ci.attributeValueListsManual{i};
                    case ci.AttributeValueBinsManual
                        valueList{i} = ci.attributeValueBinsManual{i};
                    case {ci.AttributeValueBinsAutoUniform, ci.AttributeValueBinsAutoQuantiles}
                        % the number of bins is known, so they can be specified here
                        valueList{i} = 1:ci.attributeValueBinsAutoCount(i);
                    otherwise
                         % place holder, must be determined when
                        % ConditionInfo applies it to data
                        if ci.attributeNumeric(i)
                            valueList{i} = NaN;
                        else
                            valueList{i} = {'?'};
                        end
                end
                valueList{i} = makecol(valueList{i});
            end
        end
        
        function valueList = buildAttributeValueListAsStrings(ci)
            modes = ci.attributeValueModes;
            valueList = ci.attributeValueLists;
            for i = 1:ci.nAttributes
                switch modes(i) 
                    case ci.AttributeValueListManual
                        if ci.attributeNumeric(i)
                            valueList{i} = cellfun(@num2str, valueList{i}, 'UniformOutput', false);
                        end             
                    case ci.AttributeValueBinsManual
                        bins = ci.attributeValueBinsManual{i};
                        valueList{i} = arrayfun(@(row) sprintf('%d-%d', bins(row, 1), bins(row, 2)), ...
                            1:size(bins, 2), 'UniformOutput', false);
                    case ci.AttributeValueBinsAutoUniform
                        % placeholder for actual bin limits
                        valueList{i} = arrayfun(@(bin) sprintf('uniformBin %d', bin), 1:ci.attributeValueBinsAutoCount, 'UniformOutput', false);
                    case ci.AttributeValueBinsAutoQuantiles
                        valueList{i} = arrayfun(@(bin) sprintf('quantile %d', bin), 1:ci.attributeValueBinsAutoCount(i), 'UniformOutput', false);
                    otherwise
                        valueList{i} = {};
                        % auto list leave empty, must be determined when
                        % ConditionInfo applies it to data
                end
                valueList{i} = makecol(valueList{i});
            end
        end

        function valueList = getAttributeValueList(ci, name)
            idx = ci.getAttributeIdx(name);
            valueList = makecol(ci.attributeValueLists{idx});
        end

        function valueIdx = getAttributeValueIdx(ci, attr, value)
            [tf, valueIdx] = ismember(value, ci.getAttributeValueLists(attr));
            assert(tf, 'Value not found in attribute %s valueList', attr);
        end
    end

    methods % AXES
        function n = get.nAxes(ci)
            n = numel(ci.axisAttributes);
        end
        
        function a = get.attributeAlongWhichAxis(ci)
            a = nanvec(ci.nAttributes);
            for iX = 1:ci.nAxes
                a(ci.getAttributeIdx(ci.axisAttributes{iX})) = iX;
            end
        end
        
        function modes = get.axisValueListMode(ci)
            modes = nanvec(ci.nAxes);
            
            for iX = 1:ci.nAxes
                if ~isempty(ci.axisValueListsManual{iX})
                    modes(iX) = ci.AxisValueListManual;
                elseif ci.axisValueListsOccupiedOnly(iX)
                    modes(iX) = ci.AxisValueListAutoOccupied;
                else
                    modes(iX) = ci.AxisValueListAutoAll;
                end
            end
        end
        
        function desc = get.axisDescriptions(ci)
            desc = cellvec(ci.nAxes);
            
            for iX = 1:ci.nAxes
                attr = ci.axisAttributes{iX};
                nv = ci.conditionsSize(iX);
                switch ci.axisValueListMode(iX)
                    case ci.AxisValueListAutoAll
                        vlStr = ' auto';
                    case ci.AxisValueListAutoOccupied
                        vlStr = ' autoOccupied';
                    case ci.AxisValueListManual
                        vlStr = ' manual';
                end   
                        
                switch ci.axisRandomizeMode(iX)
                    case ci.AxisOriginal
                        randStr = '';
                    case ci.AxisShuffled
                        randStr = ' shuffled';
                    case ci.AxisResampled
                        randStr = ' resampled';
                    case ci.AxisResampledFromFirst
                        randStr = ' resampledFromFirst';
                end
                
                desc{iX} = sprintf('axis %d : %s%s (%d%s)', iX, ...
                    strjoin(attr, ' x '), randStr, nv, vlStr);
            end
        end
        
        function ci = addAxis(ci, varargin)
            ci.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addOptional('attributes', {}, @(x) ischar(x) || iscellstr(x));
            p.addParamValue('name', '', @ischar);
            p.addParamValue('valueList', [], @(x) true);
            p.parse(varargin{:});

            if ~iscell(p.Results.attributes)
                attr = {p.Results.attributes};
            else
                attr = p.Results.attributes;
            end
            ci.assertHasAttribute(attr);
            
            ci = ci.removeAttributesFromAxes(attr);

            % create a grouping axis
            idx = ci.nAxes + 1; 
            ci.axisAttributes{idx} = attr;
            ci.axisValueListsManual{idx} = p.Results.valueList;
            ci.axisRandomizeMode(idx) = ci.AxisOriginal;
            ci.axisValueListsOccupiedOnly(idx) = true;

            ci = ci.invalidateCache();
        end
        
        function ci = maskAxes(ci, mask)
            ci.warnIfNoArgOut(nargout);
            
            ci.axisAttributes = ci.axisAttributes(mask);
            ci.axisValueListsManual = ci.axisValueListsManual(mask);
            ci.axisRandomizeMode = ci.axisRandomizeMode(mask);
            ci.axisValueListsOccupiedOnly = ci.axisValueListsOccupiedOnly(mask);
            
            ci = ci.invalidateCache();
        end

        % wipe out existing axes and creates simple auto axes along each 
        function ci = groupBy(ci, varargin)
            ci.warnIfNoArgOut(nargout);
            if iscell(varargin{1})
                attributes = varargin{1};
            else
                attributes = varargin;
            end

            for i = 1:numel(attributes)
                ci = ci.addAxis(attributes{i});
            end
        end

        function ci = groupByAll(ci)
            ci.warnIfNoArgOut(nargout);
            ci = ci.groupBy(ci.attributeNames);
        end

        % remove all axes
        function ci = clearAxes(ci)
            ci.warnIfNoArgOut(nargout);

            ci = ci.maskAxes([]);

            ci = ci.invalidateCache();
        end
        
        function ci = removeAttributesFromAxes(ci, namesOrIdx)
            ci.warnIfNoArgOut(nargout);
            attrIdx = ci.getAttributeIdx(namesOrIdx);
            attrNames = ci.attributeNames(attrIdx);
            
            if ci.nAxes == 0
                return;
            end
            
            whichAxis = ci.attributeAlongWhichAxis;
            removeAxisMask = falsevec(ci.nAxes);
            for iAI = 1:numel(attrIdx)
                iA = attrIdx(iAI);
                iX = whichAxis(iA);
                if isnan(iX)
                    continue;
                end
                
                % remove this attribute from axis iX
                maskInAxis = strcmp(ci.axisAttributes{iX}, attrNames{iAI});
                if all(maskInAxis)
                    removeAxisMask(iX) = true;
                else
                    ci.axisAttributes{iX} = ci.axisAttributes{iX}(~maskInAxis);
                    % clear out manual value list as it's likely invalid now
                    ci.axisValueListsManual{iX} = [];
                    % and reset the randomization
                    ci.axisRandomizeMode(iX) = ci.AxisOriginal;
                end
            end
            
            ci = ci.maskAxes(~removeAxisMask);
        end
    end

    methods % General methods, setters and getters
        
        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function ci = invalidateCache(ci)
            ci.warnIfNoArgOut(nargout);

            % here we precompute these things to save time, 
            % but each of these things also has a get method that will
            % recompute this for us
            ci.odc  = ci.odc.copy();
            ci.odc.flush();
        end

        % Manually freeze the condition appearances so that they don't
        % change when we downsample the conditions
        function ci = freezeAppearances(ci)
            ci.warnIfNoArgOut(nargout);

            % freeze current appearance information
            ci.frozenAppearanceConditions = ci.conditions;
            ci.frozenAppearanceData = ci.appearances;
            ci.appearanceFn = @ConditionDescriptor.frozenAppearanceFn;
        end

        function ci = set.nameFn(ci, fn)
            ci.nameFn = fn;
            ci = ci.invalidateCache();
        end

        function ci = set.appearanceFn(ci, fn)
            ci.appearanceFn = fn;
            ci = ci.invalidateCache();
        end
        
        function nv = get.conditionsSize(ci)
            nv = size(ci.conditions);
        end

        function linearInds = get.conditionsAsLinearInds(ci)
            linearInds = TensorUtils.containingLinearInds(ci.conditionsSize);
        end

        function n = get.nConditions(ci)
            n = prod(ci.conditionsSize);
        end

        function printDescription(ci) 
            tcprintf('yellow', '%s:', class(ci));
            tcprintf('inline', '{bright blue}Attributes: {white}%s\n', strjoin(ci.attributeDescriptions));
            axisDesc = ci.axisDescriptions;
            for iX = 1:ci.nAxes
                tcprintf('inline', '\t%s\n', axisDesc{iX});
            end
        end

        function disp(ci)
            ci.printDescription();
            fprintf('\n');
            builtin('disp', ci);
        end
    end

    methods % Attribute manipulations
        function [tf, idx] = hasAttribute(ci, name)
            [tf, idx] = ismember(name, ci.attributeNames);
        end

        function idx = assertHasAttribute(ci, name)
            [tf, idx] = ci.hasAttribute(name);
            if ~all(tf)
                if iscell(name)
                    name = strjoin(name(~tf));
                end
                error('Attribute(s) %s not found', name);
            end
        end

        function na = get.nAttributes(ci)
            na = length(ci.attributeNames);
        end

                function idxList = getAttributeIdx(ci,name)
            if isempty(name)
                idxList = [];
                return;
            end
            
            if isnumeric(name)
                % already idx, just return
                idxList = floor(name);
                idxList(idxList < 0 | idxList > ci.nAttributes) = NaN;
                return;
            end
            
            if ~iscell(name)
                name = {name};
            end

            idxList = nan(length(name), 1);
            for i = 1:length(name)
                if ischar(name{i})
                    idx = find(strcmp(ci.attributeNames, name{i}), 1);
                else
                    idx = name{i};
                end
                if isempty(idx)
                    error('Cannot find attribute named %s', name{i});
                end
                idxList(i) = idx;
            end
        end
        
        function tf = getIsAttributeNumeric(ci, name)
            idx = ci.getAttributeIdx(name);
            tf = ci.attributeNumeric(idx);
        end
        
        % return an A x 1 numeric array of constants in the AttributeValue*
        % set listed above, describing how this attribute's values are
        % determined
        function modes = get.attributeValueModes(ci)
            % check for manual value list, then manual bins, then auto
            % bins, otherwise auto value list
            modes = nanvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                if ~isempty(ci.attributeValueListsManual{i})
                    modes(i) = ci.AttributeValueListManual;
                elseif ~isempty(ci.attributeValueBinsManual{i})
                    modes(i) = ci.AttributeValueBinsManual;
                elseif ~isnan(ci.attributeValueBinsAutoCount(i))
                    modes(i) = ci.attributeValueBinsAutoMode(i);
                else
                    modes(i) = ci.AttributeValueListAuto;
                end
            end
        end

        % determine the number of attributes, where possible, otherwise
        % leave as NaN. returns A x 1 numeric array
        function nv = get.nValuesByAttribute(ci)
            modes = ci.attributeValueModes;
            nv = nanvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                switch modes(i) 
                    case {ci.AttributeValueListManual, ci.AttributeValueBinsManual}
                        nv(i) = numel(ci.attributeValueLists{i});
                    case {ci.AttributeValueBinsAutoQuantiles, ci.AttributeValueBinsAutoUniform}
                        nv(i) = ci.attributeValueBinsAutoCount(i);
                    otherwise
                        % will be automatically determined in conditioninfo
                        % not actually one value, but want to put '?' there
                        nv(i) = 1;
                end
            end
        end

        function desc = get.attributeDescriptions(ci)
            desc = cellvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                name = ci.attributeNames{i};  
                nValues = ci.nValuesByAttribute(i);
                nAutoBins = ci.attributeValueBinsAutoCount(i);

                switch ci.attributeValueModes(i)
                    case ci.AttributeValueListManual
                        suffix = sprintf('(%d)', nValues);
                    case ci.AttributeValueListAuto
                        suffix = '(? auto)';
                    case ci.AttributeValueBinsManual
                        suffix = sprintf('(%d bins)', nValues);
                    case ci.AttributeValueBinsAutoUniform
                        suffix = sprintf('(%d unif bins)', nAutoBins);
                    case ci.AttributeValueBinsAutoQuantiles
                        suffix = sprintf('(%d quantiles)', nAutoBins);
                end

                desc{i} = sprintf('%s %s', name, suffix);
            end
        end 

        % add a new attribute
        function ci = addAttribute(ci, name, varargin)
            ci.warnIfNoArgOut(nargout);

            p = inputParser;
            p.addRequired('name', @ischar);
            % is this attribute always numeric?
            % list of allowed values for this value (other values will be ignored)
            p.addParamValue('requestAs', '', @ischar);
            p.addParamValue('valueList', {}, @(x) isnumeric(x) || iscell(x)); 
            p.addParamValue('valueBins', {}, @(x) isnumeric(x) || iscell(x));
            p.parse(name, varargin{:});
            valueList = p.Results.valueList;
            requestAs = p.Results.requestAs;
            if isempty(requestAs)
                requestAs = name;
            end

            if ci.hasAttribute(name)
                error('ConditionDescriptor already has attribute %s', name);
            end

            iAttr = ci.nAttributes + 1;
            ci.attributeNames{iAttr} = name;
            ci.attributeNumeric(iAttr) = isnumeric(valueList) || islogical(valueList); 
            ci.attributeRequestAs{iAttr} = requestAs;

            if isempty(valueList)
                ci.attributeValueListsManual{iAttr} = {};
            else
                if ~iscell(valueList)
                    valueList = num2cell(valueList);
                end
                ci.attributeValueListsManual{iAttr} = valueList;
            end

            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeValueBinsAutoCount(iAttr) = NaN;
            ci.attributeValueBinsAutoMode(iAttr) = NaN;

            ci = ci.invalidateCache();
        end

        % remove an existing attribute
        function ci = removeAttribute(ci, varargin)
            ci.warnIfNoArgOut(nargout);

            if iscell(varargin{1})
                attributes = varargin{1};
            else
                attributes = varargin;
            end

            ci.warnIfNoArgOut(nargout);

            if ~isnumeric(attributes)
                % check all exist
                ci.getAttributeIdx(attributes);
            else
                attributes = ci.attributeNames(attributes);
            end

            if ~ci.hasAttribute(attributes)
                error('ConditionDescriptor has no attribute %s', name);
            end

            iAttr = ci.getAttributeIdx(attributes);
            maskOther = true(ci.nAttributes, 1);
            maskOther(iAttr) = false;

            ci = ci.maskAttributes(maskOther);
        end
        
        function ci = maskAttributes(ci, mask)
            ci.warnIfNoArgOut(nargout);
            
            idxRemove = find(~mask);
            if ~any(idxRemove)
                return;
            end
            
            % first remove the attributes from any axes they are on
            ci = ci.removeAttributeFromAxes(idxRemove); 
            
            % then remove it from the attribute lists
            ci.attributeNames = ci.attributeNames(mask);
            ci.attributeRequestAs = ci.attributeRequestAs(mask);
            ci.attributeNumeric = ci.attributeNumeric(mask);
            ci.attributeValueList = ci.attributeValueList(mask);
            ci.attributeValueBinsAutoCount = ci.attributeValueBinsAutoCount(mask);
            ci.attributeValueBinsAutoMode = ci.attributeValueBinsAutoMode(mask);
            ci.attributeValueBinsManual = ci.attributeValueBinsManual(mask);
        end 

        % manually set the attribute value list
        function ci = setAttributeValueList(ci, name, valueList)
            ci.warnIfNoArgOut(nargout);

            iAttr = ci.getAttributeIdx(name);           
            if isempty(valueList)
                ci.attributeValueListManual{iAttr} = {};
            else
                ci.attributeValueListManual{iAttr} = valueList;                
            end
            ci.attributeNumeric(iAttr) = isnumeric(valueList) || islogical(valueList);

            ci = ci.invalidateCache();
        end

        % manually set attribute bins
        function ci = binAttribute(ci, name, bins)
            ci.warnIfNoArgOut(nargout);

            if isvector(bins)
                assert(issorted(bins), 'Bins specified as vector must be in sorted order');
                binsMat = nan(numel(bins)-1, 2);
                binsMat(:, 1) = bins(1:end-1);
                binsMat(:, 2) = bins(2:end);
            else
                assert(ismatrix(bins) && size(bins, 2) == 2, 'Bins matrix must be nBins x 2');
                assert(all(bins(:, 2) >= bins(:, 1)), 'Bins matrix must have larger value in second column than first');

                binsMat = bins;
            end 

            iAttr = ci.getAttributeIdx(name); 
            ci.attributeValueBinsManual{iAttr} = binsMat;
            ci.attributeNumeric(iAttr) = true;
            ci.attributeValueListManual{iAttr} = {};
            ci.attributeValueBinsAutoCount = NaN;
            ci.attributeValueBinsAutoMode = NaN;

            ci = ci.invalidateCache();
        end

        % automatically set attribute binned uniformly by range
        function ci = binAttributeUniform(ci, name, nBins)
            ci.warnIfNoArgOut(nargout);
            
            iAttr = ci.getAttributeIdx(name);

            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeNumeric(iAttr) = true;
            ci.attributeValueListManual{iAttr} = {};
            ci.attributeValueBinsAutoCount = nBins;
            ci.attributeValueBinsAutoMode = ci.AttributeValueBinsAutoUniform;

            ci = ci.invalidateCache();
        end

        % automatically set attribute binned into quantiles
        function ci = binAttributeQuantiles(ci, name, nQuantiles)
            ci.warnIfNoArgOut(nargout);

            iAttr = ci.assertHasAttribute(name);
            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeNumeric(iAttr) = true;
            ci.attributeValueListsManual{iAttr} = {};
            ci.attributeValueBinsAutoCount(iAttr) = nQuantiles;
            ci.attributeValueBinsAutoMode(iAttr) = ci.AttributeValueBinsAutoQuantiles;

            ci = ci.invalidateCache();
        end
    end

    methods(Static) % Default nameFn and appearanceFn
        function name = defaultNameFn(ci, attrValues, conditionSubs) %#ok<INUSD>
            name = '';
            attr = fieldnames(attrValues);
            attrIsNumeric = ci.getIsAttributeNumeric(attr);

            for iAttr = 1:length(attr)
                include = false;
                val = attrValues.(attr{iAttr});

                if ~ischar(val) && ~isscalar(val)
                    % skip attributes where more than one value is
                    % specified, they won't be part of the condition name
                    continue;
                end

                if isnumeric(val)
                    if isscalar(val)
                        val = num2str(val);
                    else
                        val = mat2str(val);
                    end
                    include = true;
                elseif ischar(val)
                    % okay as is
                    include = true;
                elseif length(val) > 1
                    include = false;
                end

                % include attribute name if its numeric
                if attrIsNumeric(iAttr)
                    val = num2str(val);
                end

                if include
                    name = [name attr{iAttr} '=' val ' ']; %#ok<AGROW>
                end
            end

            name = strtrim(name);
        end

        function a = defaultAppearanceFn(ci, varargin)
            % returns a struct specifying the default set of appearance properties 
            % for the given group. indsGroup is a length(ci.groupByList) x 1 array
            % of the inds where this group is located in the high-d array, and dimsGroup
            % gives the full dimensions of the list of groups.
            %
            % We vary color along all axes simultaneously, using the linear
            % inds.
            %
            % Alternatively, if no arguments are passed, simply return a set of defaults

            nConditions = ci.nConditions;

            a = emptyStructArray(ci.conditionsSize, {'color', 'lineWidth'});

            if nConditions == 1
                cmap = [0.3 0.3 1];
            else
                %cmap = jet(nConditions);
                cmap =pmkmp(nConditions, 'isol');
            end

            for iC = 1:nConditions
                a(iC).lineWidth = 2;
                a(iC).color = cmap(iC, :);
            end
        end

        function a = frozenAppearanceFn(ci, a, varargin)
            % this function looks at ci.frozenAppearanceConditions and
            % .frozenAppearanceData and does a lookup of the stored
            % appearance for each condition, essentially allowing you to
            % freeze the condition appearance through filtering,
            % regrouping, etc.
            %
            % Call .freezeAppearance() to activate

            % for each condition in ci, search
            % ci.frozenAppearanceConditions for the first match

            matchIdx = nan(ci.nConditions, 1);
            fieldsCurrent = fieldnames(ci.conditions);
            if ~isempty(ci.frozenAppearanceConditions)
                fieldsFrozen = fieldnames(ci.frozenAppearanceConditions);
            else
                fieldsFrozen = {};
            end
            fieldsCheck = intersect(fieldsCurrent, fieldsFrozen);
            nFrozenConditions = numel(ci.frozenAppearanceConditions);

            for iC = 1:ci.nConditions
                for iCFrozen = 1:nFrozenConditions
                    isMatch = true;
                    for iF = 1:numel(fieldsCheck)
                        fld = fieldsCheck{iF};
                        if ~isequal(ci.conditions(iC).(fld), ci.frozenAppearanceConditions(iCFrozen).(fld))
                            isMatch = false;
                            break;
                        end                        
                    end

                    if isMatch
                        matchIdx(iC) = iCFrozen;
                        break;
                    end
                end
            end

            mask = ~isnan(matchIdx);
            a(mask) = ci.frozenAppearanceData(matchIdx(mask));
        end
    end

    methods(Static) % construct from another condition descriptor, used primarily by ConditionInfo
        function cdNew = fromConditionDescriptor(cd, cdNew)
            cd.warnIfNoArgOut(nargout);
            
            if nargin < 2
                cdNew = ConditionDescriptor();
            end

            meta = ?ConditionDescriptor;
            props = meta.PropertyList;

            for iProp = 1:length(props)
                prop = props(iProp);
                if prop.Dependent || prop.Constant || prop.Transient
                    continue;
                else
                    name = prop.Name;
                    cdNew.(name) = cd.(name);
                end
            end

            cdNew = cdNew.invalidateCache();
        end
    end

    methods % copy if handle
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
            end
        end
    end
end

