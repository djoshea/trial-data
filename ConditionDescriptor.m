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
    
    properties(SetAccess=protected)
        attributeNumeric = []; % A x 1 logical array : is this attribute a numeric value? 
        attributeValueListsManual = {}; % A x 1 cell array of permitted values (or cells of values) for this attribute
        attributeValueBinsManual = {}; % A x 1 cell array of value Nbins x 2 value bins to use for numeric lists
        attributeValueBinsAutoCount % A x 1 numeric array of Nbins to use when auto computing the bins, NaN if not in use
        attributeValueBinsAutoModes % A x 1 numeric array of either AttributeValueBinsAutoUniform or AttributeValueBinsAutoQuantiles
        
        axisValueListsManual % G x 1 cell of cells: each contains a struct specifying an attribute specification for each element along the axis
        axisValueListsOccupiedOnly % G x 1 logical indicating whether to constrain the combinatorial valueList to only occupied elements (with > 0 trials)
        axisRandomizeModes % G x 1 numeric of constants beginning with Axis* (see below)
    end
    
    properties(Hidden, SetAccess=protected)
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
        
        % these are X-dimensional objects where X is nAxes
        conditions % X-dimensional struct where values(...idx...).attribute is the value of that attribute on that condition
        conditionsAsStrings % includes attribute values as strings rather than numeric 
        conditionsAxisAttributesOnly % includes only the attributes actively selected for
        
        appearances % A-dimensional struct of appearance values
        names % A-dimensional cellstr array with names of each condition 
        attributeValueLists % A x 1 cell array of values allowed for this attribute
                           % here just computed from attributeValueListManual, but in ConditionInfo
                           % can be automatically computed from the data
        attributeValueListsAsStrings % same as above, but everything is a string
        
        axisValueLists % G dimensional cell array of structs which select attribute values for that position along an axis
        axisValueListsAsStrings
        axisValueListModes % G dimensional array of 
    end
    
    % how are attribute values determined for a given attribute?
    properties(Constant, Hidden)
        % for attributeValueListModes
        AttributeValueListManual = 1;
        AttributeValueListAuto = 2;
        AttributeValueBinsManual = 3;
        AttributeValueBinsAutoUniform = 4;
        AttributeValueBinsAutoQuantiles = 5;
        
        % for axisRandomizeModes
        AxisOriginal = 1;
        AxisShuffled = 2;
        AxisResampled = 3;
        AxisResampledFromFirst = 4;
        
        % for axisValueListModes
        AxisValueListAutoAll = 1;
        AxisValueListAutoOccupied = 2;
        AxisValueListManual = 3;
    end
        
    properties(Dependent, Transient)
        nAttributes % how many attributes: ndims(values)
        attributeDescriptions
        nValuesByAttribute % how many values per attribute: size(values)
        attributeAlongWhichAxis % A x 1 array indicating which axis an attribute contributes to (or NaN)
        attributeValueModes % A x 1 array of AttributeValue* constants above
        attributeActsAsFilter % A x 1 logical array : does this attribute have a
                % value list or manual bin setup that would invalidate trials?
        
        nAxes % how many dimensions of grouping axe
        nValuesAlongAxes % X x 1 array of number of elements along the axis
        axisDescriptions % strcell describing each axis
        
        nConditions % how many total conditions
        conditionsSize 
        conditionsAsLinearInds % linear index corresponding to each condition if flattened 
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
        
        function printDescription(ci) 
            tcprintf('yellow', '%s:\n', class(ci));
            tcprintf('inline', '\t{bright blue}Attributes: {white}%s\n', strjoin(ci.attributeDescriptions));
            tcprintf('inline', '\t{bright blue}Axes: {white}%s\n', strjoin(ci.axisDescriptions, ' , '));
        end
        
        function printOneLineDescription(ci)           
            if ci.nAxes == 0
                axisStr = 'no grouping axes';
            else
                axisStr = strjoin(ci.axisDescriptions, ' , ');
            end
            
            attrFilter = ci.attributeNames(ci.attributeActsAsFilter);
            if isempty(attrFilter)
                filterStr = 'no filtering';
            else
                filterStr = sprintf('filtering by %s', strjoin(attrFilter));
            end
            
            tcprintf('inline', '{yellow}%s: {none}%s, %s\n', ...
                class(ci), axisStr, filterStr);
        end

        function disp(ci)
            ci.printDescription();
            fprintf('\n');
            builtin('disp', ci);
        end
    end

    methods % Axis related 
        function n = get.nAxes(ci)
            n = numel(ci.axisAttributes);
        end
        
        function a = get.attributeAlongWhichAxis(ci)
            a = nanvec(ci.nAttributes);
            for iX = 1:ci.nAxes
                a(ci.getAttributeIdx(ci.axisAttributes{iX})) = iX;
            end
        end
        
        function modes = get.axisValueListModes(ci)
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
        
        % determine whether each attribute acts to filter valid trials
        function tf = get.attributeActsAsFilter(ci)
            modes = ci.attributeValueModes;
            tf = ismember(modes, [ci.AttributeValueListManual, ci.AttributeValueBinsManual]);
        end
        
        function desc = get.axisDescriptions(ci)
            desc = cellvec(ci.nAxes);
            
            for iX = 1:ci.nAxes
                attr = ci.axisAttributes{iX};
                nv = ci.conditionsSize(iX);
                switch ci.axisValueListModes(iX)
                    case ci.AxisValueListAutoAll
                        vlStr = ' auto';
                    case ci.AxisValueListAutoOccupied
                        vlStr = ' autoOccupied';
                    case ci.AxisValueListManual
                        vlStr = ' manual';
                end
                        
                switch ci.axisRandomizeModes(iX)
                    case ci.AxisOriginal
                        randStr = '';
                    case ci.AxisShuffled
                        randStr = ' shuffled';
                    case ci.AxisResampled
                        randStr = ' resampled';
                    case ci.AxisResampledFromFirst
                        randStr = ' resampledFromFirst';
                end
                
                desc{iX} = sprintf('%s%s (%d%s)', ...
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
            ci.axisRandomizeModes(idx) = ci.AxisOriginal;
            ci.axisValueListsOccupiedOnly(idx) = false;

            ci = ci.invalidateCache();
        end
        
        function ci = maskAxes(ci, mask)
            ci.warnIfNoArgOut(nargout);
            
            ci.axisAttributes = ci.axisAttributes(mask);
            ci.axisValueListsManual = ci.axisValueListsManual(mask);
            ci.axisRandomizeModes = ci.axisRandomizeModes(mask);
            ci.axisValueListsOccupiedOnly = ci.axisValueListsOccupiedOnly(mask);
            
            ci = ci.invalidateCache();
        end

        % wipe out existing axes and creates simple auto axes along each 
        function ci = groupBy(ci, varargin)
            ci.warnIfNoArgOut(nargout);
            ci = ci.clearAxes();
            
            for i = 1:numel(varargin)
                ci = ci.addAxis(varargin{i});
            end
        end

        function ci = groupByAll(ci)
            ci.warnIfNoArgOut(nargout);
            ci = ci.groupBy(ci.attributeNames{:});
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
                    ci.axisRandomizeModes(iX) = ci.AxisOriginal;
                end
            end
            
            ci = ci.maskAxes(~removeAxisMask);
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
    end

    methods % Attribute related 
        function [tf, idx] = hasAttribute(ci, name)
            if isnumeric(name)
                [tf, idx] = ismember(name, 1:ci.nAttributes);
            else
                [tf, idx] = ismember(name, ci.attributeNames);
            end
        end

        function idx = assertHasAttribute(ci, name)
            [tf, idx] = ci.hasAttribute(name);
            if ~all(tf)
                if isnumeric(name)
                    name = strjoin(name(~tf), ', ');
                elseif iscell(name)
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
                    modes(i) = ci.attributeValueBinsAutoModes(i);
                else
                    modes(i) = ci.AttributeValueListAuto;
                end
            end
        end

        % determine the number of attributes, where possible, otherwise
        % leave as NaN. returns A x 1 numeric array
        function nv = get.nValuesByAttribute(ci)
            nv = nanvec(ci.nAttributes);
            for i = 1:ci.nAttributes
                nv(i) = numel(ci.attributeValueLists{i});
            end
        end

        function desc = get.attributeDescriptions(ci)
            desc = cellvec(ci.nAttributes);
            isFilter = ci.attributeActsAsFilter;
            modes = ci.attributeValueModes;
            for i = 1:ci.nAttributes
                name = ci.attributeNames{i};  
                nValues = ci.nValuesByAttribute(i);
                nAutoBins = ci.attributeValueBinsAutoCount(i);

                switch modes(i)
                    case ci.AttributeValueListManual
                        suffix = sprintf('(%d)', nValues);
                    case ci.AttributeValueListAuto
                        suffix = sprintf('(%d auto)', nValues);
                    case ci.AttributeValueBinsManual
                        suffix = sprintf('(%d bins)', nValues);
                    case ci.AttributeValueBinsAutoUniform
                        suffix = sprintf('(%d bins)', nAutoBins);
                    case ci.AttributeValueBinsAutoQuantiles
                        suffix = sprintf('(%d quantiles)', nAutoBins);
                end
                
                if isFilter(i)
                    filterStr = ' [filter]';
                else
                    filterStr = '';
                end

                if ci.attributeNumeric(i)
                    numericStr = '#';
                else
                    numericStr = '';
                end
                desc{i} = sprintf('%s %s%s%s', name, numericStr, suffix, filterStr);
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
                assert(isnumeric(valueList) || iscell(valueList), 'ValueList must be numeric or cell');
                % filter for unique values or 
                ci.attributeValueListsManual{iAttr} = unique(valueList, 'stable');
            %    if ~iscell(valueList)
             %       valueList = num2cell(valueList);
             %   end
            end

            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeValueBinsAutoCount(iAttr) = NaN;
            ci.attributeValueBinsAutoModes(iAttr) = NaN;

            ci = ci.invalidateCache();
        end

        function ci = addAttributes(ci, names)
            ci.warnIfNoArgOut(nargout);
            for i = 1:numel(names)
                ci = ci.addAttribute(names{i});
            end
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
            ci = ci.removeAttributesFromAxes(idxRemove); 

            % then remove it from the attribute lists
            ci.attributeNames = ci.attributeNames(mask);
            ci.attributeRequestAs = ci.attributeRequestAs(mask);
            ci.attributeNumeric = ci.attributeNumeric(mask);
            ci.attributeValueListsManual = ci.attributeValueListsManual(mask);
            ci.attributeValueLists = ci.attributeValueLists(mask);
            ci.attributeValueListsAsStrings = ci.attributeValueListsAsStrings(mask);
            ci.attributeValueBinsAutoCount = ci.attributeValueBinsAutoCount(mask);
            ci.attributeValueBinsAutoModes = ci.attributeValueBinsAutoModes(mask);
            ci.attributeValueBinsManual = ci.attributeValueBinsManual(mask);
        end 
        
        % set all attribute value lists to auto
        function ci = setAllAttributeValueListsAuto(ci)
            ci.warnIfNoArgOut(nargout);
            for i = 1:ci.nAttributes
                ci = ci.setAttributeValueListAuto(i);
            end
        end
        
        % restore value list to automatically include all values, with no
        % binning
        function ci = setAttributeValueListAuto(ci, attr)
            ci.warnIfNoArgOut(nargout);
            iAttr = ci.assertHasAttribute(attr);
            ci.attributeValueListsManual{iAttr} = [];
            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeValueBinsAutoCount(iAttr) = NaN;
            ci.attributeValueBinsAutoModes(iAttr) = NaN;
            ci = ci.invalidateCache();
        end
        
        function ci = setAttributeNumeric(ci, attr, tf)
            ci.warnIfNoArgOut(nargout);
            iAttr = ci.assertHasAttribute(attr);
            ci.attributeNumeric(iAttr) = tf;
        end  

        % manually set the attribute value list
        function ci = setAttributeValueList(ci, name, valueList)
            ci.warnIfNoArgOut(nargout);

            iAttr = ci.getAttributeIdx(name);           
            if isempty(valueList)
                ci.attributeValueListsManual{iAttr} = {};
            else
                ci.attributeValueListsManual{iAttr} = valueList;                
            end
            
            %ci.attributeNumeric(iAttr) = isnumeric(valueList) || islogical(valueList);

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
            ci.attributeValueBinsAutoCount(iAttr) = NaN;
            ci.attributeValueBinsAutoModes(iAttr) = NaN;

            ci = ci.invalidateCache();
        end

        % automatically set attribute binned uniformly by range
        function ci = binAttributeUniform(ci, name, nBins)
            ci.warnIfNoArgOut(nargout);
            
            iAttr = ci.getAttributeIdx(name);

            ci.attributeValueBinsManual{iAttr} = [];
            ci.attributeNumeric(iAttr) = true;
            ci.attributeValueListsManual{iAttr} = {};
            ci.attributeValueBinsAutoCount(iAttr) = nBins;
            ci.attributeValueBinsAutoModes(iAttr) = ci.AttributeValueBinsAutoUniform;

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
            ci.attributeValueBinsAutoModes(iAttr) = ci.AttributeValueBinsAutoQuantiles;

            ci = ci.invalidateCache();
        end
    end

    % get, set data stored inside odc
    methods 
        % NOTE: all of these should copy odc before writing to it
        
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
        
        function v = get.conditionsAsStrings(ci)
            v = ci.odc.conditionsAsStrings;
            if isempty(v)
                ci.odc.conditionsAsStrings = ci.buildConditionsAsStrings();
                v = ci.odc.conditionsAsStrings;
            end
        end
        
        function ci = set.conditionsAsStrings(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionsAsStrings = v;
        end
        
        function v = get.conditionsAxisAttributesOnly(ci)
            v = ci.odc.conditionsAxisAttributesOnly;
            if isempty(v)
                ci.odc.conditionsAxisAttributesOnly = ci.buildConditionsAxisAttributesOnly();
                v = ci.odc.conditionsAxisAttributesOnly;
            end
        end
        
        function ci = set.conditionsAxisAttributesOnly(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.conditionsAxisAttributesOnly = v;
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
            v = ci.odc.attributeValueLists;
            if isempty(v)
                ci.odc.attributeValueLists = ci.buildAttributeValueLists();
                v = ci.odc.attributeValueLists;
            end
        end
        
        function ci = set.attributeValueLists(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.attributeValueLists = v;
        end
        
        function v = get.attributeValueListsAsStrings(ci)
            v = ci.odc.attributeValueListsAsStrings;
            if isempty(v)
                ci.odc.attributeValueListsAsStrings = ci.buildAttributeValueListsAsStrings();
                v = ci.odc.attributeValueListsAsStrings;
            end
        end
        
        function ci = set.attributeValueListsAsStrings(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.attributeValueListsAsStrings = v;
        end

        function v = get.axisValueLists(ci)
            v = ci.odc.axisValueLists;
            if isempty(v)
                ci.odc.axisValueLists = ci.buildAxisValueLists();
                v = ci.odc.axisValueLists;
            end
        end

        function ci = set.axisValueLists(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.axisValueLists = v;
        end
        
        function v = get.axisValueListsAsStrings(ci)
            v = ci.odc.axisValueListsAsStrings;
            if isempty(v)
                ci.odc.axisValueListsAsStrings = ci.buildAxisValueListsAsStrings();
                v = ci.odc.axisValueListsAsStrings;
            end
        end

        function ci = set.axisValueListsAsStrings(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.axisValueListsAsStrings = v;
        end
    end

    % build data stored inside odc (used by getters above)
    methods 
        function values = buildConditionsAxisAttributesOnly(ci)
            if ci.nAxes == 0
                values = struct();
            else
                valueLists = ci.axisValueLists; 
                values = TensorUtils.mapFromAxisLists(@structMergeMultiple,...
                    valueLists, 'asCell', false);
            end
        end
        
        function values = buildConditions(ci)
            values = ci.conditionsAxisAttributesOnly;
            
            % and add "wildcard" match for all other attributes that act as
            % filter (i.e. have manual value list or bins specified)
            whichAxis = ci.attributeAlongWhichAxis;
            isFilter = ci.attributeActsAsFilter;
            valueLists = ci.attributeValueLists;
            for iA = 1:ci.nAttributes
                if isnan(whichAxis(iA)) && isFilter(iA)
                    valueList = valueLists{iA};
                    
                    % flatten any subgroupings of values in the value list
                    if ci.attributeNumeric(iA) && iscell(valueList)
                        valueList = [valueList{:}];
                    elseif ~ci.attributeNumeric(iA) && ~iscellstr(valueList) && ~ischar(valueList)
                        valueList = [valueList{:}];
                    end
                    % wrap in cell to avoid scalar expansion
                    values = assignIntoStructArray(values, ci.attributeNames{iA}, {valueList});
                end
            end
        end
        
        % build a wildcard search struct where each .attribute field is the
        % value list for that attribute
        function values = buildStructAllAttributeValueLists(ci)
            values = struct();
            for iA = 1:ci.nAttributes
                values = assignIntoStructArray(values, ci.attributeNames{iA}, ...
                    ci.attributeValueLists(iA));
            end
        end
        
        function values = buildStructNonAxisAttributeValueLists(ci)
            whichAxis = ci.attributeAlongWhichAxis;
            values = struct();
            for iA = 1:ci.nAttributes
                if isnan(whichAxis(iA))
                    values = assignIntoStructArray(values, ci.attributeNames{iA}, ...
                        ci.attributeValueLists(iA));
                end
            end
        end
        
        function values = buildConditionsAsStrings(ci)
            if ci.nAxes == 0
                values = {structToString(ci.conditions)};
            else
                valueLists = ci.axisValueListsAsStrings; 
                values = TensorUtils.mapFromAxisLists(@(varargin) strjoin(varargin, ' '),...
                    valueLists, 'asCell', true);
            end
        end
        
        function valueListByAxes = buildAxisValueLists(ci)
            valueListByAxes = cellvec(ci.nAxes);
            for iX = 1:ci.nAxes
                % build a cellstr of descriptions of the values along this axis
               
                % G x 1 cell of cells: each contains a struct specifying an attribute specification for each element along the axis
                if isempty(ci.axisValueListsManual{iX})
                    % build auto list of attributes
                    valueListByAxes{iX} = makecol(buildAutoValueListForAttributeSet(ci.axisAttributes{iX}));
                else
                    valueListByAxes{iX} = makecol(ci.axisValueListsManual{iX});
                end
            end

            function values = buildAutoValueListForAttributeSet(attributes)
                % build a struct array for a set of attributes that walks all possible combinations of the attribute values 
                if ischar(attributes)
                    attributes = {attributes};
                end
                attrIdx = ci.getAttributeIdx(attributes);
                valueLists = ci.attributeValueLists(attrIdx);

                % convert bin edges value lists to the cell vectors
                for i = 1:numel(attrIdx)
                    switch ci.attributeValueModes(attrIdx(i))
                        case {ci.AttributeValueBinsManual, ci.AttributeValueBinsAutoUniform, ...
                                ci.AttributeValueBinsAutoQuantiles}
                            % convert valueList from Nbins x 2 matrix to
                            % Nbins x 1 cellvector so that it gets mapped
                            % correctly
                            valueLists{i} = mat2cell(valueLists{i}, ones(size(valueLists{i}, 1), 1));
                    end
                end
                            
                values = TensorUtils.mapFromAxisLists(@buildStruct, valueLists, ...
                    'asCell', false);

                function s = buildStruct(varargin)
                    for j = 1:numel(varargin)
                        s.(attributes{j}) = varargin{j};
                    end
                end

            end
        end
        
        function strCell = buildAxisValueListsAsStrings(ci)
            strCell = cellvec(ci.nAxes);
            valueLists = ci.axisValueLists;
            
            % describe the list of values selected for along each position on each axis
            for iX = 1:ci.nAxes  
                
                % start with axisValueLists
                attr = ci.axisAttributes{iX};
                attrIdx = ci.getAttributeIdx(attr);
                
                % replace binned values with strings
                for iA = 1:numel(attrIdx)
                    switch ci.attributeValueModes(attrIdx(iA))
                        case {ci.AttributeValueBinsManual, ci.AttributeValueBinsAutoUniform, ...
                                ci.AttributeValueBinsAutoQuantiles}
                            % convert valueList from 1 x 2 vector to '#-#' string
                            for iV = 1:numel(valueLists{iX})
                                valueLists{iX}(iV).(attr{iA}) = sprintf('%g-%g', valueLists{iX}(iV).(attr{iA}));
                            end
                    end
                end
                
                strCell{iX} = arrayfun(@structToString, makecol(valueLists{iX}), ...
                   'UniformOutput', false);
            end
        end

        function names = buildNames(ci)
            % pass along values(i) and the subscripts of that condition in case useful 
            if ci.nConditions > 0
                fn = ci.nameFn;
                if isempty(fn)
                    fn = @ConditionDescriptor.defaultNameFn;
                end
                names = fn(ci);
                assert(iscellstr(names) && isequal(size(names), ci.conditionsSize), ...
                    'nameFn must return cellstr with same size as .conditions');
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
                    case ci.AttributeValueBinsAutoUniform
                        % placeholder string to be replaced by actual bins
                        % matrix
                        valueList{i} = arrayfun(@(bin) sprintf('bin%d', bin), ...
                            1:ci.attributeValueBinsAutoCount(i), 'UniformOutput', false);
                    case ci.AttributeValueBinsAutoQuantiles
                        % the number of bins is known, so they can be specified here
                        valueList{i} = arrayfun(@(bin) sprintf('quantile%d', bin), ...
                            1:ci.attributeValueBinsAutoCount(i), 'UniformOutput', false);
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
        
        function valueList = buildAttributeValueListsAsStrings(ci)
            modes = ci.attributeValueModes;
            valueList = ci.attributeValueLists;
            for i = 1:ci.nAttributes
                switch modes(i) 
                    case ci.AttributeValueListManual
                        if ci.attributeNumeric(i)
                            if iscell(valueList{i})
                                % could have multiple attribute values
                                % grouped together as one element
                                valueList{i} = cellfun(@(i) sprintf('%.3g', i), valueList{i}, 'UniformOutput', false);
                                valueList{i} = cellfun(@(vals) strjoin(vals, ','), valueList{i}, 'UniformOutput', false);
                            else
                                valueList{i} = arrayfun(@(i) sprintf('%.3g', i), valueList{i}, 'UniformOutput', false);
                            end
                        else
                            % non-numeric, can leave as is unless...
                            if ~iscellstr(valueList{i})
                                % could have multiple attribute values
                                % grouped together as one element
                                valueList{i} = cellfun(@(vals) strjoin(vals, ','), valueList{i}, 'UniformOutput', false);
                            end
                        end
                                
                    case {ci.AttributeValueBinsManual, ci.AttributeValueBinsAutoUniform, ci.AttributeValueBinsAutoQuantiles}
                        if ~iscell(valueList{i})
                            bins = valueList{i};
                            valueList{i} = arrayfun(@(row) sprintf('%g-%g', bins(row, 1), bins(row, 2)), ...
                                1:size(bins, 1), 'UniformOutput', false);
                        else
                            % already cellstr for auto bins, leave as is
                        end
                    case ci.AttributeValueListAuto
                        % auto list leave empty, must be determined when
                        % ConditionInfo applies it to data
                        valueList{i} = {'?'};   
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
    
    methods(Static) % Default nameFn and appearanceFn
        function nameCell = defaultNameFn(ci, varargin) 
            % receives the condition descriptor itself and returns a
            %  a cell tensor specifying the names of each condition
            
            nameCell = ci.conditionsAsStrings;
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
                if nConditions > 256
                    cmap = jet(nConditions);
                else
                    cmap =pmkmp(nConditions, 'isol');
                end
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
        
        % construct condition descriptor from a struct of attribute values
        % for numeric attributes, if there are more than 10 different
        % values, the attribute will be binned into quintiles
        function cd = fromStruct(s)
            cd = ConditionDescriptor();
            cd = cd.addAttributes(fieldnames(s));
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                warning('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
            end
        end
        
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy(); %#ok<MCNPN>
            end
        end
    end
end

