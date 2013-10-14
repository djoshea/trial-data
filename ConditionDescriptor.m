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
        attributeNumeric = []; % A x 1 logical array : is this attribute a numeric value? 
    end
    
    properties(SetAccess=protected)
        % When re-arranging the axes, default condition appearances can get shuffled around
        % which can make comparison across figures difficult. These cache
        % the condition appearances to make things easier. Call
        % freezeAppearances to activate()
        appearanceFrozen = false;
        frozenAppearanceConditions
        frozenAppearanceData

        % don't access these directly, call .getAttributeValueList
        attributeValueListManual = {}; % A x 1 cell array of permitted values for this attribute
        attributeValueListSpecified = []; % A x 1 logical array: was a non-empty value list provided for this attribute

        axisSpec % specification for how to form the group axes
    end
    
    % END OF STORED TO DISK PROPERTIES
    
    properties(Transient, Hidden)
        noUpdateCache = false; % if true, invalidateCache does nothing, useful for loading or building 
                               % be sure to set to true and call invalidateCache afterwards                        
    end
    
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
        appearances % A-dimensional struct of appearance values
        names % A-dimensional cellstr array with names of each condition 
        attributeValueList % A x 1 cell array of values allowed for this attribute
                           % here just computed from attributeValueListManual, but in ConditionInfo
                           % can be automatically computed from the data
    end

    properties(Dependent, Transient)
        nAttributes % how many attributes: ndims(values)

        nValuesByAttribute % how many values per attribute: size(values)

        nConditions % how many total conditions
        
        %groupByListAttributeIdx % idx into .attributeNames of each attribute in groupByList
        
        %isAttributeInGroupByList % mask indicating which .attributeNames{i} is in groupByList

        %attributeValueListGroupBy % same as attributeValueList but for grouped attributes only
        
        %nAttributesGroupBy % how many attributes in group by list
        
        conditionsSize 

        conditionsAsLinearInds % linear index corresponding to each condition if flattened 
        
        %nValuesByAttributeGroupBy % same as conditionsSize (except won't auto expand to be N x 1)
        
        %attributeNamesGroupBy % attributeNames(isAttributeInGroupByList)
        
        %attributeRequestAsGroupBy % attributeRequestAs(isAttributeInGroupByList)
        
        %conditionsWithGroupByFieldsOnly % same as conditions, but with only fields in groupByList specified
    end
    
    % Constructor, load, save methods
    methods
        function ci = ConditionDescriptor()
            ci.odc = ci.buildOdc();
        end
        
        function odc = buildOdc(ci)
            odc = ConditionDescriptorOnDemandCache();
        end
    end
   
    % Specify and manipulate group axes
    methods
        function ci = addAxis(ci, varargin)
            ci.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addParamValue('name', '', @ischar);
            p.addParamValue('attributes', {}, @iscellstr);
            p.addParamValue('values', [], @(x) true);
            p.parse(varargin{:});
            
            ax.attributes = p.Results.attributes;
            ci.assertHasAttribute(ax.attributes);
            
            
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
        
        function v = get.attributeValueList(ci)
            v = ci.odc.attributeValueList;
            if isempty(v)
                ci.odc.attributeValueList = ci.buildAttributeValueList();
                v = ci.odc.attributeValueList;
            end
        end
        
        function ci = set.attributeValueList(ci, v)
            ci.odc = ci.odc.copy();
            ci.odc.attributeValueList = v;
        end
    end

    % build data stored inside odc
    methods 
        function values = buildConditions(ci)
            if ci.nConditions > 0
                nAttr = ci.nAttributes;
                nAttrGroupBy = ci.nAttributesGroupBy;
                if nAttrGroupBy > 0
                    nValuesByAttr = ci.nValuesByAttributeGroupBy;
                    values = TensorUtils.mapToSizeFromSubs(nValuesByAttr, @getInfoFn);
                else
                    values = getInfoFn();
                end
            else
                values = struct([]);
            end
            
            function val = getInfoFn(varargin)
                % varargin{i} is the index into the value list of groupByAttribute i
                % return a struct where .attribute = value for both
                % groupByAttributes and non-groupByAttributes with only one
                % value allowed
                
                val = struct();
                % loop over group by attributes
                for iA = 1:nAttrGroupBy
                    val.(ci.groupByList{iA}) = ci.attributeValueListGroupBy{iA}{varargin{iA}};
                end
                
                % loop over everything else
                for iA = 1:nAttr
                    attr = ci.attributeNames{iA};
                    if ~ismember(attr, ci.groupByList)
                        list = ci.attributeValueList{iA};                        
                        if length(list) == 1
                            val.(attr) = list{1};
                        else
                            val.(attr) = list;
                        end
                    end
                end
            end
        end

        function names = buildNames(ci)
            % pass along values(i) and the subscripts of that condition in case useful 
            if ci.nConditions > 0
                nameFn = ci.nameFn;
                if isempty(nameFn)
                    nameFn = @ConditionDescriptor.defaultNameFn;
                end

                wrapFn = @(varargin) nameFn(ci, varargin{:});
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
        
        function valueList = buildAttributeValueList(ci)
            % just pull the manual lists
            valueList = ci.attributeValueListManual;
        end
        
        function valueList = getAttributeValueList(ci, name)
            idx = ci.getAttributeIdx(name);
            valueList = makecol(ci.attributeValueList{idx});
        end

        function valueIdx = getAttributeValueIdx(ci, attr, value)
            [tf valueIdx] = ismember(value, ci.getAttributeValueList(attr));
            assert(tf, 'Value not found in attribute %s valueList', attr);
        end
    end
    
    methods
        % flush the contents of odc as they are invalid
        % call this at the end of any methods which would want to
        % regenerate these values
        function ci = invalidateCache(ci)
            ci.warnIfNoArgOut(nargout);
            
            % here we precompute these things to save time, 
            % but each of these things also has a get method that will
            % recompute this for us
            ci.conditions = [];
            ci.appearances = [];
            ci.names = [];
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
            ci.names = [];
        end
        
        function ci = set.appearanceFn(ci, fn)
            ci.appearanceFn = fn;
            ci.appearances = [];
        end
        
        function tf = hasAttribute(ci, name)
            tf = ismember(name, ci.attributeNames);
        end
        
        function assertHasAttribute(ci, name)
            tf = ci.hasAttribute(name);
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

        function nv = get.nValuesByAttribute(ci)
            nv = cellfun(@length, ci.attributeValueList); 
        end

        function nv = get.conditionsSize(ci)
            if ci.nAttributes == 0
                nv = [1];
            elseif ci.nAttributesGroupBy == 0
                nv = [1 1];
            else
                nv = ci.nValuesByAttributeGroupBy; 
            end
            if isscalar(nv)
                nv(2) = 1;
            end
        end
        
        function linearInds = get.conditionsAsLinearInds(ci)
            linearInds = TensorUtils.containingLinearInds(ci.conditionsSize);
        end

        function n = get.nConditions(ci)
            n = prod(ci.conditionsSize);
        end

        function idxList = getAttributeIdx(ci,name)
            if isempty(name)
                idxList = [];
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
        
        function str = getDescriptionGroupBy(cd)
            % get description of attributes being grouped by
            desc = cd.description;
            attrDesc = strjoin(cellfun(@(name, vals) sprintf('%s (%d)', name, vals), ... 
                cd.attributeNamesGroupBy, num2cell(cd.nValuesByAttributeGroupBy), 'UniformOutput', false), ' x ');
            if isempty(desc)
                str = attrDesc;
            else
                str = sprintf('%s : %s', desc, attrDesc);
            end
        end

        function printDescription(ci) 
            tcprintf('yellow', '%s: ', class(ci));
            nAttr = ci.nAttributes;
            % print full list
            str = sprintf('{bright blue}Group by ');
            if nAttr > 0
                % print attribute list with value list counts
                for i = 1:ci.nAttributesGroupBy
                    if i > 1
                        str = [str '{dark gray} x '];
                    end
                    str = [str sprintf('{white}%s {bright blue}(%d) ', ...
                        ci.groupByList{i}, ci.nValuesByAttributeGroupBy(i))];
                end
                
            else
                str = '{bright blue}no group by attributes';
            end
            str = [str '\n'];
            tcprintf('inline', str);
            
%             fprintf('All attributes: ');
             if nAttr > 0
%                 % print attribute list with value list counts
%                 for i = 1:ci.nAttributes
%                     % bright color if we're grouping on this
%                     if ci.isAttributeInGroupByList(i) 
%                         nameColorStr = '{bright blue}';
%                     else
%                         nameColorStr = '{blue}';
%                     end
%                     % is the requestAs the same as the attribute name?
%                     if strcmp(ci.attributeNames{i}, ci.attributeRequestAs{i})
%                         nameStr = sprintf('%s%s', nameColorStr, ci.attributeNames{i});
%                     else
%                         nameStr = sprintf('%s%s as %s', nameColorStr, ci.attributeNames{i}, ci.attributeRequestAs{i});
%                     end
%                     tcprintf('inline', strcat(nameStr, '{gray} ({white}%d{gray}) '), ci.nValuesByAttribute(i));
% 
%                     if i < ci.nAttributes
%                         tcprintf('dark gray', ' x ');
%                     end
%                 end
%                 
%                 fprintf('\n');

                % print attribute value lists on each line
                for i = 1:ci.nAttributes
                    if ci.isAttributeInGroupByList(i) 
                        tcprintf('inline', '\t{gray}%s: {white}%s\n', ...
                            ci.attributeNames{i}, strjoin(ci.attributeValueList{i}, ', '));
                    else
                        tcprintf('inline', '\t{gray}%s: {none}%s\n', ...
                            ci.attributeNames{i}, strjoin(ci.attributeValueList{i}, ', '));
                    end
                end
            else
                %tcprintf('dark gray', 'no attributes\n');
            end

        end

        function disp(ci)
            ci.printDescription();
            fprintf('\n');
            builtin('disp', ci);
        end
    end
    
    methods % Adding attributes, setting value lists
        function ci = addAttribute(ci, name, varargin)
            ci.warnIfNoArgOut(nargout);
            
            p = inputParser;
            p.addRequired('name', @ischar);
            % is this attribute always numeric?
            % list of allowed values for this value (other values will be ignored)
            p.addParamValue('requestAs', '', @ischar);
            p.addParamValue('valueList', {}, @(x) isnumeric(x) || iscell(x)); 
            % list of names to substitute for each value in the list
            p.addParamValue('groupBy', true, @islogical);
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
                ci.attributeValueListManual{iAttr} = {};
                ci.attributeValueListSpecified(iAttr) = false;
            else
                if ~iscell(valueList)
                    valueList = num2cell(valueList);
                end
                ci.attributeValueListManual{iAttr} = valueList;
                ci.attributeValueListSpecified(iAttr) = true;
            end
            
            if p.Results.groupBy
                ci.groupByList{end+1} = name;
            end
            
            ci = ci.invalidateCache();
        end

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
        
        function ci = setValueList(ci, name, valueList)
            ci.warnIfNoArgOut(nargout);
            
            iAttr = ci.getAttributeIdx(name);
            
            ci.attributeNumeric(iAttr) = isnumeric(valueList) || islogical(valueList);
            if isempty(valueList)
                ci.attributeValueListManual{iAttr} = {};
                ci.attributeValueListSpecified(iAttr) = false;
            else
                if ~iscell(valueList)
                    valueList = num2cell(valueList);
                end
                ci.attributeValueListManual{iAttr} = valueList;
                ci.attributeValueListSpecified(iAttr) = true;
            end
                       
            ci = ci.invalidateCache();
        end
        
        function ci = groupBy(ci, varargin)
            if iscell(varargin{1})
                attributes = varargin{1};
            else
                attributes = varargin;
            end
            
            if ~isnumeric(attributes)
                % check all exist
                ci.getAttributeIdx(attributes);
            else
                attributes = ci.attributeNames(attributes);
            end
            
            ci.warnIfNoArgOut(nargout);
            ci.groupByList = attributes;
            ci = ci.invalidateCache();
        end

        function ci = groupByAll(ci)
            ci.warnIfNoArgOut(nargout);
            ci = ci.groupBy(ci.attributeNames);
        end
        
        function ci = dontGroupBy(ci, varargin)
            if iscell(varargin{1})
                attributes = varargin{1};
            else
                attributes = varargin;
            end
            
            if ~isnumeric(attributes)
                % check all exist
                ci.getAttributeIdx(attributes);
            else
                attributes = ci.attributeNames(attributes);
            end
            
            ci.warnIfNoArgOut(nargout);
            ci.groupByList = setdiff(ci.groupByList, attributes);
            ci = ci.invalidateCache();
        end
        
        function ci = ungroup(ci)
            ci.groupByList = {};
            ci = ci.invalidateCache();
        end
    end

    methods % Filtering
        function [ci, mask] = filteredByAttribute(ci, attributeName, valueListKeep, varargin)
            ciOrig = ci;
            ci.warnIfNoArgOut(nargout);
            ci = ci.copyIfHandle(); 

            % filter attribute by the intersection of its current value list and valueListKeep 
            p = inputParser;
            p.addRequired('attributeName', @ischar);
            p.addRequired('valueListKeep', @(x) true);
            p.addParamValue('removeFromGroupBy', false, @islogical); 
            p.parse(attributeName, valueListKeep, varargin{:});

            idx = ci.getAttributeIdx(attributeName);
            valueList = ci.attributeValueList{idx};

            if ~any(ismemberCell(valueListKeep, valueList))
                error('No conditions will be kept by this filter');
            end

            if ischar(valueListKeep)
                valueListKeep = {valueListKeep};
            end
            
            if ~isempty(ci.attributeValueList{idx})
                 % maintain the original sort order when filtering
                [ci.attributeValueList{idx} indKeep] = intersectCell(valueList, valueListKeep, 'stable');
            else
                ci.attributeValueList{idx} = valueListKeep;
                indKeep = [];
            end
            
            if p.Results.removeFromGroupBy
                ci.groupByList = setdiff(ci.groupByList, attributeName);
            end
            ci = ci.invalidateCache();

            % generate mask of conditionsKept
            if nargout > 1
                mask = TensorUtils.maskSelectAlongDimension(ciOrig.conditionsSize, idx, indKeep);
            end
        end

        function [ci, mask] = filteredByAttributeStruct(ci, attributeValues, varargin)
            ci.warnIfNoArgOut(nargout);
            ci = ci.copyIfHandle(); 
            
            attributes = fieldnames(attributeValues);
            
            for i = 1:length(attributes)
                [ci mask] = ci.filteredByAttribute(attributes{i}, attributeValues.(attributes{i}), varargin{:});
            end
        end
        
        function ci = withoutAttribute(ci, name)
            % remove an attribute from the list
            ci.warnIfNoArgOut(nargout);
            ci = ci.copyIfHandle();
            idx = ci.getAttributeIdx(name);
            mask = true(ci.nAttributes, 1);
            mask(idx) = false;
            
            ci = ci.maskAttributes(mask);
            
            % also remove from groupByList
            ci.groupByList = setdiff(ci.groupByList, name);
            ci = ci.invalidateCache();
        end
    end
    
    methods(Access=protected)
       function ci = maskAttributes(ci, mask)
           ci.warnIfNoArgOut(nargout);
            ci.attributeNames = ci.attributeNames(mask);
            ci.attributeRequestAs = ci.attributeRequestAs(mask);
            ci.attributeNumeric = ci.attributeNumeric(mask);
            ci.attributeValueList = ci.attributeValueList(mask);
            ci.attributeValueListSpecified = ci.attributeValueListSpecified(mask);
            
            ci.groupByList = setdiff(ci.groupByList, ci.attributeNames);
        end 
    end

    methods % Comparison axis building
        function varargout = compareAlong(ci, attrNames, varargin)
            % identical to compareSlices(attrNames) except attrNames must be a single attr
            assert(ischar(attrNames) || length(attrNames) == 1, 'compareAlong accepts only single attribute comparisons. Use compareSlices');
            [varargout{1:nargout}] = ci.compareSlices(attrNames, varargin{:});
        end 

        function varargout = compareWithin(ci, attrNames, varargin)
            assert(ischar(attrNames) || length(attrNames) == 1, 'compareWithin accepts only single attribute comparisons. Use compareSlicesWithin');
            % identical to compareSlicesWithin(attrNames) except attrNames must be a single attr
            [varargout{1:nargout}] = ci.compareSlicesWithin(attrNames, varargin{:});
        end

        function varargout = compareSlicesWithin(ci, attrNames, varargin)
            % shortcut for .compareSlices( otherAttrNames )
            % build slices across all other attribute so that each inner comparison
            % has conditions which share the same value for attrNames

            otherAttrNames = setdiff(ci.groupByList, attrNames);
            [varargout{1:nargout}] = ci.compareSlices(otherAttrNames, varargin{:});
        end

        function [conditionIdxCell conditionDescriptorOuter conditionDescriptorInnerCell ...
                conditionDescriptorInnerCommon] = ...
                compareSlices(ci, attrNames, varargin)
            % compareSlice is used to generate comparisons where we consider collectively
            % conditions with each set of values of some subset of attributes (attrNamesOrInds)
            % repeating this analysis for each set of values of all the other attributes.
            % In other words, generate a set of slices over attrNamesOrInds, for each set of 
            % other attribute values.
            %
            % Suppose we have four attributes (A, B, C, D) with value counts nA, nB, nC, nD
            % Calling compareSlices({'A', 'B'}) will generate a cell tensor of size 
            % nC x nD. Within each element is a nA x nB tensor of condition indices. 
            % The indices in conditionIdxCell{iC, iD}{iA, iB} represent conditions having value 
            % iA, iB for attributes A, B and values iC, iD for attribute C, D. 
            % The purpose of this reorganization is that it makes it easy to run a comparison
            % involving every condition along the A, B axes, while holding attributes C, D
            % constant, and then repeating this for each value of C, D.
            %
            % Define szOuter as the sizes of the dimensions not in attrNamesOrInds.
            % Define szInner as the sizes of the dimensions in attrNamesOrInds.
            % In the example: szOuter == [nC nD], szInner = [nA nB]
            % 
            % conditionIdxCell : szOuter cell tensor of szInner numeric tensors.
            %   conditionIdxCell{iC, iD}{iA, iB} has the conditionIdx for iA,iB,iC,iD
            %
            % conditionDescriptorOuter : scalar ConditionDescriptor instance, formed by grouping
            %   on attributes not in attrNamesOrInds. This describes the layout of conditions selected
            %   in the outer tensor over C, D. Each inner tensor of conditionIdx will have the corresponding
            %   iC, iD values for C, D in all conditions within.
            %
            % conditionDescriptorInnerCell : szOuter cell tensor of Condition Descriptor instances.
            %   Each instance is similar to conditionDescriptor, but it also filters for the single
            %   attribute values for C, D, and thus perfectly describes the conditions within the corresponding 
            %   conditionIdxCell inner tensor
            %
            % conditionDescriptorInnerCommon: scalar ConditionDescriptor instance, formed
            %   grouping on attrNamesOrInds only. This condition descriptor is common to 
            %   the structure of each inner conditionIdx tensor's comparisons, i.e. it
            %   describes the layout of conditions over A, B. 
            %
            
            p = inputParser;
            p.parse(varargin{:});
        
            % ensure attributes are in groupByList
            attrIdx = makecol(ci.getAttributeIdxInGroupByList(attrNames));

            [otherAttrNames otherAttrIdx] = setdiff(ci.groupByList, attrNames);

            % generate the regrouped conditionInd tensor 
            tInds = ci.conditionsAsLinearInds;
            conditionIdxCell = TensorUtils.regroupAlongDimension(tInds, otherAttrIdx);

            sz = ci.conditionsSize;
            szOuter = TensorUtils.expandScalarSize(sz(otherAttrIdx));
            szInner = TensorUtils.expandScalarSize(sz(attrIdx));

            if nargout > 1
                conditionDescriptorOuter = ci.copyIfHandle().groupBy(otherAttrNames);
            end

            if nargout > 2
                conditionDescriptorInnerCell = cell(szOuter);
                for iOuter = 1:prod(szOuter)
                    conditionDescriptorInnerCell{iOuter} = ci.filteredByAttributeStruct(...
                        conditionDescriptorOuter.conditionsWithGroupByFieldsOnly(iOuter), ...
                        'removeFromGroupBy', true);
                end
            end
            
            if nargout > 3
                conditionDescriptorInnerCommon = ci.copyIfHandle();
                conditionDescriptorInnerCommon = conditionDescriptorInnerCommon.groupBy(attrNames);
            end
            
        end
    end

    
    methods(Static) % Default nameFn and appearanceFn
        function name = defaultNameFn(ci, attrValues, conditionSubs)
            name = '';
            attr = fieldnames(attrValues);
            attrIsNumeric = ci.getIsAttributeNumeric(attr);

            for iAttr = 1:length(attr)
                include = false;
                val = attrValues.(attr{iAttr});
                
                if ~ischar(val) && isscalar(val)
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
                    name = [name attr{iAttr} '=' val ' '];
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

            conditionsSize = ci.conditionsSize;
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
            if nargin < 2
                cdNew = ConditionDescriptor();
            end
            cdNew.noUpdateCache = true;

            cdNew.attributeNames = cd.attributeNames;
            cdNew.attributeNumeric = cd.attributeNumeric; 
            cdNew.attributeValueList = cd.attributeValueList; 

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
            
            cdNew.noUpdateCache = false;
            cdNew = cdNew.invalidateCache();
        end
    end
    
    methods 
        function obj = copyIfHandle(obj)
            if isa(obj, 'handle')
                obj = obj.copy();
            end
        end
    end

    methods(Access=protected) % Utility methods
        function warnIfNoArgOut(obj, nargOut)
            if nargOut == 0 && ~isa(obj, 'handle')
                message = sprintf('WARNING: %s is not a handle class. If the instance handle returned by this method is not stored, this call has no effect.\\n', ...
                    class(obj));
                expr = sprintf('debug(''%s'')', message);
                evalin('caller', expr); 
            end
        end
    end
end

