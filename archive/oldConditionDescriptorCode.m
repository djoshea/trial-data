    
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

        function [conditionIdxCell, conditionDescriptorOuter, conditionDescriptorInnerCell, ...
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
            %attrIdx = makecol(ci.getAttributeIdxInGroupByList(attrNames));

            [otherAttrNames, otherAttrIdx] = setdiff(ci.groupByList, attrNames);

            % generate the regrouped conditionInd tensor 
            tInds = ci.conditionsAsLinearInds;
            conditionIdxCell = TensorUtils.regroupAlongDimension(tInds, otherAttrIdx);

            sz = ci.conditionsSize;
            szOuter = TensorUtils.expandScalarSize(sz(otherAttrIdx));
            %szInner = TensorUtils.expandScalarSize(sz(attrIdx));

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