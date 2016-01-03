classdef TensorUtils
    % set of classes for building, manipulating, and computing on
    % high-d (or arbitrary-d) matrices easily
    methods(Static) % Mapping, construction via callback
        function varargout = mapToSizeFromSubs(sz, varargin)
            % t = mapTensor(sz, contentsFn = @(varargin) NaN, asCell = false)
            % build a tensor with size sz by passing subscripts inds to
            % contentsFn(sub1, sub2, ...) maps subscript indices as a vector to the contents
            % asCell == true --> returns cell, asCell == false returns matrix, defaults to false
            
            p = inputParser;
            p.addRequired('size', @(x) isempty(x) || (isvector(x) && isnumeric(x)));
            p.addOptional('contentsFn', [], @(x) isa(x, 'function_handle'));
            p.addOptional('asCell', false, @islogical);
            p.parse(sz, varargin{:});
            asCell = p.Results.asCell;
            contentsFn = p.Results.contentsFn;
            
            if isempty(sz)
                varargout = cellvec(nargout);
                for i = 1:nargout
                    if asCell
                        varargout{i} = {};
                    else
                        varargout{i} = [];
                    end
                end
                return
            end
            
            sz = TensorUtils.expandScalarSize(sz);
            nDims = length(sz);
            idxEachDim = arrayfun(@(n) 1:n, sz, 'UniformOutput', false);
            [subsGrids{1:nDims}] = ndgrid(idxEachDim{:});
            
            if isempty(contentsFn)
                if asCell
                    contentsFn = @(varargin) {};
                else
                    contentsFn = @(varargin) NaN;
                end
            end
            
            [varargout{1:nargout}] = arrayfun(contentsFn, subsGrids{:}, 'UniformOutput', ~asCell);
        end
        
        function varargout = map(fn, varargin)
            % works just like cellfun or arrayfun except auto converts each arg
            % to a cell so that cellfun may be used. Returns a cell arrays with
            % the same size as the tensor, although cell arrays containing
            % the same 
            for iArg = 1:length(varargin)
                if ~iscell(varargin{iArg})
                    varargin{iArg} = num2cell(varargin{iArg});
                end
            end
            [varargout{1:nargout}] = cellfun(fn, varargin{:}, 'UniformOutput', false);
            % convert scalar numeric cells back to matrices
%             for iArg = 1:numel(varargout);
%                 if all(cellfun(@(x) isnumeric(x) && isscalar(x), varargout{iArg}));
%                     varargout{iArg} = cell2mat(varargout{iArg});
%                 end
%             end
        end
        
        function results = mapIncludeSubs(fn, varargin)
            % mapWithInds(fn, t1, t2, ...) calls fn(t1(subs), t2(subs), ..., subs) with subs
            % being the subscript indices where the element of t1, t2, etc.
            % was extracted
            
            for iArg = 1:length(varargin)
                if ~iscell(varargin{iArg})
                    varargin{iArg} = num2cell(varargin{iArg});
                end
            end
            tSubs = TensorUtils.containingSubscripts(size(varargin{1}));
            results = cellfun(fn, varargin{:}, tSubs, 'UniformOutput', false);
        end
        
        function varargout = mapIncludeSubsAndSize(fn, varargin)
            % mapWithInds(fn, t1, t2, ...) calls fn(t1(subs), t2(subs), ..., subs, sz) with subs
            % being the subscript indices where the element of t1, t2, etc.
            % was extracted and sz being size(t1) == size(t2).
            
            sz = size(varargin{1});
            fnWrap = @(varargin) fn(varargin{:}, sz);
            [varargout{1:nargout}] = TensorUtils.mapIncludeSubs(fnWrap, varargin{:});
        end
        
        function varargout = mapSlices(fn, spanDim, varargin)
            % varargout = mapSlices(fn, spanDims, varargin)
            %
            % this acts like map, calling fn(varargin{1}(ind),varargin{2}(ind))
            % except rather than being called on each element of varargin{:}
            % individually, it is called on slices of the tensor(s) at once. These slices
            % are created by selecting all elements along the dimensions in spanDims and
            % repeating this over each set of subscripts along the other dims.
            % The slices passed to fn will not be squeezed, so they will
            % have singleton dimensions for dim in spanDim. Call
            % .squeezeDims(in, spanDim) to obtain a squeezed slice.
            %
            % The result will be reassembled into a tensor, whose size is determined by
            % the sizes of dimensions not in spanDim. Because the output
            % values will be stored as cell tensor elements, there are no
            % constraints on what these outputs look like
            
            %sz = size(varargin{1});
            nd = ndims(varargin{1});
            %nArgs = length(varargin);
            
            % we select individual slices by selecting each along the non-spanned dims
            dim = setdiff(1:nd, spanDim);
            
            % slice through each of the varargin
            tCellArgs = cellfun(@(t) TensorUtils.selectEachAlongDimension(t, dim), ...
                varargin, 'UniformOutput', false);
            
            % run the function on each slice
            [resultCell{1:nargout}] = cellfun(fn, tCellArgs{:}, 'UniformOutput', false);
            
            varargout = resultCell;
            
            % (old) reassemble the result
            % varargout = cellfun(@(r) TensorUtils.reassemble(r, dim), resultCell, 'UniformOutput', false);
        end
        
        function varargout = mapSlicesInPlace(fn, spanDim, varargin)
            % This function acts very similarly to mapSlices. The only
            % difference is that fn must return outputs that has the same
            % size and shape as it's inputs. Provided this constraint is
            % met, the output will be a tensor that has the same shape as
            % the input tensor. Effectively you loop through the
            % input tensors one slice at a time, transform that slice in
            % place via slice = fn(slice) and rebuild the output tensor one
            % slice at a time.
            
            [resultCell{1:nargout}] = TensorUtils.mapSlices(fn, spanDim, varargin{:});
            nonSpanDims = TensorUtils.otherDims(size(resultCell{1}), spanDim, ndims(varargin{1}));
            varargout = cellfun(@(r) TensorUtils.reassemble(r, nonSpanDims, ndims(varargin{1})), ...
                resultCell, 'UniformOutput', false);
        end
        
        function varargout = mapFromAxisLists(fn, axisLists, varargin)
            % varargout = mapToSizeFromAxisElements(sz, [fn = @(varargin) varargin], axis1, axis2, ...)
            % build a tensor with size sz by setting
            % tensor(i,j,k,...) = fn(axis1{i}, axis2{j}, axis3{k})
            % (or axis(i) for numeric arrays)
            % if fn is omitted returns a cell array of
            % {axis1{i}, axis2{j}, axis3{j}}
            
            p = inputParser;
            p.addParamValue('asCell', true, @islogical);
            p.parse(varargin{:});
            
            if isempty(fn)
                fn = @(varargin) varargin;
            end
            if iscell(axisLists)
                sz = cellfun(@numel, axisLists);
            else
                sz = arrayfun(@numel, axisLists);
            end
            
            [varargout{1:nargout}] = TensorUtils.mapToSizeFromSubs(sz, @indexFn, 'asCell', p.Results.asCell);
            
            function varargout = indexFn(varargin)
                inputs = cellvec(numel(axisLists));
                for iAx = 1:numel(axisLists)
                    if iscell(axisLists{iAx})
                        inputs{iAx} = axisLists{iAx}{varargin{iAx}};
                    else
                        inputs{iAx} = axisLists{iAx}(varargin{iAx});
                    end
                end
                [varargout{1:nargout}] = fn(inputs{:});
            end
            
        end
        
        function t = buildCombinatorialStructTensor(varargin)
            % Given all struct vector arguments passed in, constructs a 
            % tensor struct array containing the merged version of each
            % struct array
            
            nByArg = cellfun(@numel, varargin);
            fieldsByArg = cellfun(@fieldnames, varargin, 'UniformOutput', false);
            allFields = cat(1, fieldsByArg{:});
                 
            asCells = cell(nargin, 1);
            nArg = numel(varargin);
            for iArg = 1:nArg
                % get as vector with size nByArg(iArg) entries
                flatVector = varargin{iArg}(:);
                % struct2cell returns nFields x nByArg(iArg) cell array,
                % which we transpose
                thisAsCell = struct2cell(flatVector)';
                % which we orient to have size nFields along dim nargin+1
                % and size nByArg(iArg) along dim iArg
                thisAsCellOriented = TensorUtils.orientSliceAlongDims(thisAsCell, [iArg, nArg+1]);
                % and expand to be the same size as the full combinatorial tensor
                asCells{iArg} = TensorUtils.singletonExpandToSize(thisAsCellOriented, nByArg);
            end
                
            % combined will be size [nByArg(:) nFieldsTotal)
            combined = cat(nargin+1, asCells{:});
            t = cell2struct(combined, allFields, nArg+1);
        end
        
        function t = mapCatToTensor(fn, varargin)
            % like arrayfun, except with UniformOutput false, and the
            % results will be concatenated together using cellmat to form a
            % tensor
            t = TensorUtils.map(fn, varargin{:});
            if iscell(t)
                t = cell2mat(t);
            end
        end
    end
    
    methods(Static) % Dimensions and sizes
        function out = emptyWithSameType(in, szOut)
            % create out as size szOut with same type as in (either cell or
            % nans)
            if iscell(in)
                out = cell(szOut);
            else
                out = nan(szOut);
            end
        end
        
        function tf = isvector(vec)
            % like isvector except works for high-d vectors as well
            sz = size(vec);
            tf = nnz(sz ~= 1) == 1 && ~any(sz == 0);
        end
        
        function sz = sizeMultiDim(t, dims)
            % sz = sizeMultiDim(t, dims) : sz(i) = size(t, dims(i))
            szAll = TensorUtils.expandSizeToNDims(size(t), max(dims));
            sz = arrayfun(@(d) szAll(d), dims);
        end
        
        function sz = sizeOtherDims(t, excludeDims)
            % sz = sizeOtherDims(t, excludeDims)
            otherDims = TensorUtils.otherDims(size(t), excludeDims);
            sz = TensorUtils.sizeMultiDim(t, otherDims);
        end
        
        function sz = expandScalarSize(sz)
            % if sz (size) is a scalar, make it into a valid size vector by
            % appending 1 to the end. i.e. 3 --> [3 1]
            if isempty(sz)
                sz = [0 0];
            elseif isscalar(sz)
                sz = [sz 1];
            end
        end
        
        function sz = sizeNDims(t, nDims)
            sz = TensorUtils.expandSizeToNDims(size(t), nDims);
        end
        
        % pads sz with 1s to make it nDims length
        function sz = expandSizeToNDims(sz, nDims)
            szPad = ones(nDims - numel(sz), 1);
            sz(end+1:nDims) = szPad;
        end
        
        function tf = compareSizeVectors(sz1, sz2)
            nDims = max(numel(sz1), numel(sz2));
            tf = isequal(TensorUtils.expandSizeToNDims(sz1, nDims), ....
                TensorUtils.expandSizeToNDims(sz2, nDims));
        end
        
        function other = otherDims(sz, dims, ndims)
            % otherDims(t, dims, [ndims]) returns a list of dims in t NOT in dims
            % e.g. if ndims(t) == 3, dims = 2, other = [1 3]
            if nargin < 3
                ndims = numel(sz);
            end
            allDims = 1:ndims;
            other = makerow(setdiff(allDims, dims));
        end
        
        function d = firstNonSingletonDim(t)
            sz = size(t);
            d = find(sz > 1, 1, 'first');
            if isempty(d)
                d = 1;
            end
        end
    end
    
    methods(Static) % Mask generation and mask utilities
       function idx = vectorMaskToIndices(mask)
            if islogical(mask);
                idx = makecol(find(mask));
            else
                idx = makecol(mask);
            end
        end
        
        function mask = vectorIndicesToMask(idx, len)
            if islogical(idx)
                mask = makecol(idx);
            else
                mask = falsevec(len);
                mask(idx) = true;
            end
        end
        
        function maskNew = subselectVectorMask(maskOrig, maskAnd)
            % given logical maskOrig
            existingKeep = find(maskOrig);
            assert(numel(maskAnd) == numel(existingKeep), 'Subselection mask does not match nnz of original mask');
            keep = TensorUtils.vectorIndicesToMask(existingKeep(maskAnd), numel(maskOrig));
            maskNew = maskOrig;
            maskNew(~keep) = false;
        end
        
        function inflated = inflateMaskedTensor(maskedTensor, dims, masks, fillWith)
            % inflated = inflateMaskedTensor(maskedTensor, dims, masks, fillWith=NaN)
            % takes maskedTensor, which has been formed by slicing some
            % original tensor by selecting with mask along each dimension
            % in dims, and returns the original tensor where the masked out
            % rows, cols, etc. are filled with fill (defaults to NaN)
            % if dims is a scalar, masks is a logical or numeric vector to use
            % for selecting along dim. If dim is a vector, select is a cell array of
            % vectors to be used for selecting along dim(i)
            
            if ~iscell(masks)
                masks = {masks};
            end
            assert(all(cellfun(@(x) islogical(x) && isvector(x), masks)), ...
                'Mask must be (cell of) logical vectors. Use vectorIndicesToMask to convert.');
           
            if numel(masks) == 1
                masks = repmat(masks, numel(dims), 1);
            end
            assert(numel(masks) == numel(dims), 'Number of dimensions must match number of masks provided');
            masks = makecol(masks);
                
            nInflatedVec = cellfun(@numel, masks);
            nSelectedVec = cellfun(@nnz, masks);
            
            % check size along selected dims
            assert(isequal(makecol(TensorUtils.sizeMultiDim(maskedTensor, dims)), nSelectedVec), ...
                'Size along each masked dimension must match nnz(mask)');
            
            % compute inflated size
            inflatedSize = size(maskedTensor);
            inflatedSize(dims) = nInflatedVec;
            
            if nargin < 4
                if iscell(maskedTensor)
                    fillWith = {[]};
                else
                    fillWith = NaN;
                end
            end
            inflated = repmat(fillWith, inflatedSize);
            
            maskByDim = TensorUtils.maskByDimCellSelectAlongDimension(inflatedSize, dims, masks);
            inflated(maskByDim{:}) = maskedTensor;
        end
        
        function idxFull = indicesIntoMaskToOriginalIndices(idxIntoMasked, mask)
            maskInds = find(mask);
            idxFull = maskInds(idxIntoMasked);
        end

    end
    
    methods(Static) % Indices and subscripts
      
        function t = containingLinearInds(sz)
            % build a tensor with size sz where each element contains the linear
            % index it would be accessed at, e.g. t(i) = i
            sz = TensorUtils.expandScalarSize(sz);
            t = TensorUtils.mapToSizeFromSubs(sz, @(varargin) sub2ind(sz, varargin{:}), false);
        end
        
        function t = containingSubscripts(sz, asCell)
            sz = TensorUtils.expandScalarSize(sz);
            
            % asCell == true means each element is itself a cell rather then a vector of
            % subscripts
            if nargin < 2
                asCell = false;
            end
            
            % build a tensor with size sz where each element contains the vector
            % of subscripts it would be accessed at, e.g. t(i) = i
            if asCell
                t = TensorUtils.mapToSizeFromSubs(sz, @(varargin) varargin, true);
            else
                t = TensorUtils.mapToSizeFromSubs(sz, @(varargin) [varargin{:}]', true);
            end
        end

        function subsCell = ndgridCell(sz)
            % works like ndgrid except captures all outputs
            
            args = arrayfun(@(s) 1:s, sz, 'UniformOutput', false);
            [subsCell{1:ndims(sz)}] = ndgrid(args{:});
        end
        
        function t = containingSubscriptsCatAlongFirstDimension(sz)
            % t is a tensor whose size is [numel(sz) sz]
            % t(dim, i, j, k) is the subscript of i,j,k on dimension dim.
            % e.g. if dim == 2, t(2, i, j, k, ...) == j
            t = cell2mat(shiftdim(TensorUtils.containingSubscripts(sz), -1));
        end
        
        function mat = ind2subAsMat(sz, inds)
            sz = TensorUtils.expandScalarSize(sz);
            
            % sz is the size of the tensor
            % mat is length(inds) x length(sz) where each row contains ind2sub(sz, inds(i))
            
            ndims = length(sz);
            subsCell = cell(ndims, 1);
            
            [subsCell{:}] = ind2sub(sz, makecol(inds));
            
            mat = [subsCell{:}];
        end
        
        function inds = subMat2Ind(sz, mat)
            % inds = subMat2Ind(sz, mat)
            % sz is the size of the tensor
            % mat is length(inds) x length(sz) where each row contains a set of subscripts
            % as would be returned by ind2sub(sz, inds(i))
            % subMat2Ind essentially converts back to linear indices using
            % sub2ind and is the inverse of ind2subAsMat
            
            sz = TensorUtils.expandScalarSize(sz);
            
            ndims = length(sz);
            % DO NOT UNCOMMENT. THIS WILL BREAK THINGS SINCE SZ CAN HAVE
            % SZ(1) == 1.
             if ndims == 2 && sz(2) == 1 && size(mat, 2) == 1
                 ndims = 1;
             end
            subsCell = arrayfun(@(dim) mat(:, dim), 1:ndims, 'UniformOutput', false);
            
            inds = sub2ind(sz, subsCell{:});
        end
    end
    
    methods(Static) % Selection Mask generation
        function maskByDim = maskByDimCell(sz)
            sz = TensorUtils.expandScalarSize(sz);
            
            % get a cell array of selectors into each dim that would select
            % every element if used via t(maskByDim{:})
            maskByDim = arrayfun(@(n) true(n, 1), sz, 'UniformOutput', false);
        end
        
        % the next few methods accept a dim and select argument
        % if dim is a scalar, select is a logical or numeric vector to use
        % for selecting along dim. If dim is a vector, select is a cell array of
        % vectors to be used for selecting along dim(i)
        function maskByDim = maskByDimCellSelectAlongDimension(sz, dim, select)
            sz = TensorUtils.expandScalarSize(sz);
            
            % get a cell array of selectors into each dim that effectively select
            % select{i} along dim(i). These could be used by indexing a tensor t
            % via t(maskByDim{:}) --> se selectAlongDimension
            if ~iscell(select)
                select = {select};
            end
            
            assert(length(dim) == length(select), 'Number of dimensions must match length of select mask cell array');
            maskByDim = TensorUtils.maskByDimCell(sz);
            maskByDim(dim) = select;
        end
        
        function mask = maskSelectAlongDimension(sz, dim, select)
            sz = TensorUtils.expandScalarSize(sz);
            
            % return a logical mask where for tensor with size sz
            % we select t(:, :, select, :, :) where select acts along dimension dim
            
            mask = false(sz);
            maskByDim = TensorUtils.maskByDimCellSelectAlongDimension(sz, dim, select);
            mask(maskByDim{:}) = true;
        end
        
        function maskByDim = maskByDimCellSelectPartialFromOrigin(sz, dim, select)
            maskByDim = TensorUtils.maskByDimCell(sz);
            maskByDim(dim) = arrayfun(@(n) 1:n, select, 'UniformOutput', false);
        end
        
        function mask = maskSelectPartialFromOrigin(sz, dim, select)
            % select the first select(i) items from dim(i)
            % for i = 1:numel(sz)
            selectCell = arrayfun(@(n) 1:n, select, 'UniformOutput', false);
            mask = TensorUtils.maskSelectAlongDimension(sz, dim, selectCell);
        end
    end
    
    methods(Static) % Multi-dim extensions of any, all, find etc.
        function t = anyMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = any(t, dims(iD));
            end
        end
        
        function t = allMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = all(t, dims(iD));
            end
        end

        function idxTensor = findNAlongDim(t, dim, N, direction)
            % idxTensor = findNAlongDim(t, dim, N, direction)
            % finds the first/last N non-zero values in t along dimension t
            % proceeds in direction 'first' (default) or 'last'
            % if less than N non-zero values are found, idxTensor will be
            % NaN. idxTensor is the same size as t along all dimensions
            % except dim, where it has size N.
            
            if nargin < 4
                direction = 'first';
            end
            if ~islogical(t)
                t = t ~= 0;
            end
            
            % create a template for findInner to use with the correct size
            % and orientation
            sizeSlice = onesvec(max(dim, ndims(t)))';
            sizeSlice(dim) = N;
            emptySlice = nan(sizeSlice);
            
            idxTensor = cell2mat(TensorUtils.mapSlices(@findInner, dim, t));
            
            function w = findInner(v)
                w = emptySlice;
                idx = find(v, N, direction);
                w(1:numel(idx)) = idx;
            end
        end
    end
    
    methods(Static) % Squeezing along particular dimensions
        function tsq = squeezeDims(t, dims)
            % like squeeze, except only collapses singleton dimensions in list dims
            siz = size(t);
            dims = dims(dims <= ndims(t));
            dims = dims(siz(dims) == 1);
            siz(dims) = []; % Remove singleton dimensions.
            siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
            tsq = reshape(t,siz);
        end
        
        function newDimIdx = shiftDimsPostSqueeze(sz, squeezeDims, dimsToShift)
            assert(isvector(sz), 'First arg must be size');
            % when squeezing along squeezeDims, the positions of dims in
            % dimsToShift will change. The new dim idx will be returned
            origDims = 1:length(sz);
            squeezeDims = squeezeDims(squeezeDims <= length(sz));
            remainDims = setdiff(origDims, squeezeDims);
            
            [~, newDimIdx] = ismember(dimsToShift, remainDims);
            newDimIdx(newDimIdx==0) = NaN;
        end
        
        function tsq = squeezeOtherDims(t, dims)
            other = TensorUtils.otherDims(size(t), dims);
            tsq = TensorUtils.squeezeDims(t, other);
        end
    end
    
    methods(Static) % Regrouping, Nesting, Selecting, Reshaping
        function tCell = regroupAlongDimension(t, dims)
            % tCell = regroupAlongDimension(t, dims)
            % returns a cell tensor of tensors, where the outer tensor is over
            % the dims in dims. Each inner tensor is formed by selecting over
            % the dims not in dims.
            %
            % e.g. if size(t) = [nA nB nC nD] and dims is [1 2],
            % size(tCell) = [nA nB] and size(tCell{iA, iB}) = [nC nD]
            
            tCell = TensorUtils.squeezeSelectEachAlongDimension(t, dims);
            tCell = TensorUtils.squeezeOtherDims(tCell, dims);
        end
        
        function tCell = nestedRegroupAlongDimension(t, dimSets)
            assert(iscell(dimSets), 'dimSets must be a cell array of dimension sets');
            
            dimSets = makecol(cellfun(@makecol, dimSets, 'UniformOutput', false));
            allDims = cell2mat(dimSets);
            
            assert(length(unique(allDims)) == length(allDims), ...
                'A dimension was included in multiple dimension sets');
            
            otherDims = TensorUtils.otherDims(size(t), allDims);
            if ~isempty(otherDims)
                dimSets{end+1} = otherDims;
            end
            
            tCell = inner(t, dimSets);
            return;
            
            function tCell = inner(t, dimSets)
                if length(dimSets) == 1
                    % special case, no grouping, just permute dimensions and
                    % force to be cell
                    tCell = permute(t, dimSets{1});
                    if ~iscell(tCell)
                        tCell = num2cell(tCell);
                    end
                elseif length(dimSets) == 2
                    % last step in recursion, call final regroup
                    tCell = TensorUtils.regroupAlongDimension(t, dimSets{1});
                else
                    % call inner on each slice of dimSets{1}
                    remainingDims = TensorUtils.otherDims(size(t), dimSets{1});
                    tCell = TensorUtils.mapSlices(@(t) mapFn(t, dimSets), ...
                        remainingDims, t);
                    tCell = TensorUtils.squeezeOtherDims(tCell, dimSets{1});
                end
            end
            
            function tCell = mapFn(t, dimSets)
                remainingDimSets = cellfun(...
                    @(dims) TensorUtils.shiftDimsPostSqueeze(size(t), dimSets{1}, dims), ...
                    dimSets(2:end), 'UniformOutput', false);
                tCell = inner(TensorUtils.squeezeDims(t, dimSets{1}), remainingDimSets);
            end
            
        end
        
        function [res, maskByDim] = selectAlongDimension(t, dim, select, squeezeResult)
            if nargin < 4
                squeezeResult = false;
            end
            sz = size(t);
            maskByDim = TensorUtils.maskByDimCellSelectAlongDimension(sz, dim, select);
            res = t(maskByDim{:});
            
            if squeezeResult
                % selectively squeeze along dim
                res = TensorUtils.squeezeDims(res, dim);
            end
        end
        
                
        function sel = selectAlongDimensionWithNaNs(t, dim, select, varargin)
            assert(numel(dim) == 1, 'Must be single dimension');
            nanMask = isnan(select(:));
            selNoNan = TensorUtils.selectAlongDimension(t, dim, select(~nanMask), varargin{:});
            sel = TensorUtils.inflateMaskedTensor(selNoNan, dim, ~nanMask);
        end
        
        function [res, maskByDim] = squeezeSelectAlongDimension(t, dim, select)
            % select ind along dimension dim and squeeze() the result
            % e.g. squeeze(t(:, :, ... ind, ...))
            
            [res, maskByDim] = TensorUtils.selectAlongDimension(t, dim, select, true);
        end
        
        function tCell = selectEachAlongDimension(t, dim, squeezeEach)
            % returns a cell array tCell such that tCell{i} = selectAlongDimension(t, dim, i)
            % optionally calls squeeze on each element
            if nargin < 3
                squeezeEach = false;
            end
            
            sz = size(t);
            
            % generate masks by dimension that are equivalent to ':'
            %maskByDimCell = TensorUtils.maskByDimCell(sz);
            
            dimMask = true(ndims(t), 1);
            dimMask(dim) = false;
            szResult = sz;
            szResult(dimMask) = 1;
            
            % oh so clever
            tCell = TensorUtils.mapToSizeFromSubs(szResult, 'asCell', true, ...
                'contentsFn', @(varargin) TensorUtils.selectAlongDimension(t, dim, ...
                varargin(dim), squeezeEach));
        end
        
        function tCell = squeezeSelectEachAlongDimension(t, dim)
            % returns a cell array tCell such that tCell{i} = squeezeSelectAlongDimension(t, dim, i)
            tCell = TensorUtils.selectEachAlongDimension(t, dim, true);
        end
        
        function t = reassemble(tCell, dim, nd)
            % given a tCell in the form returned by selectEachAlongDimension
            % dim is the dimension embedded within each slice
            % return the original tensor
            
            %             if nargin < 3
            %                 nd = ndims(tCell);
            %             end
            szOuter = size(tCell);
            szOuter = [szOuter ones(1, nd - length(szOuter))];
            szInner = size(tCell{1});
            szInner = [szInner ones(1, nd - length(szInner))];
            
            % dimMask(i) true i
            dimMask = false(nd, 1);
            dimMask(dim) = true;
            
            % compute size of result t
            % use outerDims when its in dim, innerDims when it isn't
            szT = nan(1, ndims(tCell));
            szT(dimMask) = szOuter(dimMask);
            szT(~dimMask) = szInner(~dimMask);
            
            % rebuild t by grabbing the appropriate element from tCell
            %subs = TensorUtils.containingSubscripts(szT);
            t = TensorUtils.mapToSizeFromSubs(szT, @getElementT, true);
            
            function el = getElementT(varargin)
                [innerSubs, outerSubs] = deal(varargin);
                % index with dim into tt, non-dim into tt{i}
                [outerSubs{~dimMask}] = deal(1);
                [innerSubs{dimMask}] = deal(1);
                tEl = tCell{outerSubs{:}};
                if iscell(tEl)
                    el = tEl{innerSubs{:}};
                else
                    el = tEl(innerSubs{:});
                end
            end
        end
        
        function vec = flatten(t)
            vec = makecol(t(:));
        end
        
        function mat = flattenAlongDimension(t, dim)
            % returns a 2d matrix where mat(i, :) is the flattened vector of tensor
            % values from each t(..., i, ...) where i is along dim
            
            nAlong = size(t, dim);
            nWithin = numel(t) / nAlong;
            if iscell(t)
                mat = cell(nAlong, nWithin);
            else
                mat = nan(nAlong, nWithin);
            end
            
            sqMask = TensorUtils.maskByDimCell(size(t));
            for iAlong = 1:nAlong
                sqMask{dim} = iAlong;
                within = t(sqMask{:});
                mat(iAlong, :) = within(:);
            end
        end
        
        function tCell = flattenAlongDimensionAsCell(t, dim)
            % returns a cell array of length size(t, dim)
            % where each element is the flattened vector of tensor
            % values from each t(..., i, ...) where i is along dim
            tCell = TensorUtils.regroupAlongDimension(t, dim);
            for iAlong = 1:length(tCell)
                tCell{iAlong} = tCell{iAlong}(:);
            end
        end
        
        function [out, which] = catWhich(dim, varargin)
            % works like cat, but returns a vector indicating which of the
            % inputs each element of out came from
            out = cat(dim, varargin{:});
            if nargout > 1
                which = cell2mat(makecol(cellfun(@(in, idx) idx*ones(size(in, dim), 1), varargin, ...
                    num2cell(1:numel(varargin)), 'UniformOutput', false)));
            end
        end
        
%         function out = catPad(dim, with, varargin)
%             % works like cat, except pads each element to be the same size
%             % along dimension dim before concatenating
%             szAlongDim = cellfun(@(x) size(x, dim), varargin);
%             if numel(unique(szAlongDim)) == 1
%                 % no padding needed
%                 out = cat(dim, varargin{:});
%             else
%                 newSzAlongDim = max(szAlongDim);
%                 paddedCell = cell(size(varargin));
%                 for iA = 1:numel(varargin)
%                     szVec = size(varargin{iA});
%                     szVec(dim) = newSzAlongDim - szVec(dim);
%                     paddedCell{iA} = cat(dim, varargin{iA}, repmat(with, szVec));
%                 end
%                 out = cat(dim, paddedCell{:});
%             end
%         end
        
        function [out, which] = catWhichIgnoreEmpty(dim, varargin)
            % works like cat, but returns a vector indicating which of the
            % inputs each element of out came from
            
            isEmpty = cellfun(@isempty, varargin);
            out = cat(dim, varargin{~isEmpty});
            if isempty(out)
                which = [];
                return;
            end
            if nargout > 1
                whichMasked = cell2mat(makecol(cellfun(@(in, idx) idx*ones(size(in, dim), 1), varargin(~isEmpty), ...
                    num2cell(1:nnz(~isEmpty)), 'UniformOutput', false)));
                
                % whichMasked indexes into masked varargin, reset these to
                % index into the original varargin
                which = TensorUtils.indicesIntoMaskToOriginalIndices(whichMasked, ~isEmpty);
            end
        end
        
        function [out, which] = catInnerDimOverOuterDim(t, innerDim, outerDim)
            % for each entry in t along dimension outerDim, concatenates
            % these entries along dimension, innerDim
            
            if nargout == 2
                [out, whichCell]  = TensorUtils.mapSlices(@(slice) TensorUtils.catWhich(innerDim, slice{:}), outerDim, t);
                which = whichCell{1};
            else
                out = TensorUtils.mapSlices(@(slice) cat(innerDim, slice{:}), outerDim, t);
            end
        end
        
        function [out, labelsByDimOut] = reshapeByConcatenatingDims(in, whichDims, labelsByDim)
            % reshapes a tensor by concatenating dims.
            % whichDims is a cell vector indicating which dimensions of in to
            % capture along each dimension of out.
            %
            % For example, suppose in has size s1 x s2 x s3 x s4 x s5.
            % And whichDims = { [ 1 2 ], [3 5], 4 }
            % Then out will be 3-dimensional with size:
            %    (s1*s2) x (s3*s5) x (s4)
            % such that in(i1, i2, i3, i4, i5) will end up at
            %   out(o1, o2, o3), where:
            %       o1 = sub2ind([s1 s2], i1, i2);
            %       o2 = sub2ind([s3 s5], i3, i5);
            %       o3 = sub2ind(s4, i4) == s4;
            % This ensures that along out dim 1, all elements of in-dim 1
            % will remain together for dim2=1, then all elements along
            % in-dim 1 for dim2=2, etc. That is, the first dimension will
            % be swept with the later dimensions constant. This is the
            % reverse of how you would implement this using nested for loops,
            % i.e. [1 2] means for j = 1:size(in, 2), for i = 1:size(in, 1)
            % and is ultimately because Matlab uses column major matrices
            %
            % LabelsByDimOut allows the user to keep track of what pieces
            % of in end up at which positions in out. LabelsByDimIn
            % represent a set of labels along each of the "axes" or
            % dimensions of in. LabelsByDimOut represent a set of labels
            % along each of the 'axes" of dimension out. Those labels are
            % taken to match the labels of in and are concatenated in the
            % same fashion.
            %
            % labelsByDim (optional) is a ndims(in) cell vector, each
            % element i contains a vector or matrix with size(in, i) rows and
            % an arbitrary number of colums. If labelsByDim{i} is a vector,
            % it MUST be a column vector.
            %
            % labelsByDimOut will have a similar format: a ndims(out) cell
            % where each element i contains a matrix with size(out, i)
            % rows. The contents are formed by concatenating the columns
            % taken from each of the labelsByDim{...} cells.
            %
            % The format of labelsByDimIn and labelsByDimOut is that the
            % output can be used as a input to this function again.
            
            if ~iscell(whichDims)
                whichDims = {whichDims};
            end
            whichDims = makecol(whichDims);
            
            allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
            allDims = cat(1, allDims{:});
            
            szIn = TensorUtils.sizeNDims(in, max(allDims));
            ndimsIn = max(max(allDims), ndims(in));
            
            assert(all(ismember(1:numel(allDims), allDims)), ...
                'whichDims must contain each dim in 1:length(whichDims)');
            % add any trailing dimensions which are missing from the list
            % automatically
            if max(allDims) < ndimsIn;
                whichDims = cat(1, whichDims, num2cell(max(allDims)+1:ndimsIn));
                % recompute allDims in case trailing dims were added
                allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
                allDims = cat(1, allDims{:});
            end

            ndimsOut = length(whichDims);
           
            szOut = nan(1, ndimsOut);
            for iDimOut = 1:ndimsOut
                szOut(iDimOut) = prod(szIn(whichDims{iDimOut}));
            end
            
            % order dimensions of out according to whichDims
            out = reshape(permute(in, allDims), szOut);
            
            % build labels for output dimensions
            if nargout > 1
                szInExpand = TensorUtils.expandScalarSize(szIn);
                labelsByDimTemplate = arrayfun(@(dim) 1:szInExpand(dim), 1:numel(szInExpand), 'UniformOutput', false);
                if ~exist('labelsByDim', 'var')
                    labelsByDim = labelsByDimTemplate;
                else
                    % extend labelsByDim to the full length ndims if necessary
                    labelsByDimTemplate(1:numel(labelsByDim)) = labelsByDim;
                    labelsByDim = labelsByDimTemplate;
                end
                
                szOut = size(out);
                labelsByDimOut = cellvec(ndimsOut);
                for iDimOut = 1:ndimsOut
                    dimsFromIn = whichDims{iDimOut};
                    % temporary storage of the columns of
                    % labelsByDimOut{iDimOut}, before concatenation
                    labelsThisOut = cellvec(numel(dimsFromIn));
                    subsCell = cellvec(numel(dimsFromIn));
                    
                    % get subscripts for each element in the dimsFromIn slice
                    [subsCell{1:numel(dimsFromIn)}] = ind2sub(szIn(dimsFromIn), 1:szOut(iDimOut));
                    for iDimFromIn = 1:numel(dimsFromIn)
                        % check whether it's a coorectly sized row vector and convert to
                        % column vector
                        labelsThisIn = labelsByDim{dimsFromIn(iDimFromIn)};
                        if isvector(labelsThisIn) && length(labelsThisIn) == size(in, dimsFromIn(iDimFromIn))
                            labelsThisIn = makecol(labelsThisIn);
                        end
                        % and pull the correct labels by selecting the
                        % correct rows (typically multiple times)
                        labelsThisOut{iDimFromIn} = labelsThisIn(subsCell{iDimFromIn}, :);
                    end
                    % then concatenate these columns together
                    labelsByDimOut{iDimOut} = cat(2, labelsThisOut{:});
                end
            end
            
        end
        
        function [out, newDims] = reshapeDimsInPlace(in, whichDims, newSizeInThoseDims)
            % reshapes a tensor by taking the consecutive dimension in
            % whichDims and making these a new size, e.g. tensorizes a
            % flattened dim or set of dims
            sz = size(in);
            assert(isequal(makecol(whichDims), makecol(min(whichDims):max(whichDims))), 'Dims must be consecutive');
            assert(prod(newSizeInThoseDims) == prod(TensorUtils.sizeMultiDim(in, whichDims)), 'Size must not change during reshape');
            mind = min(whichDims);
            maxd = max(whichDims);
            
            szNew = [sz(1:mind-1), newSizeInThoseDims, sz(maxd+1:end)];
            out = reshape(in, szNew);
            newDims = (mind : mind+numel(newSizeInThoseDims)-1)';
        end
        
        function out = flattenDimsInPlace(in, whichDims)
            newLen = prod(TensorUtils.sizeMultiDim(in, whichDims));
            out = TensorUtils.reshapeDimsInPlace(in, whichDims, newLen);
        end
        
        function out = selectSetsAlongDimension(in, dim, selectIdxCell)
            % given a tensor with size s1 x s2 x s3, where lets say dim is 2, 
            % selects multiple elements along dimension 2 for every
            % position in dimensions 1, 3. if selectIdxCell has length
            % S2, returns cell with size s1 x S2 x s3
            
            otherDims = TensorUtils.otherDims(size(in), dim);
            each = cellfun(@(idx) TensorUtils.squeezeSelectEachAlongDimension(...
                TensorUtils.selectAlongDimension(in, dim, idx), otherDims), selectIdxCell, ...
                'UniformOutput', false);
            out = cat(dim, each{:});
        end
        
        function out = selectSpecificIndicesAlongDimensionEachPosition(in, dim, idxForEachOtherDim)
            % typically we do a selection along dim e.g. where dim = 2, using in(:, idx, :) 
            % this function enables idx to vary with the position along the
            % other dimensions. If size(in) = in szIn, and
            % idxForEachOtherDim has szIdx, which matches szIn except on
            % dimension dim.
            % for example, if in = [1 2; 3 4], dim = 1, and
            % idxForEachOtherDim = [1 2], out is [1 4]
        
            out = cell2mat(TensorUtils.mapSlices(@(slice, idx) slice(idx), dim, in, idxForEachOtherDim));
        end
        
    end
    
    methods(Static) % Slice orienting and repmat
        % A slice is a selected region of a tensor, in which each dimension
        % either selects all of the elements (via :), or 1 of the elements.
        % This is the 2-d matrix equivalent of a row (slice with spanDim = 1)
        % or a column (slice with spanDim = 2). For a 3-d tensor, a slice
        % along spanDim = [2 3] would look like t(1, :, :);
        function out = orientSliceAlongDims(slice, spanDim)
            % given a slice with D dimensions, orient slice so that it's
            % dimension(i) becomes dimension spanDim(i)
            
            %szSlice = size(slice);
            ndimsSlice = length(spanDim);
            if ndimsSlice == 1
                % when spanDim is scalar, expect column vector or makecol
                slice = makecol(squeeze(slice));
            end
            
            ndimsOut = max(2, max(spanDim));
            sliceHigherDims = ndimsSlice+1:ndimsOut;
            nonSpanDims = setdiff(1:ndimsOut, spanDim);
            
            permuteOrder = nan(1, ndimsOut);
            permuteOrder(spanDim) = 1:ndimsSlice;
            permuteOrder(nonSpanDims) = sliceHigherDims;
            
            out = permute(slice, permuteOrder);
        end
        
        function out = orientSliceAlongSameDimsAs(slice, refSlice)
            % out = orientSliceAlongSameDimsAs(slice, refSlice)
            % orient slice such that it has the same shape as refSlice
            spanDim = find(size(refSlice) > 1);
            if ~isempty(spanDim);
                out = TensorUtils.orientSliceAlongDims(slice, spanDim);
            else
                out = slice;
            end
        end
        
        function out = repmatSliceAlongDims(slice, szOut, spanDim)
            % given a slice with ndims(slice) == numel(spanDim),
            % orient that slice along the dimensions in spanDim and repmat
            % it so that the output is size szOut.
            
            %nSpan = 1:length(spanDim);
            %ndimsOut = length(szOut);
            
            sliceOrient = TensorUtils.orientSliceAlongDims(slice, spanDim);
            
            repCounts = szOut;
            repCounts(spanDim) = 1;
            
            out = repmat(sliceOrient, repCounts);
        end
    end
    
    methods(Static) % Size expansion
        function out = singletonExpandToSize(in, szOut)
            % expand singleton dimensions to match szOut using repmat
            szIn = size(in);
            repCounts = szOut;
            repCounts(szIn ~= 1) = 1;
            out = repmat(in, repCounts);
        end
        
        function out = expandAlongDims(in, dims, by)
            % expand by adding nans or empty cells for cell array
            sz = size(in);
            sz(dims) = sz(dims) + by;
            
            out = TensorUtils.emptyWithSameType(in, sz);
            
            % build the mask over out to assign in into
            szInDims = TensorUtils.sizeMultiDim(in, dims);
            maskByDim = TensorUtils.maskByDimCellSelectPartialFromOrigin(size(out), dims, szInDims);
            
            out(maskByDim{:}) = in;
        end
    end
    
    methods(Static) % List resampling and shuffling along dimension
        function seedRandStream(seed)
            if nargin == 0
                seed = 'Shuffle';
            end
            s = RandStream('mt19937ar', 'Seed', seed);
            RandStream.setGlobalStream(s);
        end
        
        % Each of these methods apply to cell tensors containing vector
        % lists of values. They will move values around among the
        % different cells, i.e. from one list to another, but generally
        % they each preserve the total count within each list.
        % They also typically act along one or more dimensions, in that
        % they will move values among the cells that lie along a particular
        % dimension(s), but not among cells that lie at different positions
        % on the other dimensions. For example, if t is a 2-d matrix,
        % shuffling along dimension 1 would move values among cells that
        % lie along the same row, but not across rows. Here, the row would
        % be referred to as a slice of the matrix.
        %
        % Inputs are not required to have lists be column vectors, but the
        % output will have column vectors in each cell.
        %
        % The various functions differ in how they move the values around
        % and whether or not they sample with replacement, which would
        % allow some values to be replicated multiple times and others to
        % not appear in the final lists.
        
        function t = listShuffleAlongDimension(t, iA, replace)
            % t is a cell tensor of vectors. For each slice along iA, i.e.
            % t(i, j, ..., :, ..., k) where the : is in position iA,
            % shuffle the elements among all cells of that slice,
            % preserving the total number in each cell. If replace is true
            % shuffles with replacement.
            
            if nargin < 3
                replace = false;
            end
            
            t = TensorUtils.mapSlicesInPlace(@shuffleFn, iA, t);
            
            function sNew = shuffleFn(s)
                list = TensorUtils.combineListFromCells(s);
                if ~isscalar(list)
                    % randsample with scalar first input thinks we mean
                    % 1:scalar as the population rather than just [scalar]
                    list = randsample(list, numel(list), replace);
                end
                sNew = TensorUtils.splitListIntoCells(list, cellfun(@numel, s));
            end
        end
        
        function t = listResampleFromSame(t)
            % resample from replacement within each list, i.e. don't move
            % anything between lists, just take each list and resample with
            % replacement from itself.
            t = TensorUtils.map(@(list) randsample(list, numel(list), true), t);
        end
        
        function t = listResampleFromSpecifiedAlongDimension(t, from, iA, replace)
            % resample from replacement from other lists, according to a
            % specified set of cell elements along each slice. from{i} is a
            % specifies the linear index (or indices) of which cells to
            % sample from when building the list for the cell at position i
            % along that slice. If replace is true, samples with
            % replacement.
            %
            % Example: listResampleFromSpecified(t, {1, [1 2]}, 1) where t
            % is a 3 x 2 matrix. For each row i of t, the new list at cell
            % t{i, 1} will be built by sampling from list in t{i, 1}.
            % The new list at cell t{i, 2} will be built by sampling from
            % both t{i, 1} and t{i, 2}.
            
            if nargin < 4
                replace = true;
            end
            
            t = TensorUtils.mapSlicesInPlace(@resampleFn, iA, t);
            
            function sNew = resampleFn(s)
                sNew = cell(size(s));
                for i = 1:numel(s)
                    list = TensorUtils.combineListFromCells(s(from{i}));
                    sNew{i} = randsample(list, numel(s{i}), replace);
                end
            end
        end
        
        function list = combineListFromCells(t)
            % gathers all lists from t{:} into one long column vector
            t = cellfun(@makecol, t, 'UniformOutput', false);
            % make empty vectors column vectors to allow cell2mat to work
            [t{cellfun(@isempty, t)}] = deal(nan(0, 1));
            list = cell2mat(t(:));
        end
        
        function t = splitListIntoCells(list, nPerCell)
            % list is a vector list, nPerCell is a numeric tensor
            % specifying how many elements to place in each list.
            % This undoes combineListFromCells
            if isempty(list)
                % all cells will be empty, special case
                t = cellvec(numel(nPerCell));
            else
                t = mat2cell(makecol(list), nPerCell(:), 1);
            end
            t = reshape(t, size(nPerCell));
            
        end
        
        function c = splitAlongDimension(t, dim, nPerCell)
           % c = splitAlongDimension(t, dim, nPerCell)
           % wrapper around mat2cell that splits t along dimension t so
           % that size(c) == [numel(nPerCell) 1] and size(c{i}, dim) == nPerCell(i)
           % if dim is a vector, then nPerCell may be a cell vector which
           % will independently split along each dimension in dim
           
           args = cellvec(ndims(t));
           if ~iscell(nPerCell)
               nPerCell = {nPerCell};
           end
           
           assert(numel(nPerCell) == numel(dim), 'nPerCell must be vector if dim is scalar or cell with one element per dim');
           for iDim = 1:ndims(t)
               [tf, which] = ismember(iDim, dim);
               if tf
                   args{iDim} = nPerCell{which};
               else
                   args{iDim} = size(t, iDim);
               end
           end
           
           c = mat2cell(t, args{:});
        end
        
    end
   
    methods(Static) % Statistics       
        function t = meanMultiDim(t, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3]) 
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = cell2mat(TensorUtils.mapSlicesInPlace(@(slice) mean(slice(:)), dims, t));
        end
        
        function t = nanmeanMultiDim(t, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3]) 
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = cell2mat(TensorUtils.mapSlicesInPlace(@(slice) nanmean(slice(:)), dims, t));
        end
        
        function t = nansumMultiDim(t, dims)
            for iD = 1:numel(dims)
                t = nansum(t, dims(iD));
            end
        end
        
        function t = varMultiDim(t, dims, varargin)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3]) 
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = cell2mat(TensorUtils.mapSlicesInPlace(@(slice) var(slice(:), varargin{:}), dims, t));
        end
        
        function t = stdMultiDim(t, dims, varargin)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3]) 
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = cell2mat(TensorUtils.mapSlicesInPlace(@(slice) std(slice(:), varargin{:}), dims, t));
        end
    end
    
    methods(Static) % Data manipulation
        function t = centerSlicesOrthogonalToDimension(t, alongDims)
            % takes slices along dimensions specified, computes the mean 
            % of the slice and subtracts it. this ensures
            % that the mean of each slice spanning alongDims will be zero
            otherDims = TensorUtils.otherDims(size(t), alongDims);
            meanTensor =  TensorUtils.nanmeanMultiDim(t, otherDims);
            t = bsxfun(@minus, t, meanTensor);
        end
        
        
        function t = centerSlicesSpanningDimension(t, alongDims)
            % for each subscript in dimension(s) alongDims, computes the mean 
            % along all other dimensions and subtracts it. this ensures
            % that the mean along any slice in alongDims will have zero
            % mean.
            meanTensor =  TensorUtils.nanmeanMultiDim(t, alongDims);
            t = bsxfun(@minus, t, meanTensor);
        end
        
        function t = zscoreAlongDimension(t, alongDims)
            % for each subscript in dimension(s) alongDims, computes the mean 
            % along all other dimensions and subtracts it. this ensures
            % that the mean along any slice in alongDims will have zero
            % mean.
            t = TensorUtils.centerSlicesSpanningDimension(t, alongDims);
            stdTensor = TensorUtils.stdMultiDim(t, alongDims);
            t = bsxfun(@rdivide, t, stdTensor);
        end
        
        function t = makeNonDecreasing(t, dim)
            if nargin < 2
                dim = TensorUtils.firstNonSingletonDim(t);
            end
            
            % take cumulative max and replace each value that drops below
            % the cumulative max with the cumulative max at that point
            cm = cummax(t, dim);
            mask = t < cm;
            t(mask) = cm(mask);
        end
        
        function t = makeNonIncreasing(t, dim)
            if nargin < 2
                dim = TensorUtils.firstNonSingletonDim(t);
            end
            
            % take cumulative min and replace each value that rises above
            % the cumulative min with the cumulative min at that point
            cm = cummin(t, dim);
            mask = t > cm;
            t(mask) = cm(mask);
        end
        
        function reweightedTensor = linearCombinationAlongDimension(t, dim, weightsNewByOld, varargin)
            % looking along dimension dim, reweight the tensor along that
            % dimension by linearly combining. If dim==1 and the tensor
            % is a matrix, this is equivalent to matrix multiplication,
            % reweightedTensor = weightsNewByOld * t.
            
            p = inputParser();
            p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            p.addParameter('keepNaNIfAllNaNs', false, @islogical); % when replaceNaNWithZero is true, keep the result as NaN if every entry being combined is NaN
            % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis
            % this enables nanmean like computations
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical); 
            % requires that the weight matrix be square, the equivalent of adding
            % the identity matrix to the weight matrix, except that this
            % will be added after normalization. 
            p.addParameter('addToOriginal', false, @islogical); 
            p.parse(varargin{:});
            
            nOld = size(t, dim);
            assert(size(weightsNewByOld, 2) == nOld, 'Size of weight matrix must have nOld==%d columns', nOld);
            nNew = size(weightsNewByOld, 1);
            
            if p.Results.addToOriginal
                assert(size(weightsNewByOld, 1) == size(weightsNewByOld, 2), 'Weight matrix must be square for addToOriginal = true');
            end
            
            sz = size(t);
            newSz = sz;
            newSz(dim) = nNew;
            
            % put combination dimension first
            pdims = [dim, TensorUtils.otherDims(sz, dim)];
            tp = permute(t, pdims);
            
            % should be nOld x prod(size-t-other-dims)
            tpMat = tp(:, :);
            
            if p.Results.normalizeCoefficientsByNumNonNaN || p.Results.keepNaNIfAllNaNs
                % count the number of values in each row
                nValidMat = sum(~isnan(tpMat), 1);
            end
            
            if p.Results.replaceNaNWithZero
                tpMat(isnan(tpMat)) = 0;
            end
            
            % should be nNew x prod(size-t-other-dims)
            reweightMat = weightsNewByOld * tpMat;
            
            if p.Results.normalizeCoefficientsByNumNonNaN
                % do the reweighting
                reweightMat = bsxfun(@rdivide, reweightMat, nValidMat);
            end
            
            if p.Results.addToOriginal
                reweightMat = reweightMat + tpMat;
            end
            
            if p.Results.keepNaNIfAllNaNs
                % NaN out where nValidMat is zero
                mask = ones(size(nValidMat));
                mask(nValidMat == 0) = NaN;
                reweightMat = bsxfun(@times, reweightMat, mask);
            end
            
            reweightedTensor = ipermute(reshape(reweightMat, newSz(pdims)), pdims);         
        end
       
        function newTensor = linearCombinationApplyScalarFnAlongDimension(t, dim, logicalNewByOld, fn, varargin)
            % this method is conceptually similar to the linear combination
            % operation, except that it applies fn() to the set of values
            % along dimension dim that are non-zero in logicalNewByOld.
            
            nOld = size(t, dim);
            assert(size(logicalNewByOld, 2) == nOld, 'Size of weight matrix must have nOld==%d columns', nOld);
            nNew = size(logicalNewByOld, 1);
            
            % compute min over all trial counts included in each basis
            parts = cellvec(nNew);
            for iNew = 1:nNew
                fnSelect = @(c) fn(c(logicalNewByOld(iNew, :) ~= 0));
                parts{iNew} = TensorUtils.mapSlices(fnSelect, dim, t);
            end
            newTensor = cell2mat(cat(dim, parts{:}));  
        end
        
        function in = assignValueMaskedSelectionAlongDimension(in, dims, mask, value)
            % in = assignValueMaskedSelectionAlongDimension(in, dims, mask, value=NaN)
            % 
            % Assign value at each location where mask is true along dimension dim 
            % if dims is scalar assign in(..., mask, ...) = value
            % if dims is a vector of dimensions, then mask must be a cell,
            % in which case the assignment is
            % in(..., mask{1}, ..., mask{2}, ...) = value
            %
            
            if ~exist('value', 'var')
                if iscell(in)
                    value = {[]};
                elseif islogical(in)
                    value = false;
                else
                    value = NaN;
                end
            end
            if iscell(in) && ~iscell(value)
                value = {value};
            end
            
            assert(isscalar(value), 'Value must be scalar');
            maskByDim = TensorUtils.maskByDimCellSelectAlongDimension(size(in), dims, mask);
            in(maskByDim{:}) = value;
        end
    end
end
