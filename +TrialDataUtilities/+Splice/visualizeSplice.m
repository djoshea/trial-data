function visualizeSplice(dataPre, dataPost, dataSpliced, varargin)

    p = inputParser();
    p.addParameter('nOverlap', 0, @(x) isvector(x) || isscalar(x));
    p.addParameter('ndim', 1, @isscalar);
    p.addParameter('dims', 1:6, @isvector);
    p.addParameter('conditionIdx', 1, @isvector);
    p.addParameter('nPCs', 6, @isscalar);
    p.parse(varargin{:});

    [coeff, ~, ~, ~, ~, mu] =  TensorUtils.pcaAlongDim(dataSpliced,  1, 'NumComponents', p.Results.nPCs);
    
    pre = TensorUtils.linearCombinationAlongDimension(dataPre-mu, 1, coeff', 'replaceNaNWithZero', true);
    post = TensorUtils.linearCombinationAlongDimension(dataPost-mu, 1, coeff', 'replaceNaNWithZero', true);
    spliced = TensorUtils.linearCombinationAlongDimension(dataSpliced-mu, 1, coeff', 'replaceNaNWithZero', true);

    ndim = p.Results.ndim;
    ms = 20;
    dims = p.Results.dims;
    dims = dims(dims < p.Results.nPCs);
    nOverlap = p.Results.nOverlap;
    conditionIdx = p.Results.conditionIdx;
    
    postColor = 'b';
    
    clf;
    if ndim == 2 || ndim == 3
        for iC = 1:numel(conditionIdx)
            c = conditionIdx(iC);

            if ndim == 3
                plot3(pre(dims(1), :, c), pre(dims(2), :, c), pre(dims(3), :, c), 'k-');
                hold on;
                plot3(pre(dims(1), end-floor(nOverlap/2), c), pre(dims(2), end-floor(nOverlap/2), c), pre(dims(3), end-floor(nOverlap/2), c), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'k');


                plot3(post(dims(1), :, c), post(dims(2), :, c), post(dims(3), :, c), '-', 'Color', postColor);
                plot3(post(dims(1), 1+floor(nOverlap/2), c), post(dims(2), 1+floor(nOverlap/2), c), post(dims(3), 1+floor(nOverlap/2), c), 'o', 'MarkerSize', ms, 'MarkerFaceColor', postColor);

                % plot spliced
                plot3(spliced(dims(1), :, c), spliced(dims(2), :, c), spliced(dims(3), :, c), '--', 'Color', 'r');
            elseif ndim == 2
                ms = 20;
                plot(pre(dims(1), :, c), pre(dims(2), :, c), 'k-');
                hold on;
                plot(pre(dims(1), end-floor(nOverlap/2), c), pre(dims(2), end-floor(nOverlap/2), c), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'k');


                plot(post(dims(1), :, c), post(dims(2), :, c), '-', 'Color', postColor);
                plot(post(dims(1), 1+floor(nOverlap/2), c), post(dims(2), 1+floor(nOverlap/2), c), 'o', 'MarkerSize', ms, 'MarkerFaceColor', postColor);

                % plot spliced
                plot(spliced(1, :, c), spliced(dims(2), :, c), '--', 'Color', 'r');
            end
        end
    else
        if nargin < 4
            nOverlap = 0;
        end
        nPre = size(pre, 2);
        nPost = size(post, 2);
        nSpliced = size(spliced, 2);

        sz = size(spliced);
        sz(3) = numel(conditionIdx);
        if sz(2) < nPre + nPost - nOverlap
            sz(2) = nPre + nPost - nOverlap;
        end
        inflatePre = nan(sz);
        inflatePre(:, 1:nPre, conditionIdx) = pre(:, :, conditionIdx);
        inflatePost = nan(sz);
        
        for iC = 1:numel(conditionIdx)
            c = conditionIdx(iC);
            inflatePost(:, (1:nPost) + nPre - nOverlap,  c) = post(:, :, c);
        end
        
        inflateSpliced = nan(sz);
        inflateSpliced(:, 1:nSpliced, conditionIdx) = spliced(:, :, conditionIdx);

        inflateAll = cat(3, inflatePre, inflatePost, inflateSpliced);
        ptstack(2, 1, inflateAll, 'colormap', {'k', postColor, 'r'});
    end
end