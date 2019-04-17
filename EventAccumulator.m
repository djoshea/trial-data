classdef EventAccumulator
    properties
        bins % left edge
        counts
        delta = 1;
    end

    properties(Dependent)
        % Statistics are computed assuming each value takes the left edge
        % of its bin
        has_samples
        total_count
        min
        max
        mean
        median
        std
        var
        values
    end

    methods
        function ea = EventAccumulator(data, delta)
            % data is a column vector or a matrix. Each column of data will
            % generate its own EventAccumulator
            if nargin < 1
                return;
            end
            if nargin < 2
                delta = 1;
            end

            ea.delta = delta;

            if isequal(data, [])
                return;
            end

%             if isempty(data)
%                 error('No data provided to EventAccumulator constructor');
%             end
            C = size(data, 2);
            ea(C, 1) = EventAccumulator();

            binlo = floor(min(data(:) / delta)) * delta;
            binhi = floor(max(data(:) / delta)) * delta;
            bins = (binlo:delta:binhi)';
            edges = (binlo:delta:binhi+delta)';

            for c = 1:C
                ea(c).delta = delta;
                ea(c).bins = bins;
                counts = histc(data(:, c), edges)';
                if size(counts, 2) > size(counts, 1)
                    counts = counts'; % needs to be column vector
                end
                if numel(counts) > 1
                    counts(end-1) = counts(end-1) + counts(end); % include the right edge in the last bin
                    ea(c).counts = counts(1:end-1);
                end
            end
        end

        function tf = get.has_samples(ea)
            tf = ~isempty(ea.counts) && sum(ea.counts) > 0;
        end

        function v = get.total_count(ea)
            v = sum(ea.counts);
        end

        function v = get.min(ea)
            if ~ea.has_samples
                v = NaN;
            else
                v = ea.bins(find(ea.counts > 0, 1, 'first'));
            end
        end

        function v = get.max(ea)
            if ~ea.has_samples
                v = NaN;
            else
                v = ea.bins(find(ea.counts > 0, 1, 'last'));
            end
        end

        function v = get.mean(ea)
            if ~ea.has_samples
                v = NaN;
            else
                v = sum(ea.counts .* ea.bins) / sum(ea.counts);
            end
        end

        function v = get.median(ea)
            if ~ea.has_samples
                v = NaN;
            else
                cs = cumsum(ea.counts);
                v = ea.bins(find(cs >= 0.5*cs(end), 1, 'first'));
            end
        end

        function v = get.std(ea)
            if ~ea.has_samples
                v = NaN;
            else
                mu = ea.mean;
                v = sqrt(sum(ea.counts .* (ea.bins - mu).^2) / (ea.total_count - 1));
            end
        end

        function v = get.var(ea)
            if ~ea.has_samples
                v = NaN;
            else
                mu = ea.mean;
                v = sum(ea.counts .* (ea.bins - mu).^2);
            end
        end

        function out = kde(ea, varargin)
            out = EventAccumulator.build_kde(ea.bins, ea.counts, varargin{:});
        end

        function values = get.values(ea)
            values = cell2mat(arrayfun(@(b, c) repmat(b, c, 1), ea.bins, ea.counts, 'UniformOutput', false));
        end

        function h = plotDistribution(ea, varargin)
            if ea.has_samples
                h = histogram(ea.values, ea.bins, varargin{:});
            end
        end

        function h = plotKDE(ea, bw, varargin)
            if nargin < 2
                bw = [];
            end

            density = ea.kde(bw);
            h = plot(ea.bins, density * ea.total_count, varargin{:});
        end

        function plotDistributionWithKDE(ea, bw)
            if nargin < 2
                bw = [];
            end

            tf = ishold();
            ea.plotDistribution('FaceColor', [0 0.4470 0.7410]);
            hold on;
            ea.plotKDE(bw, 'Color', 'k', 'LineWidth', 2);
            if ~tf, hold off; end
        end

    end

    methods(Static)
        function ea = constructForEachColumn(data, delta, aggregateOverColumnsAlso)
            % like the constructor, except supports empty data matrix (in
            % which case this returns [])

            % data is a column vector or a matrix.
            %
            % When aggregateOverColumnsAlso == true, each column of data will
            % generate its own EventAccumulator, so this will return a
            % column vector of EventAccumulator objects, one for each
            % column in data
            %
            % When aggregateOverColumnsAlso == false, only a scalar
            % EventAccumulator will be returned that flattens data(:).

            if nargin < 1
                return;
            end
            if nargin < 2
                delta = 1;
            end
            if nargin < 3
                aggregateOverColumnsAlso = false;
            end

            if isempty(data) && size(data, 2) == 0
                ea = [];
            elseif aggregateOverColumnsAlso
                ea = EventAccumulator(data(:), delta);
            else
                ea = EventAccumulator(data, delta);
            end
        end

        function out = build_kde(bins, counts, bw)
            if numel(bins) < 2
                out = counts;
                return;
            end

            if nargin < 3
                bw = [];
            end

            data = cell2mat(arrayfun(@(b, c) repmat(b, c, 1), bins, counts, 'UniformOutput', false));
            minmax = [min(data) max(data)];
            binMask = bins >= minmax(1) & bins <= minmax(2);
            if nnz(binMask) == 1
                out = counts;
            else
                if ~isempty(bw)
                    args = {'Bandwidth', bw};
                else
                    args = {};
                end
                out = ksdensity(data, bins, args{:});
            end
        end

        function out = aggregate(varargin)
            % varargin is the cell over which to aggregate. Each can be a
            % matrix of EventAccumulators and the aggregation will happen
            % among elements at the same location in each matrix, so that the result
            % will be a matrix of the same size.
            % If the sizes of these matrices don't match,
            % the result will be the largest size along each dimension

            emptyMask = cellfun(@isempty, varargin);
            args = varargin(~emptyMask);
            nA = numel(args);
            delta = cellfun(@(ag) ag(1).delta, args);
            assert(all(delta == delta(1)), 'Deltas must match');

            rows = cellfun(@(ag) size(ag, 1), args);
            cols = cellfun(@(ag) size(ag, 2), args);
            sz(1) = max(rows);
            sz(2) = max(cols);
            out(sz(1), sz(2)) = EventAccumulator();

            for i = 1:sz(1)
                for j = 1:sz(2)
                    % only look at args that are at least matrices with
                    % size i x j
                    argMask = rows >= i & cols >= j;
                    binlo = cellfun(@(ag) ag(i,j).bins(1), args(argMask), 'ErrorHandler', @(varargin) NaN);
                    binhi = cellfun(@(ag) ag(i,j).bins(end), args(argMask), 'ErrorHandler', @(varargin) NaN);
                    bins = (min(binlo):delta:max(binhi))';

                    counts = zeros(size(bins));
                    for iA = 1:nA
                        if argMask(iA) && ~isempty(args{iA}(i,j).counts)
                            mask = bins >= args{iA}(i,j).bins(1) & bins <= args{iA}(i,j).bins(end);
                            counts(mask) = counts(mask) + args{iA}(i,j).counts;
                        end
                    end

                    out(i,j) = EventAccumulator();
                    out(i,j).bins = bins;
                    out(i,j).counts = counts;
                    out(i,j).delta = delta(1);
                end
            end
        end

        function v = min_array(eav)
            v = arrayfun(@(ea) ea.min, eav);
        end

        function v = max_array(eav)
            v = arrayfun(@(ea) ea.max, eav);
        end

        function v = mean_array(eav)
            v = arrayfun(@(ea) ea.mean, eav);
        end

        function v = median_array(eav)
            v = arrayfun(@(ea) ea.median, eav);
        end
    end

end
