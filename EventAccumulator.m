classdef EventAccumulator
    properties
        bins % left edge
        counts
        delta = 1;
    end
    
    properties(Dependent)
        total_count
        min
        max
        mean
        median
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
            
            C = size(data, 2);
            ea(C, 1) = EventAccumulator();
            
            binlo = floor(min(data(:) / delta)) * delta;
            binhi = floor(max(data(:) / delta)) * delta;
            bins = (binlo:delta:binhi)';
            edges = (binlo:delta:binhi+delta)';
            
            for c = 1:C
                ea(c).bins = bins;
                counts = histc(data(:, c), edges);
                counts(end-1) = counts(end-1) + counts(end); % include the right edge in the last bin
                ea(c).counts = counts(1:end-1);
            end
        end
        
        function v = get.total_count(ea)
            v = sum(ea.counts);
        end
        
        function v = get.min(ea)
            if isempty(ea.counts)
                v = NaN;
            else
                v = ea.bins(find(ea.counts > 0, 1, 'first'));
            end
        end
        
        function v = get.max(ea)
            if isempty(ea.counts)
                v = NaN;
            else
                v = ea.bins(find(ea.counts > 0, 1, 'last'));
            end
        end
        
        function v = get.mean(ea)
            if isempty(ea.counts)
                v = NaN;
            else
                v = sum(ea.counts .* ea.bins) / sum(ea.counts);
            end
        end
        
        function v = get.median(ea)
            if isempty(ea.counts)
                v = NaN;
            else
                cs = cumsum(ea.counts);
                v = ea.bins(find(cs >= 0.5*cs(end), 1, 'first'));
            end
        end   
        
        function out = kde(ea, varargin)
            out = EventAccumulator.build_kde(ea.bins, ea.counts, varargin{:});
        end
    end 
    
    methods(Static)
        function out = build_kde(bins, counts, bw)
            if numel(bins) < 2
                out = counts;
                return;
            end
            
            if nargin < 3 || isempty(bw)
                bw = (bins(2) - bins(1)) * 10;
            end
            
            data = cell2mat(arrayfun(@(b, c) repmat(b, c, 1), bins, counts, 'UniformOutput', false));
            out = ksdensity(data, bins, 'Bandwidth', bw);
        end
        
        function out = aggregate(varargin)
            % varargin is the cell over which to aggregate. Each can be a
            % matrix of EventAccumulators and the aggregation will happen
            % along that axis. If the sizes of these matrices don't match,
            % the result will be the largest size along each dimension
            
            args = varargin;
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
    
        