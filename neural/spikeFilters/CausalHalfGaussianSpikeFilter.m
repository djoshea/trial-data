classdef CausalHalfGaussianSpikeFilter < GaussianSpikeFilter
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    methods
        function sf = CausalHalfGaussianSpikeFilter(varargin)
            sf = sf@GaussianSpikeFilter('truncateFuture', 0, 'delayPeak', 0, varargin{:});
        end
    end
    
    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f sigma', sf.sigma);
        end
    end

end
