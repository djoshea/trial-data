classdef DelayedCausalGaussianSpikeFilter < GaussianSpikeFilter

    methods
        function sf = DelayedCausalGaussianSpikeFilter(varargin)
            sf = sf@GaussianSpikeFilter(varargin{:});
            sf.binAlignmentMode = BinAlignmentMode.Causal;
            sf.delayPeak = sf.halfWidthSigmas * sf.sigma;
            sf.truncateFuture = 0;
        end
    end
    
    methods(Access=protected)
        function str = subclassGetDescription(sf)
            str = sprintf('%.0f sigma', sf.sigma);
        end
    end

end
