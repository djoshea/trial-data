classdef SoftRangeNormalization < StateSpaceTranslationNormalization
    % normalize by range + K
    
    methods(Access=protected) % don't call constructor, use factory methods
        function obj = SoftRangeNormalization()
            % should be built using factory methods below
        end
    end
    
    methods(Static)
        function obj = buildFromPopulationTrajectorySet(pset, K)
            % measure the range of each neuron
            rangeByBasis = pset.computeRangeByBasis();
            minByBasis = pset.computeMinByBasis();
            
            normalization = rangeByBasis + K;
            
%             % don't touch bases whose range is less than M
%             rangeByBasis(maskLow) = 1;
%             % normalize everything else down to M
%             rangeByBasis(~maskLow) = rangeByBasis(~maskLow) / M;            

            obj = StateSpaceTranslationNormalization.buildManual(...
                -minByBasis, normalization, ...
                'translationDescription', 'min-subtract', ...
                'normalizationDescription', 'soft-range-normalized');
        end
    end
    
end