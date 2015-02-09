classdef RangeNormalization < StateSpaceTranslationNormalization
    % normalize each basis to the range 0 to 1

    methods(Static)
        function obj = buildFromPopulationTrajectorySet(pset)
            % measure the range of each neuron
            rangeByBasis = pset.computeRangeByBasis();
            minByBasis = pset.computeMinByBasis();

            obj = StateSpaceTranslationNormalization.buildManual(...
                'translationByBasis', -minByBasis, 'translationByBasis', 'min-subtract', ...
                'normalizationDescription', rangeByBasis, 'normalizationByBasis', 'range-normalized');
        end
    end
    
end