classdef SpikeFilter < handle & matlab.mixin.Copyable
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    % Dependent properties are provided as a convenience, override the underlying methods
    properties(Dependent)
        preWindow
        postWindow
        isCausal
    end

    methods(Abstract, Access=protected)
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        t = getPreWindow(sf)

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        t = getPostWindow(sf)

        % spikeCell is nTrains x 1 cell array of time points
        %
        % tWindowPerTrial is nTrials x 2 matrix of tMin and tMax to grab
        % per trial, excluding any padding which should be accommodated for
        % when grabbing spikeCell
        %
        % multiplierToSpikesPerSec is the factor by which rates in the
        % spikeCell time units should be multiplied to get sec (e.g. 1000
        % for spikeCell in ms)
        [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowPerTrial, multiplierToSpikesPerSec)
    end
    
    methods(Access=protected) % Subclasses may wish to override these 
        % return a string describing this filter's particular parameters,
        % do not include the classname as this will be added automatically
        function str = subclassGetDescription(sf) %#ok<MANU>
            str = '';
        end 
    end

    methods % Dependent properties
        function t = get.preWindow(sf)
            t = sf.getPreWindow();
        end

        function t = get.postWindow(sf)
            t = sf.getPostWindow();
        end

        function tf = get.isCausal(sf)
            tf = sf.getPreWindow() >= 0;
        end
        
        function str = getDescription(sf)
            % build description using subclass implemented method including
            % causality and classname
            subclassStr = sf.subclassGetDescription();
            classname = class(sf);
            if sf.isCausal()
                causalStr = 'Causal';
            else
                causalStr = 'Acausal';
            end
            if isempty(subclassStr)
                str = sprintf('%s %s', causalStr, classname);
            else
                str = sprintf('%s %s(%s)', causalStr, classname, subclassStr);
            end
        end
    end

    methods
        function [rateCell, timeCell] = filterSpikeTrainsWindowByTrial(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec)
            % filters each trial individually, using a cell array to return
            % filtered rates and timeCell 
            tWindowMat = [makecol(tMinByTrial), makecol(tMaxByTrial)];
            [rateCell, timeCell] = sf.subclassFilterSpikeTrains(spikeCell, tWindowMat, multiplierToSpikesPerSec);
        end
        
        function [rates, tvec] = filterSpikeTrainsWindowByTrialAsMatrix(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, varargin)
            % filters each trial individually and then embeds each filtered
            % trace in a nTrials x nTime matrix, where missing samples are left as NaN
            % before and after each trial. tvec is the time vector that
            % indicates time along the columns.
            import TrialDataUtilities.Data.embedTimeseriesInMatrix;
            p = inputParser;
            p.addParamValue('timeDelta', 1, @isscalar);
            p.parse(varargin{:});
            
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec);
            
            % convert to matrix
            [rates, tvec] = embedTimeseriesInMatrix(rateCell, timeCell, ...
                'timeDelta', p.Results.timeDelta, 'interpolate', false, ...
                'fixDuplicateTimes', false);
        end
    end
    
    methods(Static)
        function sf = getDefaultFilter()
            sf = GaussianSpikeFilter();
        end
    end
end
