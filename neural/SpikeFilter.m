classdef SpikeFilter < handle & matlab.mixin.Copyable
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    % Dependent properties are provided as a convenience, override the underlying methods
    properties(Dependent)
        preWindowMs
        postWindowMs
        isCausal
    end

    methods(Abstract, Access=protected)
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        t = getPreWindowMs(sf)

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time 
        t = getPostWindowMs(sf)

        % spikeCell is nTrains x 1 cell array of time points
        % tvec = timeMin:timeSpacing:timeMax is an evenly spaced time vector
        rates = subclassFilterSpikeTrains(sf, spikeCell, tvec)
    end
    
    methods(Access=protected) % Subclasses may wish to override these 
        % return a string describing this filter's particular parameters,
        % do not include the classname as this will be added automatically
        function str = subclassGetDescription(sf) %#ok<MANU>
            str = '';
        end 
    end

    methods % Dependent properties
        function t = get.preWindowMs(sf)
            t = sf.getPreWindowMs();
        end

        function t = get.postWindowMs(sf)
            t = sf.getPostWindowMs();
        end

        function tf = get.isCausal(sf)
            tf = sf.getPreWindowMs() >= 0;
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
        function [rates, tvec] = filterSpikeTrains(sf, spikeCell, tWindow)
            rates = sf.subclassFilterSpikeTrains(spikeCell, tWindow);
            tvec = tWindow(1):tWindow(2);
        end
        
        function [rates, tvec] = filterSpikeTrainsWindowByTrial(sf, spikeCell, tMinByTrial, tMaxByTrial)
            tWindow(1) = nanmin(tMinByTrial);
            tWindow(2) = nanmax(tMaxByTrial);
            
            [rates, tvec] = sf.filterSpikeTrains(spikeCell, tWindow);
                
            % go through and mark as NaN any time outside each trial's
            % valid window, in case tMin / tMax are larger than that
            for iTrial = 1:numel(tMinByTrial)
                mask = falsevec(numel(tvec));
                mask(tvec < tMinByTrial(iTrial)) = true;
                mask(tvec > tMaxByTrial(iTrial)) = true;
                rates(iTrial, mask) = NaN;
            end
        end
    end
    
    methods(Static)
        function sf = getDefaultFilter()
            sf = GaussianSpikeFilter();
        end
    end
end
