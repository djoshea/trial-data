classdef SpikeFilter % < handle & matlab.mixin.Copyable
% SpikeFilters take spike trains and provide rate estimates
% They also provide information about the amount of pre and post window timepoints 
% they require in order to estimate the rate at a given time point

    properties
        timeDelta = 1; % sampling interval for the filtered output
        binAlignmentMode = BinAlignmentMode.Centered;
        resampleMethod = 'filter';
    end
    
    % Dependent properties are provided as a convenience, override the underlying methods
    properties(Dependent)
        preWindow
        postWindow
        
        padWindow % includes pre / post window plus binWidth time for the edge bins
        isCausal
    end

    methods(Abstract, Access=protected) 
        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time, including both the spike filtering and
        % the timeDelta binning
        t = getPreWindow(sf)

        % return the time window of preceding spike data in ms required to estimate
        % the rate at any particular time, including both the spike filtering and
        % the timeDelta binning
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
        % 
        % this subclass method should respect sf.timeDelta,
        % sf.binAlignmentMode, and sf.resampleMethod. We could do this work
        % here, but not doing so leaves flexibility to the subclass
        [rateCell, timeCell] = subclassFilterSpikeTrains(sf, spikeCell, tWindowPerTrial, multiplierToSpikesPerSec)
    end
    
    methods(Access=protected) % Subclasses may wish to override these 
        % return a string describing this filter's particular parameters,
        % do not include the classname as this will be added automatically
        function str = subclassGetDescription(sf) %#ok<MANU>
            str = '';
        end
        
        function tf = getIsCausal(sf) % allows subclasses to override
            tf = sf.getPreWindow() >= 0 && sf.binAlignmentMode == BinAlignmentMode.Causal;
        end
    end

    methods % Dependent properties
        function checkOkay(sf)
            % throw an error if the settings are wrong
            return;
        end
        
        function t = get.preWindow(sf)
            t = sf.getPreWindow();
        end

        function t = get.postWindow(sf)
            t = sf.getPostWindow();
        end
        
        function w = get.padWindow(sf)
            w = [sf.preWindow + sf.binAlignmentMode.getBinStartOffsetForBinWidth(sf.binWidthMs), ...
                sf.postWindow + sf.binAlignmentMode.getBinStopOffsetForBinWidth(sf.binWidthMs)];
        end
                
        function tf = get.isCausal(sf)
            tf = sf.getIsCausal();
        end
        
        % these allow subclasses to override the behavior
        function sf = set.timeDelta(sf, v)
            sf.timeDelta = v;
            sf = sf.postSetTimeDelta();
        end
        
        function sf = postSetTimeDelta(sf)
        end
        
        % these allow subclasses to override the behavior
        function sf = set.binAlignmentMode(sf, v)
            sf.binAlignmentMode = v;
            sf = sf.postSetBinAlignmentMode();
        end
        
        function sf = postSetBinAlignmentMode(sf)
        end
        
        function sf = set.resampleMethod(sf, v)
            sf.resampleMethod = v;
            sf = sf.postSetResampleMethod();
        end
        
        function sf = postSetResampleMethod(sf)
            
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
                str = sprintf('%s %s (%s)', causalStr, classname, subclassStr);
            end
        end
        
        function str = char(sf)
            str = sf.getDescription();
        end
    end

    methods
        function sf = SpikeFilter(varargin)
            p = inputParser();
            p.addParameter('timeDelta', 1, @isscalar);
            p.addParameter('binAlignmentMode', BinAlignmentMode.Centered, @(x) isa(x, 'BinAlignmentMode'));
            p.addParameter('resampleMethod', 'filter', @isscalar);
            p.parse(varargin{:});
            
            sf.timeDelta = p.Results.timeDelta; % sampling interval for the filtered output
            sf.binAlignmentMode = p.Results.binAlignmentMode;
            sf.resampleMethod = p.Results.resampleMethod;
        end
       
        function [rateCell, timeCell] = filterSpikeTrainsWindowByTrial(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, varargin)
            % filters each trial individually, using a cell array to return
            % filtered rates and timeCell 
            sf.checkOkay();
            tWindowMat = [makecol(tMinByTrial), makecol(tMaxByTrial)];
            [rateCell, timeCell] = sf.subclassFilterSpikeTrains(spikeCell, tWindowMat, multiplierToSpikesPerSec, varargin{:});
        end
        
        function [rates, tvec] = filterSpikeTrainsWindowByTrialAsMatrix(sf, spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec, varargin)
            % filters each trial individually and then embeds each filtered
            % trace in a nTrials x nTime x nUnits matrix, where missing samples are left as NaN
            % before and after each trial. tvec is the time vector that
            % indicates time along the columns.=
            p = inputParser;
            p.addParameter('timeReference', 0, @isscalar);
            p.addParameter('interpolateMethod', 'linear', @ischar);
            p.addParameter('showProgress', true, @islogical);    
            p.addParameter('tMinByTrialExcludingPadding', [], @isvector);
            p.addParameter('tMaxByTrialExcludingPadding', [], @isvector);
            p.parse(varargin{:});
            
            % calls checkOkay
            [rateCell, timeCell] = sf.filterSpikeTrainsWindowByTrial(spikeCell, tMinByTrial, tMaxByTrial, multiplierToSpikesPerSec);
            
            % convert to matrix
            [rates, tvec] = TrialDataUtilities.Data.embedTimeseriesInMatrix(rateCell, timeCell, ...
                'assumeUniformSampling', true, ...
                'fixDuplicateTimes', false, ... % no need for this since we know the time vectors are monotonic
                'tMinExcludingPadding', p.Results.tMinByTrialExcludingPadding, ...
                'tMaxExcludingPadding', p.Results.tMaxByTrialExcludingPadding); 
        end
    end
    
    methods(Static)
        function sf = getDefaultFilter()
            sf = GaussianSpikeFilter();
        end
    end
end
