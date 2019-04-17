classdef RepeatedUnitUtils
    methods(Static)
        function [chList, electrodeNums] = getChListByElectrodeForArray(td, array)
            aeTable = td.listArrayElectrodesWithSpikeChannelsAsTable();
            mask = strcmp(aeTable.array, array);
            aeTable = sortrows(aeTable(mask, :), {'electrode'});
            chList = aeTable.channelList;
            electrodeNums = aeTable.electrode;
        end 

        function [xc, edges] = computeCrossCorrelation(td, unitA, unitB, varargin)
            p = inputParser();
            p.addParameter('edges', 0:0.05:10, @isvector);
            p.addParameter('showPlot', false, @islogical);
            p.addParameter('normalize', false, @islogical);
            p.parse(varargin{:});

            timesA = td.getSpikeTimes(unitA, 'combine', true);
            timesB = td.getSpikeTimes(unitB, 'combine', true);

            edges = p.Results.edges;

            xc = zeros(numel(edges)-1, 1);
            for t = 1:td.nTrials
                tA = timesA{t};
                tB = timesB{t};

                if isempty(tB) || isempty(tA), continue; end
                dist = pdist2(tA, tB, 'cityblock');

                xc = xc + histcounts(dist(:), edges)';
            end
            
            if p.Results.normalize
                xc = xc ./ nansum(xc);
            end
            
            if p.Results.showPlot
                stairs(edges(1:end-1), xc, 'LineWidth', 1);
                xlabel('{\Delta}Time (ms)');
                ylabel('Density');
                AutoAxis.replace();
            end
        end

        function xcRatioMat = computeCrossCorrPeakRatioMatrix(td, unitNameSets, varargin)
            p = inputParser;
            p.addParameter('peakWidth', 0.2, @isscalar);
            p.addParameter('offPeak', [10 40], @isscalar);
            p.addParameter('softNormSpikes', 100, @isscalar);
            p.parse(varargin{:});

            edges = [0 p.Results.peakWidth p.Results.offPeak];
            span = [p.Results.peakWidth, diff(p.Results.offPeak)];

            C = numel(unitNameSets);
            xcRatioMat = nan(C, C);
            prog = ProgressBar(C*(C-1)/2, 'Computing cross-correlation matrix for %d electrodes', C);
            for a = 1:numel(unitNameSets)
                for b = (a+1):numel(unitNameSets)
                    prog.increment('Computing cross-correlation matrix for %s vs. %s', strjoin(unitNameSets{a}), strjoin(unitNameSets{b}));
                    xc = TrialDataUtilities.SpikeData.RepeatedUnitUtils.computeCrossCorrelation(...
                        td, unitNameSets{a}, unitNameSets{b}, 'edges', edges);
                    ratio = xc(1)/span(1) / (xc(3)/span(2) + p.Results.softNormSpikes);

                    xcRatioMat(a,b) = ratio;
                    xcRatioMat(b,a) = ratio;
               end
            end
            prog.finish();
        end
        
        function xcRatioMat = computeCrossCorrPeakRatioMatrixForArray(td, array, varargin)
            p = inputParser();
            p.addParameter('showPlot', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [chList, electrodeNums] = TrialDataUtilities.SpikeData.RepeatedUnitUtils.getChListByElectrodeForArray(td, array);
            xcRatioMatReduced = TrialDataUtilities.SpikeData.RepeatedUnitUtils.computeCrossCorrPeakRatioMatrix(td, chList, p.Unmatched);
            
            N = max(electrodeNums);
            xcRatioMat = TensorUtils.inflateMaskedTensor(xcRatioMatReduced, [1 2], electrodeNums);
            
            if p.Results.showPlot
                pmat(xcRatioMat);
                xlabel('Electrode #');
                ylabel('Electrode #');
                title('Peak to non-peak cross-corr ratio');
            end
        end
        
        function [td, nRemovedMat] = removeSharedSpikes(td, setWillKeep, setWillRemove, varargin)
            % remove spikes that occur in channels in setWillRemove if they also occur in setWillKeep
            % within a window of +/- timeWindow/2 ms
            p = inputParser;
            p.addParameter('timeWindow', 0.5, @isscalar);
            p.addParameter('keepRemovedSpikesAs', '', @ischar);
            p.parse(varargin{:});
            
            setWillKeep = TrialDataUtilities.Data.wrapCell(setWillKeep);
            setWillRemove = TrialDataUtilities.Data.wrapCell(setWillRemove);
            
            nR = numel(setWillRemove);

            % nTrials x 1
            dataKeep = td.getSpikeTimes(setWillKeep, 'combine', true);
            timeWindowHalf = p.Results.timeWindow / 2;
            
            % nTrials x nR
            dataRemove = td.getSpikeTimes(setWillRemove);
            maskKeep = cellfun(@(times) true(size(times)), dataRemove, 'UniformOutput', false);
            
            for c = 1:nR
                unitName = setWillRemove{c};
                prog = ProgressBar(td.nTrials, 'Removing shared spikes by trial');
                for t = 1:td.nTrials
                    prog.update(t);
                    thisTrial = dataRemove{t, c};

                    if ~isempty(thisTrial)
                        thisTrialOC = dataKeep{t};

                        % detect my matches on this other channel, i.e. whether i have
                        % a matching spike within timeWindow/2 of each of my spikes
                        if ~isempty(thisTrialOC)
                            minDist = pdist2(thisTrialOC, thisTrial, 'cityblock', 'Smallest', 1);
                            maskHasMatch = minDist' <= timeWindowHalf;
                            
                            % mark spikes with close match for removal
                            maskKeep{t, c} = ~maskHasMatch;
                        end
                    end
                end
                prog.finish();

                td = td.maskSpikeChannelSpikes(unitName, maskKeep(:, c), 'keepRemovedSpikesAs', p.Results.keepRemovedSpikesAs);    
            end
            
            nRemovedMat = cellfun(@(x) nnz(~x), maskKeep);
        end
        
        function [td, nRemoved] = removeSharedSpikesUsingAdjacencyMatrix(td, chNames, adjacencyMatrix, varargin)
            % chNames is nChannels x 1 cellstr of spike channel names in td
            % adjacencyMatrix should be constructed so that spikes on
            % channel i will be deleted if they are coincident with spikes
            % on channel j if adjacencyMatrix(i, j) == true. Thus a lower triangular matrix will preserve
            % shared spikes on the channel with the lower index but delete
            % them from the channel with the higher index.
            %
            % nRemoved is nTrials x nChannels
            
            import TrialDataUtilities.SpikeData.RepeatedUnitUtils;
            
            p = inputParser;
            p.addParameter('timeWindow', 0.1, @isscalar);
            p.addParameter('electrodeSpan', 1, @isscalar);
            p.addParameter('keepRemovedSpikes', false, @islogical);
            p.parse(varargin{:});
            
            assert(all(ismember(chNames, td.listSpikeChannels())), 'Some spike channels were not found');
           
            N = numel(chNames);
            nRemoved = zeros(td.nTrials, N);
            
            assert(size(adjacencyMatrix, 1) == N && size(adjacencyMatrix, 2) == N, 'Size of adjacency matrix must be nChannels x nChannels');
            adjacencyMatrix = logical(adjacencyMatrix);
            for c = 1:N
                adjacencyMatrix(c, c) = false; % set diagonal false
            end
            
            % loop from max electrode number to min, remove spikes from the
            % max of the two, so that the effects cascade up the 
            prog = ProgressBar(N, 'Removing shared spikes');
            for cRemoveFrom = 1:N
                chRemoveFrom = chNames{cRemoveFrom};
                
                % we'll remove spike from cRemoveFrom that are coincident
                % with spikes on these channels
                refChannels = chNames(adjacencyMatrix(cRemoveFrom, :))';

                if isempty(refChannels)
                    continue;
                end
                prog.increment('Removing spikes from electrode %s', chRemoveFrom);
                    
                if p.Results.keepRemovedSpikes
                    [a,e] = SpikeChannelDescriptor.parseArrayElectrodeUnit(chRemoveFrom);
                    keepAs = SpikeChannelDescriptor.generateNameFromArrayElectrodeUnit(a, e, 255);
                else
                    keepAs = '';
                end
                
                % keep in low, remove from high
                [td, nRemoved(:, cRemoveFrom)] = RepeatedUnitUtils.removeSharedSpikes(td, ...
                    refChannels, chRemoveFrom, 'timeWindow', p.Results.timeWindow, ...
                    'keepRemovedSpikesAs', keepAs);
            end
            prog.finish();
        end
        
        function [td, nRemoved] = removeSharedSpikesAdjacentElectrodesOnArray(td, array, varargin)
            % remove spikes that occur in channels in setWillRemove if they also occur in setWillKeep
            % within a window of +/- timeWindow/2 ms
            import TrialDataUtilities.SpikeData.RepeatedUnitUtils;
            
            p = inputParser;
            p.addParameter('timeWindow', 0.1, @isscalar);
            p.addParameter('electrodeSpan', 1, @isscalar);
            p.addParameter('keepRemovedSpikes', false, @islogical);
            p.parse(varargin{:});
            
            [chList, electrodeNums] = RepeatedUnitUtils.getChListByElectrodeForArray(td, array);
            
            N = max(electrodeNums);
            nRemoved = zeros(N, 1);
            
            % loop from max electrode number to min, remove spikes from the
            % max of the two, so that the effects cascade up the 
            span = min(p.Results.electrodeSpan, N);
            prog = ProgressBar(N-1, 'Removing shared spikes');
            for cHigh = N:-1:2
                if ~ismember(cHigh, electrodeNums), continue; end
                chListHigh = chList{electrodeNums == cHigh};

                % concatenate all chListLow
                chListLow = cellvec(span);
                for delta = 1:span
                    cLow = cHigh-delta;
                    if ~ismember(cLow, electrodeNums), continue; end
                    chListLow{delta} = chList{electrodeNums == cLow};
                end
                chListLow = cat(1, chListLow{:});
                if isempty(chListLow), continue, end
                cLowList = intersect(cHigh-span:cHigh-1, electrodeNums);
                
                prog.increment('Removing spikes from electrode %d also in %s', cHigh, TrialDataUtilities.String.strjoin(cLowList));
                    
                if p.Results.keepRemovedSpikes
                    keepAs = sprintf('%s%02d_255', array, cHigh);
                else
                    keepAs = '';
                end
                
                % keep in low, remove from high
                [td, nRemovedMat] = RepeatedUnitUtils.removeSharedSpikes(td, ...
                    chListLow, chListHigh, 'timeWindow', p.Results.timeWindow, ...
                    'keepRemovedSpikesAs', keepAs);
                
                nRemoved(cHigh) = nRemoved(cHigh) + sum(nRemovedMat(:));
            end
            prog.finish();
        end     
    end
end