classdef BinAlignmentMode < int32
    enumeration
        Centered (1)
        Causal (2) % align bin so that all spikes came before the timepoint
        Acausal (3) % align bin so that all spikes came before the timepoint
    end
   
    methods
        function offsetStart = getBinStartOffsetForBinWidth(mode, binWidth)
            switch mode
                case BinAlignmentMode.Centered
                    offsetStart = -binWidth/2;
                case BinAlignmentMode.Causal
                    offsetStart = -binWidth;
                case BinAlignmentMode.Acausal
                    offsetStart = 0;
            end
        end
        
        function offsetStop = getBinStopOffsetForBinWidth(mode, binWidth)
            offsetStop = mode.getBinStartOffsetForBinWidth(binWidth) + binWidth;
        end
        
        function [tMinNew, tMaxNew] = getTimeLimitsForRebinning(mode, tMinOrig, tMaxOrig, origDelta, newDelta, timeReference)
            if nargin < 6
                timeReference = 0;
            end
            
            switch mode
                % the assumption here is that we treat the existing time
                % bins (with width origDelta) as having the same
                % binAlignmentMode as mode.
                case BinAlignmentMode.Causal
                    tMinNew = timeReference + ceil((tMinOrig - origDelta - timeReference + newDelta) ./ newDelta) .* newDelta;
                    tMaxNew = timeReference + floor((tMaxOrig - timeReference) / newDelta) * newDelta;
                case BinAlignmentMode.Acausal
                    tMinNew = timeReference + ceil((tMinOrig - timeReference) / newDelta) .* newDelta;
                    tMaxNew = timeReference + floor((tMaxOrig + origDelta - timeReference - newDelta) ./ newDelta) .* newDelta;
                case BinAlignmentMode.Centered
                    tMinNew = timeReference + ceil((tMinOrig - origDelta/2 - timeReference + newDelta/2) ./ newDelta) .* newDelta;
                    tMaxNew = timeReference + floor((tMaxOrig + origDelta/2 - timeReference - newDelta/2) ./ newDelta) .* newDelta;
            end
                    
            mask = tMinNew > tMaxNew;
            tMinNew(mask) = NaN;
            tMaxNew(mask) = NaN;
        end
        
        function [tvecLabelCell, tbinsForHistcCell] = generateMultipleBinnedTimeVectors(mode, starts, stops, binWidth)
            % tvecLabel will be the time that each time bin is labeled
            % with, tbinsForHistc is the argument to histc that will select
            % times in the bin so that the output of histc matches
            % tvecLabel. The last bin output by histc should be thrown
            % away.
            %
            % actual spike zero/stop/zero will be at 0/start/stop + binAlignmentMode.getBinStartOffsetForBinWidth
                
            assert(numel(starts) == numel(stops));
            
            [tvecLabelCell, tbinsForHistcCell] = cellvec(numel(starts));
            
            binOffset = mode.getBinStartOffsetForBinWidth(binWidth);
            
            [tMin, tMax] = mode.getTimeLimitsForRebinning(starts, stops, binWidth, binWidth, 0);
            
            for i = 1:numel(starts)
                if isnan(tMin(i)) || isnan(tMax(i))
                    continue;
                end
                
                tvecLabel = (tMin(i):binWidth:tMax(i))';
                if isempty(tvecLabel)
                    continue;
                end
                
                % remove the offset from the time labels
                tvecLabelCell{i} = tvecLabel;
                
                % add the end of the last bin onto the histc bins
                tbinsForHistcCell{i} = [tvecLabel+binOffset; tvecLabel(end)+binOffset+binWidth];
            end
        end
        
        function [tvecLabel, tbinsForHistc, tbinsValidMat] = generateCommonBinnedTimeVector(mode, starts, stops, binWidth)
            % note that tbinsForHistc will have one additional entry than
            % tbinsValidMat and tvecLabel. This is for histc. Since histc returns 
            % an additional bin which captures values equal to the final
            % value, you should throw away the final bin returned by histc.
            start = nanmin(starts);
            stop = nanmax(stops);
            [tvecLabelCell, tbinsForHistcCell] = mode.generateMultipleBinnedTimeVectors(start, stop, binWidth);
            tvecLabel = tvecLabelCell{1};
            tbinsForHistc = tbinsForHistcCell{1};
            
            nTrials = numel(starts);
            tbinsValidMat = false(nTrials, numel(tvecLabel));
            for iTrial = 1:nTrials
                tbinsValidMat(iTrial, :) = tbinsForHistc(1:end-1) >= starts(iTrial) & tbinsForHistc(2:end) <= stops(iTrial);
            end
        end
    end
    
end