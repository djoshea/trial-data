classdef SpikeBinAlignmentMode < int32
    enumeration
        Centered (1)
        Causal (2) % align bin so that all spikes came before the timepoint
        Acausal (3) % align bin so that all spikes came before the timepoint
    end
   
    methods
        function offsetStart = getBinStartOffsetForBinWidth(mode, binWidth)
            switch mode
                case SpikeBinAlignmentMode.Centered
                    offsetStart = -binWidth/2;
                case SpikeBinAlignmentMode.Causal
                    offsetStart = -binWidth;
                case SpikeBinAlignmentMode.Acausal;
                    offsetStart = 0;
            end
        end
        
        function offsetStop = getBinStopOffsetForBinWidth(mode, binWidth)
            offsetStop = mode.getBinStartOffsetForBinWidth(binWidth) + binWidth;
        end
    end
    
    methods(Static)
        function [tvecLabel, tbinsForHistc, tbinsValidMat] = generateCommonBinnedTimeVector(starts, stops, binWidth, binAlignmentMode)
            % note that tbinsForHistc will have one additional entry than
            % tbinsValidMat and tvecLabel. This is for histc. Since histc returns 
            % an additional bin which captures values equal to the final
            % value, you should throw away the final bin returned by histc.
            start = nanmin(starts);
            stop = nanmax(stops);
            [tvecLabelCell, tbinsForHistcCell] = SpikeBinAlignmentMode.generateMultipleBinnedTimeVectors(start, stop, binWidth, binAlignmentMode);
            tvecLabel = tvecLabelCell{1};
            tbinsForHistc = tbinsForHistcCell{1};
            
            nTrials = numel(starts);
            tbinsValidMat = false(nTrials, numel(tvecLabel));
            for iTrial = 1:nTrials
                tbinsValidMat(iTrial, :) = tbinsForHistc(1:end-1) >= starts(iTrial) & tbinsForHistc(2:end) <= stops(iTrial);
            end
        end
        
        function [tvecLabelCell, tbinsForHistcCell] = generateMultipleBinnedTimeVectors(starts, stops, binWidth, binAlignmentMode)
            % tvecLabel will be the time that each time bin is labeled
            % with, tbinsForHistc is the argument to histc that will select
            % times in the bin so that the output of histc matches
            % tvecLabel. The last bin output by histc should be thrown
            % away.
            %
            % actual spike zero/stop/zero will be at 0/start/stop + binAlignmentMode.getBinStartOffsetForBinWidth
                
            assert(numel(starts) == numel(stops));
            
            [tvecLabelCell, tbinsForHistcCell] = cellvec(numel(starts));
            
            binOffset = binAlignmentMode.getBinStartOffsetForBinWidth(binWidth);
            
            nPre = ceil(-starts/binWidth) + 1;
            nPost = ceil(stops/binWidth) + 1;
            
            for i = 1:numel(starts)
                if isnan(nPre(i)) || isnan(nPost(i))
                    continue;
                end
                
                tvecRaw = (-nPre(i) : nPost(i))' * binWidth + binOffset;
                if isempty(tvecRaw)
                    continue;
                end
                
                % we'll now filter tvecRaw for bins that fall entirely within
                % start : stop. tvecRaw will specify the start of the bin
                tvecRaw = tvecRaw(tvecRaw+binOffset >= starts(i) & tvecRaw+binWidth+binOffset <= stops(i));
                if isempty(tvecRaw)
                    continue;
                end
                
                % remove the offset from the time labels
                tvecLabelCell{i} = tvecRaw;

                % add the end of the last bin onto the histc bins
                tbinsForHistcCell{i} = cat(1, tvecRaw+binOffset, tvecRaw(end) + binOffset + binWidth);
            end
        end
    end
end