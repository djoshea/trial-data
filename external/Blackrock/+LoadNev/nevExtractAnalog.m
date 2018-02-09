function nsxData = nevExtractAnalog(fname, varargin)
% nsxData = nevExtractAnalog(fname)
%   Loads analog data from all nsx files associated with fname (.nev)
%
% fname : name of nev file (or without extension)
% nsxData : struct array, one element for each located ns# file with the same path/name as fname
%   .fname
%   .info : header info as returned from NSX_open
%   .ext : the ns# extension
%   .samplingHz : sampling frequency converted from nsx headers
%   .data : all channels 
%   .time : in ms
%   .timeStart, .timeSamplePeriod: starting point and inter-sample interval in ms
%   .scaleFn : scaleFn(data) returns a double that has been converted into the correct units
%   .scaleLims : [ digitalMin digitalMax analogMin analogMax ] limits used to build scaleFn

par.nsxExts = {'.ns1', '.ns2', '.ns3', '.ns4', '.ns5', '.ns6'};
par.rescale = true; % convert to actual units indicated, uses more memory (double instead of int16)
par.channelIds = []; % select channels by Channel_ID to keep, throw away the rest
assignargs(par, varargin);

[path name ext] = fileparts(fname);

if isempty(path)
    error('Please provide absolute, not relative, path to file');
end

nsxPath = path;


% find existing .nsX files (which hold continuous data)
% attempt each nsx extension to find analog data
nsxCount = 0;

for iext = 1:length(nsxExts)
    fnameSearch = fullfile(nsxPath, [name nsxExts{iext}]);
    if exist(fnameSearch, 'file')
        % an nsx file with this extension exists, load it up
        nsxCount = nsxCount + 1;

        fprintf('\tLoading analog from %s\n', fnameSearch);
        
        nsxInfo = openNSx(fnameSearch);
        nsxData(nsxCount).fname = fnameSearch;
        nsxData(nsxCount).info = nsxInfo; 
        nsxData(nsxCount).ext = nsxExts{iext};
        nsxData(nsxCount).samplingHz = nsxInfo.MetaTags.SamplingFreq;
%         [data time timeStart timeSamplePeriod] = LoadNev.NSX_read(nsxInfo);
        
        if isempty(channelIds)
            % keep all channels
            nsxData(nsxCount).data = nsxInfo.Data;
            nsxData(nsxCount).channelIds = cat(1, nsxInfo.ElectrodesInfo.ElectrodeID);
            channelIndLookup = 1:size(nsxData(nsxCount).data ,1);
            found = true(numel(channelIndLookup), 1);
        else
            channelIdsInFile = cat(1, nsxInfo.ElectrodesInfo.ElectrodeID);
            [found, channelIndLookup] = ismember(channelIds, channelIdsInFile);
            
            if iscell(nsxInfo.Data)
                nsxInfo.Data = cellfun(@(d) int16(d(channelIndLookup(found), :)), nsxInfo.Data, 'UniformOutput', false);
                nsxData(nsxCount).data = cat(2, nsxInfo.Data{:});
            else
                nsxData(nsxCount).data = int16(nsxInfo.Data(channelIndLookup(found), :));
            end
            nsxData(nsxCount).channelIds = channelIds(found);
        end 
        
        % convert cerebus timestamp to ms
        nsxData(nsxCount).timeStart = nsxInfo.MetaTags.Timestamp(1) / 30;
        
        t = cell(numel(nsxInfo.MetaTags.Timestamp), 1);
        for iPart = 1:numel(nsxInfo.MetaTags.Timestamp)
            % NOTE single does not have sufficient resolution on the
            % integers for this. Convert timestamp offset to ms, then add a
            % clock with the right sampling rate for ms
            t{iPart} = double(nsxInfo.MetaTags.Timestamp(iPart)) / 30 + double(0:nsxInfo.MetaTags.DataPoints(iPart)-1)' * double((1000 / nsxInfo.MetaTags.SamplingFreq));
        end
        nsxData(nsxCount).time = cat(1, t{:}); %#ok<*AGROW>
        nsxData(nsxCount).timeSamplePeriod = 1000 / nsxInfo.MetaTags.SamplingFreq;
        % this function rescales to be in the correct units, but uses more memory
        [nsxData(nsxCount).scaleFns nsxData(nsxCount).scaleLims] = getScaleFns(nsxInfo.ElectrodesInfo, channelIndLookup(found));
        
        trim = @(s) strtok(s, char(0));
        nsxData(nsxCount).chLabels = cellfun(trim, {nsxInfo.ElectrodesInfo(channelIndLookup(found)).Label}', 'UniformOutput', false);
        nsxData(nsxCount).chUnits = cellfun(trim, {nsxInfo.ElectrodesInfo(channelIndLookup(found)).AnalogUnits}', 'UniformOutput', false);
    end
end

if nsxCount == 0
    nsxData = [];
end

end

function [scaleFns scaleLims] = getScaleFns(electrodes, channelInds)
    [scaleFns scaleLims] = deal(cell(length(channelInds), 1));
    for ich = 1:length(channelInds)
        ind = channelInds(ich);
        origLims = double([electrodes(ind).MinDigiValue electrodes(ind).MaxDigiValue]);
        newLims  = double([electrodes(ind).MinAnalogValue electrodes(ind).MaxAnalogValue]);
        scaleLims{ich} = [origLims newLims];
        scaleFns{ich} = @(signal) (double(signal) - origLims(1)) / ...
                                   diff(origLims) * diff(newLims) + newLims(1);
    end
end
