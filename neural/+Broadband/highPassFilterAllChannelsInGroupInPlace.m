function tdca = highPassFilterAllChannelsInGroupInPlace(tdca, groupName, varargin)

    p = inputParser();
    p.addParameter('hpCornerHz', 250, @isscalar);
    p.parse(varargin{:});

    Fs = 30000;
    hpCornerHz = p.Results.hpCornerHz; 
    hpCornerNormalized = hpCornerHz / (Fs/2);
    [B, A] = butter(4, hpCornerNormalized, 'high');

    tdca = tdca.reset;
    tdca = tdca.transformAnalogChannelGroupInPlace(groupName, @transformFn, 'applyScaling', false);

    function data = transformFn(data, tvec, ind) %#ok<INUSD>
        data = TrialDataUtilities.Data.filterIgnoreLeadingTrailingNaNs(B, A, data, ...
            'subtractFirstSample', true, 'filtfilt', false, 'showProgress', false);
    end
    
end

