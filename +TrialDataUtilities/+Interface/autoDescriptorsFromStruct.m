function [out, skipped] = autoDescriptorsFromStruct(S, varargin)

    p = inputParser();
    p.addParameter('timeFields', {}, @iscellstr);
    p.parse(varargin{:});

    flds = fieldnames(S);
    out = '';

    [matchMatrix, szOverTrials] = getLengthMatchesAllTrials(S, flds);
    
    isTimeField = ismember(flds, p.Results.timeFields);
    
    allScalar = all(szOverTrials == 1, 1);
    
    skipMask = falsevec(numel(flds));
    for iF = 1:numel(flds)
        fld = flds{iF};
        val = getFirstNonEmptyValue(S, fld);
        
        if ischar(val)
            % strings
            appendDesc('ParamChannelDescriptor.buildStringParam(''%s'')', fld);
            
        elseif allScalar(iF)
            % params and events
            if contains(lower(fld), 'time')
                % treat as event
                appendDesc('EventChannelDescriptor.buildSingleEvent(''%s'', ''ms'')', fld); 
            else
                appendDesc('ParamChannelDescriptor.buildScalarParam(''%s'', '''')', fld);
            end
            
        elseif isTimeField(iF)
            % skip 
            continue;
            
        elseif isvector(getFirstNonEmptyValue(S, fld))
            % try analog channel
            % find first timeField that matches
            timeFieldIndex = find(matchMatrix(iF, :) & makerow(isTimeField), 1);
            if isempty(timeFieldIndex)
                timeField = '???timeField???';
            else
                timeField = flds{timeFieldIndex};
            end
            appendDesc('AnalogChannelDescriptor.buildVectorAnalog(''%s'', ''%s'', '''');', fld, timeField);
        else
            skipMask(iF) = true;
        end
    end
    
    skipped = flds(skipMask);
    
    fprintf(out);
    
    function appendLn(varargin)
        out = [out, sprintf(varargin{:}), newline];
    end

    function appendDesc(varargin)
        appendLn('channelDescriptors(end+1) = %s;', sprintf(varargin{:}));
    end

end

function [matchMatrix, szMat] = getLengthMatchesAllTrials(S, flds)
    % nF x nF all sizes of values equivalent matrix
    nF = numel(flds);
    szMat = nan(numel(S), nF);
    for iF = 1:nF
        fld = flds{iF};
        szMat(:, iF) = arrayfun(@(s) size(makecol(s.(fld)), 1), S);
    end
    
    matchMatrix = false(nF, nF);    
    for iF = 1:nF
        for iG = 1:nF
            if iF == iG
                continue;
            end
            matchMatrix(iF, iG) = all(szMat(:, iF) == szMat(:, iG));
        end
    end
end

function val = getFirstNonEmptyValue(S, fld)
    for iS = 1:numel(S)
        if ~isempty(S(iS).(fld))
            val = S(iS).(fld);
            return;
        end
    end
    val = [];
end