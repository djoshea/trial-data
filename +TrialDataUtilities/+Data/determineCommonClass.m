function newClass = determineCommonClass(varargin)
    uclasses = unique(varargin);
    uclasses = setdiff(uclasses, '');

    if isempty(uclasses)
        newClass = '';
        return;
    end

    if ismember('cell', uclasses)
        newClass = 'cell';
    elseif ismember('double', uclasses)
        newClass = 'double';
    elseif ismember('single', uclasses)
        newClass = 'single';
    else
        % parse int types
        has_logical = ismember('logical', uclasses);
        signed_types = {'int8', 'int16', 'int32', 'int64'};
        unsigned_types = {'uint8', 'uint16', 'uint32', 'uint64'};
        which_signed = ismember(signed_types, uclasses);
        which_unsigned = ismember(unsigned_types, uclasses);

        if any(which_signed) || any(which_unsigned)
            max_signed = find(which_signed, 1, 'last');
            max_unsigned = find(which_unsigned, 1, 'last');

            if any(which_signed)
                if any(which_unsigned)
                    % return sufficiently large signed type that can
                    % handle unsigned type
                    bitidx = max(max_signed, max_unsigned + 1);
                    newClass = signed_types{bitidx};
                else
                    newClass = signed_types{max_signed};
                end

            elseif any(which_unsigned)
                newClass = unsigned_types{max_unsigned};

            end

        elseif has_logical
            newClass = 'logical';

        else
            error('Should not get here');
        end
    end
end