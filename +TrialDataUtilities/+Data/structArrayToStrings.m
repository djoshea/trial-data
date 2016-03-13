function strCell = structArrayToStrings(s, varargin)
% see structToString documentation for arguments
    strCell = arrayfun(@(el) TrialDataUtilities.Data.structToString(el, varargin{:}), s, 'UniformOutput', false);
end