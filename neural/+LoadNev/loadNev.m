function [Q info] = loadNev(fname, varargin)
% loads (either from saved Q.mat or by generating one) the Q struct for a given nev
% deferring to handler functions loadNev_ProtocolName to do the real work if necessary

def.overwrite = false; % ignore saved Q.mat and overwrite it?
def.processOnly = false; % don't load into memory, just make sure it's been processed?
def = assignargs(def, varargin);

Q = [];
info = [];

% prompt for file if not provided
if ~exist('fname', 'var') || isempty(fname)
    startDir = getPathToData('monkeyRoot');
    fname = filePicker([], 'filterList', {'*.nev', 'Cerebus Files (*.nev)'}, 'multiple', false, ...
    'prompt', 'Choose NEV File', 'selectMode', 'file', 'extensionFilter', '*.nev', 'startDir', startDir);
end

if ~exist(fname, 'file')
    error('File %s not found\n', fname);
end

fprintf('Processing %s\n', fname);
qfname = getPathToData('Q', 'nev', fname);

if(exist(qfname, 'file'))
    % the Q.mat file exists, we need to check whether the version is acceptable
    % load the info from the file
    ld = load(qfname, 'info');
    
    if ~isfield(ld, 'info')
        % info not found in file
        info = [];
    else
        info = ld.info;
    end 

    % determine if this info is okay
    infoValid = checkInfoNev(info, fname, qfname);

    if ~overwrite && infoValid
        if(~processOnly)
            % exists, loading Q(LFP).mat from disk
            fprintf('\tLoading existing %s\n', qfname);
            ld = load(qfname, 'Q', 'info');
            Q = ld.Q;
            if ~isempty(Q)
                Q = makecol(orderfields(Q));
            end
            info = ld.info;
        else
            % exists but we're not interested in loading it`
            fprintf('\tFound existing %s\n', qfname);
        end

        return;
    else
        % exists, but we're overwriting it
        fprintf('\tOverwriting existing %s\n', qfname);
    end
end

% get the loadNev delegate function 
delegates = getProtocolDelegatesNev();
delegateFn = matchDelegateByFilename(fname, delegates, 'fieldName', 'loadNevFn');
if isempty(delegateFn)
    fprintf(2, '\tError: could not find handler for %s\n', fname);
    return;
end

% Process the nev by calling the delegate function
[Q info] = delegateFn(fname, varargin{:});
assert(isstruct(info), 'Info must be a struct');

% add nev name short field
info.nevNameShort = getShortPath(fname);

if ~isempty(Q)
    Q = makecol(orderfields(Q));
end

% Save to disk!
qfname = getPathToData('Q', 'nev', fname); 
fprintf('\tSaving Q to %s\n', qfname);
saveLarge(qfname, 'Q', 'info');

end

