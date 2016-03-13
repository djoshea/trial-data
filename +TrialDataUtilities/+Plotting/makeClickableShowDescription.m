function makeClickableShowDescription(hvec, varargin)
% makeClickableShowDescription(hvec, [descCell], varargin)
%
% descCell is a cellstr of descriptions to set for each handle in hvec
% parameters 'activateFn' take handles to functions
%   which looks like 
%       handlesToDelete = activateFn(handleCopyVec, intersectionXYZ)
%   activateFn will receive a copy of the selected handles and should
%   modify them to make them standout. Because they are copies, the
%   original handle will remain intact, allowing the copies to be simply deleted
%   when a new handle is selected

    p = inputParser();
    p.addOptional('descCell', {}, @(x) isempty(x) || iscellstr(x));
    p.addParameter('activateFn', @select, @(x) isa(x, 'function_handle'));
    p.parse(varargin{:});
    descCell = p.Results.descCell;
    activateFn = p.Results.activateFn;
    
    if nargin > 1 && ~isempty(descCell)
        for i = 1:numel(hvec)
            set(hvec(i), 'Description', descCell{i});
        end
    end

    set(hvec,  'PickableParts', 'visible', 'HitTest', 'on', 'ButtonDownFcn', @clickFn, 'SelectionHighlight', 'on');

    hSelected = []; % handle to copy of that object displayed
    hText = [];

    % default activate function multiplies line width by 10
    % this will receive a copy of the original object that will later be
    % deleted, so it can modify the properties of that object without
    % worrying about restoring them
    function h = select(h, pt) %#ok<INUSD>
        if ~ishandle(h) || ~isvalid(h), return; end
        isLine = strcmp(get(h, 'Type'), 'line');
        if isLine
            set(h, 'LineWidth',  get(h, 'LineWidth') * 10);
            TrialDataUtilities.Plotting.setLineOpacity(h, 1);
        end
        set(h, 'Selected', 'on');
    end

    function clearSelection(~, ~)
        if ~isempty(hSelected) && isvalid(hSelected)
            delete(hSelected);
            hSelected = [];
        end
        if ~isempty(hText) && isvalid(hText)
            delete(hText);
            hText = [];
        end
    end

    function dispDescription(~, ~)
        % prints the contents of the text box to the command line
        desc = get(hText, 'String');
        fprintf('\n');
        if ischar(desc)
            fprintf('%s\n', desc);
        elseif iscell(desc)
            for iL = 1:numel(desc)
                fprintf('%s\n', desc{iL});
            end
        end
        fprintf('\n');
    end

    function clickFn(h, eventData)
        clearSelection();
        
        % make a copy of the object and pass to activateFn
        hSelected = copy(h);
        hSelected.Parent = h.Parent;
        pt = eventData.IntersectionPoint;
        % receive handles to delete from the function 
        if nargout(activateFn) == 0
            activateFn(hSelected, pt);
        else
            hSelected = activateFn(hSelected, pt);
        end
        
        % make the hSelected copy clickable so we can clear the selection
        set(hSelected,  'PickableParts', 'visible', 'HitTest', 'on', ...
            'ButtonDownFcn', @clearSelection, 'SelectionHighlight', 'off');
        
        % draw the text box annotation with the 'Description'
        if ~isempty(get(h, 'Description'))
            hText = text(pt(1), pt(2), pt(3), get(h, 'Description'), ...
                'Interpreter', 'none', 'VerticalAlignment', 'top', ...
                'Margin', 10, 'BackgroundColor', [1 1 1 0.9]);
        
            % make the hSelected copy clickable so we can clear the selection
            set(hText,  'PickableParts', 'visible', 'HitTest', 'on', ...
                'ButtonDownFcn', @dispDescription, 'SelectionHighlight', 'off');
        end
    end
end